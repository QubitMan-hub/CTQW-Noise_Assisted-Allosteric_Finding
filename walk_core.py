#!/usr/bin/env python3
"""walk_core.py - shared physics for the environment-assisted quantum walk.

A continuous-time quantum walk runs on a protein's residue network; a dephasing
rate gamma dials it from fully quantum (0) to strongly dephased (large). Every
other script is a thin layer over the functions here, so they all agree.

Numerics (exact reformulations of the same master equation, checked against the
original dense solver in tests/test_walk_core.py):
  * gamma = 0 is solved exactly by diagonalising H (no time stepping);
  * gamma > 0 uses a sparse H and -i[H, rho] = -i(M - M^dagger) with M = H rho;
  * several sources are one mixed initial state (the equation is linear);
  * only populations are stored, so memory grows as N^2, not N^2 x time points.
"""
import multiprocessing
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool

import numpy as np
import networkx as nx
from scipy import sparse
from scipy.integrate import RK45

_trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
_THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")

# Kyte-Doolittle hydropathy, one value per amino acid. Used as residue site
# energies: hydrophobic vs polar residues sit in different local environments,
# so they get different on-site energies. Deterministic and defensible, unlike
# arbitrary random numbers.
KYTE_DOOLITTLE = {
    "ALA": 1.8, "ARG": -4.5, "ASN": -3.5, "ASP": -3.5, "CYS": 2.5, "GLN": -3.5,
    "GLU": -3.5, "GLY": -0.4, "HIS": -3.2, "ILE": 4.5, "LEU": 3.8, "LYS": -3.9,
    "MET": 1.9, "PHE": 2.8, "PRO": -1.6, "SER": -0.8, "THR": -0.7, "TRP": -0.9,
    "TYR": -1.3, "VAL": 4.2,
}

# Default dephasing grid: gamma = 0 plus 4 rates per decade from 0.01 to 100.
DEFAULT_GAMMAS = [0.0] + [float(f"{g:.4g}") for g in np.logspace(-2, 2, 17)]


def load_network(graphml_path):
    """Return (G, node_ids, adjacency, residue_names) from a GraphML network."""
    G = nx.read_graphml(graphml_path)
    nodes = list(G.nodes())
    A = nx.to_numpy_array(G, nodelist=nodes, weight="weight")
    A = np.maximum(A, A.T)
    np.fill_diagonal(A, 0.0)
    resnames = [G.nodes[n].get("resname") for n in nodes]
    return G, nodes, A, resnames


def _zscore(v):
    """Centre and unit-scale, so 'scale' means the same disorder for any pattern."""
    v = np.asarray(v, dtype=float)
    return (v - v.mean()) / v.std() if v.std() > 1e-12 else np.zeros_like(v)


def site_energy_pattern(mode, resnames, degrees, seed):
    """Unit-spread site energies. Multiply by a scale to set disorder strength.
    none: flat. hydropathy: residue hydrophobicity (default, deterministic).
    degree: by connectivity. random: null-model control."""
    n = len(degrees)
    if mode == "none":
        return np.zeros(n)
    if mode == "hydropathy":
        return _zscore([KYTE_DOOLITTLE.get((r or "").upper(), 0.0) for r in resnames])
    if mode == "degree":
        return _zscore(degrees)
    if mode == "random":
        return _zscore(np.random.default_rng(seed).standard_normal(n))
    raise ValueError(f"unknown site-energy mode: {mode}")


def build_hamiltonian(A, resnames, mode="hydropathy", scale=3.0, seed=0):
    """H = contacts (off-diagonal) + site energies (diagonal)."""
    A = np.asarray(A, dtype=float)
    pattern = site_energy_pattern(mode, resnames, A.sum(axis=1), seed)
    H = A.astype(complex)
    H[np.diag_indices_from(H)] = scale * pattern
    return H


# ------------------------------------------------------------------ dynamics
def _unitary_populations(H, sources, tlist):
    """gamma = 0: exact populations from the eigendecomposition of H."""
    w, V = np.linalg.eigh(np.asarray(H, dtype=complex))
    phases = np.exp(-1j * np.outer(w, tlist))                  # (n, T)
    pop = np.zeros((H.shape[0], len(tlist)))
    for s in sources:
        psi = V @ (V[s].conj()[:, None] * phases)              # (n, T)
        pop += np.abs(psi) ** 2
    return pop / len(sources)


def run_walk(H, source_idx, gamma, tlist, rtol=1e-6, atol=1e-9):
    """Dephased walk -> populations, shape (N, len(tlist)).
    d rho/dt = -i[H, rho] - gamma * (off-diagonal part of rho).
    source_idx is one residue index, or several (an equal mixture of them)."""
    sources = [int(s) for s in np.atleast_1d(source_idx)]
    tlist = np.asarray(tlist, dtype=float)
    if gamma == 0:
        return _unitary_populations(H, sources, tlist)
    n = H.shape[0]
    Hs = sparse.csr_matrix(np.asarray(H, dtype=complex))
    diag = np.arange(n)

    def rhs(_t, y):
        rho = y.view(complex).reshape(n, n)
        M = Hs @ rho
        drho = -1j * (M - M.conj().T)                         # -i[H, rho], rho Hermitian
        drho -= gamma * rho
        drho[diag, diag] += gamma * rho[diag, diag]            # only coherences decay
        return drho.reshape(-1).view(float)

    rho0 = np.zeros((n, n), dtype=complex)
    rho0[sources, sources] = 1.0 / len(sources)
    # Same stepping and interpolation as solve_ivp(..., t_eval=tlist), but only the
    # populations are kept: memory O(N^2) instead of O(N^2 * len(tlist)).
    solver = RK45(rhs, tlist[0], rho0.reshape(-1).view(float).copy(), tlist[-1], rtol=rtol, atol=atol)
    pop = np.empty((n, len(tlist)))
    flat_diag = 2 * (diag * n + diag)                          # real parts of rho[i, i] in the float view
    i = 0
    while solver.status == "running":
        solver.step()
        if solver.status == "failed":
            raise RuntimeError(f"walk integration failed at gamma={gamma}")
        j = int(np.searchsorted(tlist, solver.t, side="right"))
        if j > i:
            pop[:, i:j] = solver.dense_output()(tlist[i:j])[flat_diag]
            i = j
    return pop


def populations_from_sources(H, source_indices, gamma, tlist):
    """Walk from an equal mixture of several sources (one run, by linearity)."""
    return run_walk(H, list(source_indices), gamma, tlist)


def visiting_scores(pop_time, tlist):
    """Signal each residue receives, integrated over time (higher = more reached)."""
    return _trapz(pop_time, tlist, axis=1)


_POOL_STATE = {}


def _pool_init(hams, tlist):
    _POOL_STATE.update(hams=hams, tlist=tlist)


def _pool_job(job):
    key, sources, gamma = job
    tl = _POOL_STATE["tlist"]
    return visiting_scores(run_walk(_POOL_STATE["hams"][key], sources, gamma, tl), tl)


def default_workers():
    return max(1, min(8, os.cpu_count() or 1))


def sweep_many(hams, jobs, tlist, workers=None):
    """Run many independent walks and return their visiting-score vectors in order.
    hams: {key: H}; jobs: [(key, source index or indices, gamma), ...].
    gamma = 0 jobs are exact and instant, so they run here; the rest are spread
    over worker processes (same numbers as a serial run, less wall time)."""
    tlist = np.asarray(tlist, dtype=float)
    out = [None] * len(jobs)
    slow = []
    for i, (key, src, g) in enumerate(jobs):
        if g == 0:
            out[i] = visiting_scores(run_walk(hams[key], src, 0.0, tlist), tlist)
        else:
            slow.append(i)
    workers = default_workers() if workers is None else max(1, int(workers))
    if workers == 1 or len(slow) <= 1:
        for i in slow:
            key, src, g = jobs[i]
            out[i] = visiting_scores(run_walk(hams[key], src, g, tlist), tlist)
        return out
    # one BLAS thread per worker, otherwise the processes fight over the cores
    saved = {v: os.environ.get(v) for v in _THREAD_VARS}
    os.environ.update({v: "1" for v in _THREAD_VARS})
    try:
        ctx = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(max_workers=min(workers, len(slow)), mp_context=ctx,
                                 initializer=_pool_init, initargs=(hams, tlist)) as ex:
            futures = {ex.submit(_pool_job, jobs[i]): i for i in slow}
            for f in as_completed(futures):
                out[futures[f]] = f.result()
    except BrokenProcessPool:
        # e.g. an interactive session whose main module cannot be re-imported
        print("[note] worker processes unavailable; running the sweep serially.", file=sys.stderr)
        for i in slow:
            if out[i] is None:
                key, src, g = jobs[i]
                out[i] = visiting_scores(run_walk(hams[key], src, g, tlist), tlist)
    finally:
        for v, old in saved.items():
            if old is None:
                os.environ.pop(v, None)
            else:
                os.environ[v] = old
    return out


def sweep_scores(H, source_idx, gammas, tlist, workers=None):
    """Visiting scores of every residue at every gamma, shape (len(gammas), N)."""
    return np.array(sweep_many({0: H}, [(0, source_idx, g) for g in gammas], tlist, workers))


# ------------------------------------------------------------------ metrics
def roc_auc(scores, positive_mask):
    """ROC AUC via Mann-Whitney with tie handling. 0.5 random, 1.0 perfect, nan
    if no positives/negatives."""
    scores = np.asarray(scores, dtype=float)
    pos = np.asarray(positive_mask, dtype=bool)
    n_pos, n_neg = int(pos.sum()), int((~pos).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(len(scores))
    ranks[order] = np.arange(1, len(scores) + 1)
    s = scores[order]
    i = 0
    while i < len(s):                      # average ranks within ties
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (ranks[order[i]] + ranks[order[j]]) / 2.0
        i = j + 1
    return float((ranks[pos].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def sweep_auc(H, source_indices, positive_mask, eligible_mask, gammas, tlist, workers=None):
    """AUC for recovering the positive residues at each gamma (eligible only)."""
    scores = sweep_scores(H, list(source_indices), gammas, tlist, workers)
    return np.array([roc_auc(s[eligible_mask], positive_mask[eligible_mask]) for s in scores])


def distal_mask(G, nodes, source_ids, fraction=0.5):
    """Residues in the far part of the network from the source(s): hop distance at
    least ceil(fraction * eccentricity). Averaging transport over this set is far
    less arbitrary than a single 'farthest' residue. Returns (mask, min_hops)."""
    source_ids = [source_ids] if isinstance(source_ids, str) else list(source_ids)
    dist = {}
    for s in source_ids:
        for k, d in nx.single_source_shortest_path_length(G, s).items():
            dist[k] = min(d, dist.get(k, d))
    ecc = max(dist.values()) if dist else 0
    min_hops = max(1, int(np.ceil(fraction * ecc)))
    mask = np.array([dist.get(n, -1) >= min_hops for n in nodes], dtype=bool)
    return mask, min_hops


def hump_stats(gammas, values):
    """Peak, verdict and effect size of a curve over the gamma grid.
    baseline = the better of the two ends (fully quantum, most dephased);
    gain = peak - baseline; 'helps' = grid range where the curve beats both ends."""
    g, v = np.asarray(gammas, dtype=float), np.asarray(values, dtype=float)
    peak = int(np.nanargmax(v))
    baseline = float(max(v[0], v[-1]))
    interior = 0 < peak < len(v) - 1 and v[peak] > baseline + 1e-12
    verdict = "hump" if interior else ("quantum" if peak == 0 else "dephased")
    above = np.where(v > baseline + 1e-12)[0]
    helps = [float(g[above[0]]), float(g[above[-1]])] if interior and len(above) else None
    return {"peak_index": peak, "peak_gamma": float(g[peak]), "peak_value": float(v[peak]),
            "quantum_value": float(v[0]), "dephased_value": float(v[-1]), "baseline": baseline,
            "gain_abs": float(v[peak] - baseline),
            "gain_rel": float((v[peak] - baseline) / baseline) if baseline > 0 else float("nan"),
            "helps_range": helps,
            "helps_decades": (float(np.log10(helps[1] / helps[0])) if helps and helps[0] > 0 else None),
            "verdict": verdict}


# ------------------------------------------------------------------ quantum encoding
def pad_to_power_of_two(H):
    """Pad an N x N Hermitian matrix to 2^n x 2^n (amplitude encoding needs it)."""
    n_qubits = max(1, int(np.ceil(np.log2(H.shape[0]))))
    Hp = np.zeros((2 ** n_qubits, 2 ** n_qubits), dtype=complex)
    Hp[:H.shape[0], :H.shape[0]] = H
    return Hp, n_qubits


# per-qubit change of basis (M00, M01, M10, M11) -> Pauli (I, X, Y, Z)
_PAULI_BASIS = 0.5 * np.array([[1, 0, 0, 1], [0, 1, 1, 0],
                               [0, 1j, -1j, 0], [1, 0, 0, -1]], dtype=complex)


def pauli_decompose(H, tol=1e-9):
    """2^n x 2^n Hermitian matrix -> [(pauli_string, coeff), ...] (qubit 0 first)."""
    n = int(round(np.log2(H.shape[0])))
    T = H.reshape([2] * (2 * n))
    perm = [x for q in range(n) for x in (q, n + q)]
    T = T.transpose(perm).reshape([4] * n)
    for q in range(n):
        T = np.moveaxis(np.tensordot(_PAULI_BASIS, T, axes=([1], [q])), 0, q)
    labels = np.array(["I", "X", "Y", "Z"])
    return ([("".join(labels[list(idx)]), float(np.real(T[idx])))
             for idx in np.ndindex(*([4] * n)) if abs(T[idx]) > tol], n)


def pauli_reconstruct(terms, n):
    """Rebuild the matrix from Pauli terms, for self-checking."""
    P = {"I": np.eye(2, dtype=complex), "X": np.array([[0, 1], [1, 0]], dtype=complex),
         "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
         "Z": np.array([[1, 0], [0, -1]], dtype=complex)}
    M = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for s, c in terms:
        op = np.array([[1]], dtype=complex)
        for ch in s:
            op = np.kron(op, P[ch])
        M += c * op
    return M


def _pauli_actions(terms, n):
    """For each Pauli string: (index map, phase vector) with (P psi) = phase * psi[map].
    Character 0 of a string is the most significant bit, as in pauli_decompose."""
    k = np.arange(2 ** n)
    parity = np.zeros(2 ** n, dtype=np.int8)
    for b in range(n):
        parity ^= ((k >> b) & 1).astype(np.int8)
    acts = []
    for s, c in terms:
        xm = zm = ny = 0
        for p, ch in enumerate(s):
            bit = 1 << (n - 1 - p)
            if ch in "XY":
                xm |= bit
            if ch in "ZY":
                zm |= bit
            ny += ch == "Y"
        src = k ^ xm                                           # P|j> lands on |j ^ xm>
        phase = (1j ** ny) * (1 - 2 * parity[src & zm].astype(float))
        acts.append((src, phase, c))
    return acts


def _suzuki(psi, acts, dt, order):
    """One Suzuki-Trotter step of the given (1 or even) order, term order as listed."""
    if order == 1:
        for src, phase, c in acts:
            psi = np.cos(c * dt) * psi - 1j * np.sin(c * dt) * phase * psi[src]
        return psi
    if order == 2:
        for src, phase, c in acts:
            psi = np.cos(c * dt / 2) * psi - 1j * np.sin(c * dt / 2) * phase * psi[src]
        for src, phase, c in reversed(acts):
            psi = np.cos(c * dt / 2) * psi - 1j * np.sin(c * dt / 2) * phase * psi[src]
        return psi
    p = 1.0 / (4.0 - 4.0 ** (1.0 / (order - 1)))
    for f in (p, p, 1 - 4 * p, p, p):
        psi = _suzuki(psi, acts, f * dt, order - 2)
    return psi


def trotter_check(terms, n_qubits, H_padded, source_index, time, order, repetitions, n_residues):
    """How far the product-formula circuit (standard Suzuki-Trotter, terms in the
    listed order) is from exact exp(-iHt) on the source state. Returns fidelity
    |<exact|trotter>|^2 and the largest residue-population error."""
    if order != 1 and order % 2:
        raise ValueError("Suzuki-Trotter order must be 1 or even.")
    psi0 = np.zeros(2 ** n_qubits, dtype=complex)
    psi0[source_index] = 1.0
    w, V = np.linalg.eigh(H_padded)
    exact = V @ (np.exp(-1j * w * time) * V[source_index].conj())
    acts = _pauli_actions(terms, n_qubits)
    psi = psi0
    for _ in range(repetitions):
        psi = _suzuki(psi, acts, time / repetitions, order)
    fid = float(abs(np.vdot(exact, psi)) ** 2)
    dpop = float(np.max(np.abs(np.abs(exact[:n_residues]) ** 2 - np.abs(psi[:n_residues]) ** 2)))
    return {"fidelity": fid, "max_population_error": dpop}
