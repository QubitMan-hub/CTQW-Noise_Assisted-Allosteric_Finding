#!/usr/bin/env python3
"""walk_core.py - shared physics for the environment-assisted quantum walk.

A continuous-time quantum walk runs on a protein's residue network; a dephasing
rate gamma dials it from fully quantum (0) to fully classical (large). Every
other script is a thin CLI over the functions here, so they all agree.
"""
import numpy as np
import networkx as nx
from scipy.integrate import solve_ivp

_trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

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


def run_walk(H, source_idx, gamma, tlist):
    """Dephased walk from one source residue -> populations, shape (N, len(tlist)).
    d rho/dt = -i[H, rho] - gamma * (off-diagonal part of rho)."""
    n = H.shape[0]

    def rhs(_t, y):
        rho = np.ascontiguousarray(y).view(complex).reshape(n, n)
        drho = -1j * (H @ rho - rho @ H)
        if gamma > 0:
            off = rho.copy()
            np.fill_diagonal(off, 0.0)
            drho -= gamma * off
        return drho.reshape(-1).view(float)

    rho0 = np.zeros((n, n), dtype=complex)
    rho0[source_idx, source_idx] = 1.0
    sol = solve_ivp(rhs, (tlist[0], tlist[-1]), rho0.reshape(-1).view(float).copy(),
                    t_eval=tlist, method="RK45", rtol=1e-6, atol=1e-9)
    pop = np.empty((n, len(tlist)))
    for ti in range(sol.y.shape[1]):
        rho = np.ascontiguousarray(sol.y[:, ti]).view(complex).reshape(n, n)
        pop[:, ti] = np.real(np.diag(rho))
    return pop


def populations_from_sources(H, source_indices, gamma, tlist):
    """Average the walk over several sources (= starting from an equal mixture,
    since the master equation is linear in the initial state)."""
    total = None
    for s in source_indices:
        pop = run_walk(H, s, gamma, tlist)
        total = pop if total is None else total + pop
    return total / len(source_indices)


def visiting_scores(pop_time, tlist):
    """Signal each residue receives, integrated over time (higher = more reached)."""
    return _trapz(pop_time, tlist, axis=1)


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


def sweep_auc(H, source_indices, positive_mask, eligible_mask, gammas, tlist):
    """AUC for recovering the positive residues at each gamma (eligible only)."""
    out = []
    for g in gammas:
        pop = populations_from_sources(H, source_indices, g, tlist)
        scores = visiting_scores(pop, tlist)
        out.append(roc_auc(scores[eligible_mask], positive_mask[eligible_mask]))
    return np.array(out)


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
