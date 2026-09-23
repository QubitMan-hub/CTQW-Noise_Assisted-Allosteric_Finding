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
from scipy.stats import rankdata

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
        return _unitary_populations(H.toarray() if sparse.issparse(H) else H, sources, tlist)
    n = H.shape[0]
    Hs = sparse.csr_matrix(H, dtype=complex)
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


def available_memory():
    """Bytes of memory free for new work, or None if the platform will not say."""
    try:
        with open("/proc/meminfo") as f:                      # Linux
            for line in f:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) * 1024
    except OSError:
        pass
    if sys.platform == "win32":
        import ctypes

        class _Mem(ctypes.Structure):
            _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong)] + [
                (k, ctypes.c_ulonglong) for k in ("total", "avail", "total_pf", "avail_pf",
                                                  "total_virt", "avail_virt", "avail_ext")]
        m = _Mem()
        m.dwLength = ctypes.sizeof(_Mem)
        if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m)):
            return int(m.avail)
    return None


def walk_memory(n):
    """Peak memory of one noisy walk on n residues: the RK45 state and temporaries hold
    about 20 n x n complex arrays (measured: ~98 MB at n=566, ~330 MB at n=1058),
    plus a Python/numpy baseline."""
    return 20 * 16 * n * n + 150 * 2 ** 20


def safe_workers(n, requested=None):
    """CPU cores to use, lowered so that the workers fit in the free memory."""
    w = default_workers() if requested is None else max(1, int(requested))
    avail = available_memory()
    if avail:
        w = min(w, max(1, int(0.7 * avail // walk_memory(n))))
    return w


def sweep_many(hams, jobs, tlist, workers=None):
    """Run many independent walks and return their visiting-score vectors in order.
    hams: {key: H}; jobs: [(key, source index or indices, gamma), ...].
    gamma = 0 jobs are exact and instant, so they run here; the rest are spread
    over worker processes (same numbers as a serial run, less wall time)."""
    tlist = np.asarray(tlist, dtype=float)
    out = [None] * len(jobs)

    def run_here(i):
        key, src, g = jobs[i]
        out[i] = visiting_scores(run_walk(hams[key], src, g, tlist), tlist)

    slow = [i for i, job in enumerate(jobs) if job[2] != 0]
    for i in set(range(len(jobs))) - set(slow):
        run_here(i)
    n = next(iter(hams.values())).shape[0] if hams else 0
    workers = safe_workers(n, workers)
    if workers == 1 or len(slow) <= 1:
        for i in slow:
            run_here(i)
        return out
    # one BLAS thread per worker, otherwise the processes fight over the cores
    saved = {v: os.environ.get(v) for v in _THREAD_VARS}
    os.environ.update({v: "1" for v in _THREAD_VARS})
    try:
        ctx = multiprocessing.get_context("spawn")
        # workers get sparse copies of H: a few contacts per residue instead of n x n
        light = {k: sparse.csr_matrix(np.asarray(H, dtype=complex)) for k, H in hams.items()}
        with ProcessPoolExecutor(max_workers=min(workers, len(slow)), mp_context=ctx,
                                 initializer=_pool_init, initargs=(light, tlist)) as ex:
            futures = {ex.submit(_pool_job, jobs[i]): i for i in slow}
            for f in as_completed(futures):
                out[futures[f]] = f.result()
    except BrokenProcessPool:
        # e.g. an interactive session whose main module cannot be re-imported
        print("[note] worker processes unavailable; running the sweep serially.", file=sys.stderr)
        for i in slow:
            if out[i] is None:
                run_here(i)
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
    pos = np.asarray(positive_mask, dtype=bool)
    n_pos, n_neg = int(pos.sum()), int((~pos).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    ranks = rankdata(np.asarray(scores, dtype=float))          # ties get their average rank
    return float((ranks[pos].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def distal_mask(G, nodes, source_ids, fraction=0.5):
    """Residues in the far part of the network from the source(s): hop distance at
    least ceil(fraction * eccentricity). Averaging transport over this set is far
    less arbitrary than a single 'farthest' residue. Returns (mask, min_hops)."""
    dist = hop_distance(G, nodes, [source_ids] if isinstance(source_ids, str) else source_ids)
    min_hops = max(1, int(np.ceil(fraction * max(dist.max(), 0))))
    return dist >= min_hops, min_hops


def hop_distance(G, nodes, source_ids):
    """Fewest contacts from any source to each residue (-1 if unreachable)."""
    dist = {}
    for s in source_ids:
        for k, d in nx.single_source_shortest_path_length(G, s).items():
            dist[k] = min(d, dist.get(k, d))
    return np.array([dist.get(n, -1) for n in nodes], dtype=int)


def distance_shells(dist, eligible, min_size=5):
    """Group eligible residues by hop distance, merging neighbouring distances
    until every group has at least min_size residues. Returns a list of index arrays."""
    groups, cur = [], []
    for d in sorted(set(dist[eligible].tolist())):
        cur.extend(np.where(eligible & (dist == d))[0].tolist())
        if len(cur) >= min_size:
            groups.append(cur)
            cur = []
    if cur:
        if groups:
            groups[-1].extend(cur)
        else:
            groups.append(cur)
    return [np.array(g) for g in groups]


def shell_percentile(scores, shells, n):
    """Each residue's score as a percentile among residues at the same distance
    from the source (ties averaged; 0..1; nan outside the shells). It keeps only
    what the score says beyond 'how close is this residue to the start'."""
    scores, out = np.asarray(scores, dtype=float), np.full(n, np.nan)
    for g in shells:
        out[g] = (rankdata(scores[g]) - 0.5) / len(g)          # ties get their average rank
    return out


def classical_scores(A, source_idx, rates, tmax):
    """Classical baseline: a continuous-time random walk on the same contacts,
    dp/dt = -k L p with L = D - A, started from the same (mixed) sources. Returns
    the time-integrated occupation up to tmax for each hopping rate k, exactly,
    from one eigendecomposition of L (so every rate costs almost nothing).
    Site energies do not enter: a classical walker has no phases to shift."""
    A = A.toarray() if sparse.issparse(A) else np.asarray(A, dtype=float)
    L = np.diag(A.sum(axis=1)) - A
    w, V = np.linalg.eigh(L)
    w = np.clip(w, 0.0, None)
    src = [source_idx] if np.isscalar(source_idx) else list(source_idx)
    p0 = np.zeros(len(A))
    p0[src] = 1.0 / len(src)
    c = V.T @ p0
    out = []
    for k in rates:
        x = k * w * tmax
        f = np.where(x > 1e-12, -np.expm1(-x) / np.where(w > 0, k * w, 1.0), tmax)
        out.append(V @ (c * f))
    return np.array(out)


def shell_permutation_test(score_rows, shells, positive, n_perm=10000, seed=0):
    """Significance of the distance-adjusted AUC. The known allosteric labels are
    shuffled within each distance shell (so every shuffle keeps the same number
    of allosteric residues at each distance), and the best AUC over all rows
    (noise levels, or classical rates) is recomputed each time. Taking the best
    over rows in the null too makes p honest about having picked the best gamma."""
    idx = np.concatenate(shells)
    pos = np.asarray(positive, dtype=bool)
    n_pos = int(pos[idx].sum())
    n_neg = len(idx) - n_pos
    if n_pos == 0 or n_neg == 0:
        return None
    R = np.array([rankdata(shell_percentile(row, shells, len(pos))[idx]) for row in score_rows], dtype=np.float32)
    offset = n_pos * (n_pos + 1) / 2.0
    observed = (R[:, pos[idx]].sum(axis=1) - offset) / (n_pos * n_neg)
    rng = np.random.default_rng(seed)
    P = np.zeros((len(idx), n_perm), dtype=np.float32)
    start = 0
    for g in shells:
        m, k = len(g), int(pos[g].sum())
        if k:
            pick = np.argsort(rng.random((n_perm, m)), axis=1)[:, :k]
            P[start + pick, np.arange(n_perm)[:, None]] = 1.0
        start += m
    null_best = ((R @ P - offset) / (n_pos * n_neg)).max(axis=0)
    best = float(np.nanmax(observed))
    return {"observed_best": best, "p_value": float((1 + np.sum(null_best >= best - 1e-9)) / (n_perm + 1)),
            "null_mean": float(null_best.mean()), "null_95": float(np.quantile(null_best, 0.95)),
            "n_perm": int(n_perm)}


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
