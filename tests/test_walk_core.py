#!/usr/bin/env python3
"""Checks that the fast solver reproduces the original dense solver. Run: python tests/test_walk_core.py (or pytest)."""
import os, sys
import numpy as np
import networkx as nx
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import walk_core as wc

TLIST = np.linspace(0, 30, 200)


def reference_walk(H, source_idx, gamma, tlist):
    """The original dense solver (walk_core before the speed-ups), kept verbatim."""
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
    return np.array([np.real(np.diag(np.ascontiguousarray(sol.y[:, i]).view(complex).reshape(n, n)))
                     for i in range(len(tlist))]).T


def protein_like(n=40, seed=3):
    G = nx.random_geometric_graph(n, 0.3, seed=seed)
    G = G.subgraph(max(nx.connected_components(G), key=len)).copy()
    A = nx.to_numpy_array(G)
    names = list(wc.KYTE_DOOLITTLE)
    resnames = [names[i % len(names)] for i in range(len(A))]
    return A, resnames


def test_matches_reference_solver():
    A, rn = protein_like()
    H = wc.build_hamiltonian(A, rn, "hydropathy", 3.0)
    for g in (0.0, 0.1, 1.5, 50.0):
        new = wc.visiting_scores(wc.run_walk(H, 0, g, TLIST), TLIST)
        ref = wc.visiting_scores(reference_walk(H, 0, g, TLIST), TLIST)
        # gamma = 0 is now exact; the reference integrator itself carries ~1e-6 error there
        assert np.max(np.abs(new - ref)) < (1e-5 if g == 0 else 1e-9), g


def test_mixed_sources_equal_average_of_single_sources():
    A, rn = protein_like()
    H = wc.build_hamiltonian(A, rn, "random", 3.0, seed=1)
    srcs = [0, 5, 9]
    for g in (0.0, 0.8):
        mix = wc.run_walk(H, srcs, g, TLIST)
        avg = np.mean([wc.run_walk(H, s, g, TLIST) for s in srcs], axis=0)
        assert np.max(np.abs(mix - avg)) < 1e-6, g


def test_parallel_equals_serial():
    A, rn = protein_like()
    H = wc.build_hamiltonian(A, rn)
    gammas = [0.0, 0.2, 2.0, 20.0]
    serial = wc.sweep_scores(H, 0, gammas, TLIST, workers=1)
    parallel = wc.sweep_scores(H, 0, gammas, TLIST, workers=3)
    assert np.array_equal(serial, parallel)


def test_populations_conserved():
    A, rn = protein_like()
    H = wc.build_hamiltonian(A, rn)
    for g in (0.0, 3.0):
        pop = wc.run_walk(H, [1, 2], g, TLIST)
        assert np.max(np.abs(pop.sum(axis=0) - 1.0)) < 1e-6


def test_hump_stats():
    s = wc.hump_stats([0, 0.1, 1, 10], [0.2, 0.5, 0.4, 0.1])
    assert s["verdict"] == "hump" and s["peak_gamma"] == 0.1 and abs(s["gain_abs"] - 0.3) < 1e-12
    assert wc.hump_stats([0, 1, 10], [0.5, 0.4, 0.1])["verdict"] == "quantum"


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print(f"ok  {name}")
    print("all checks passed")
