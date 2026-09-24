#!/usr/bin/env python3
"""external_baselines.py - established allosteric-site predictors under the same test.

The walks are only compared with each other in the pipeline. This script scores the
same proteins with three established, published predictors and puts them through
exactly the same distance-controlled AUC and within-shell permutation test:

  closeness    closeness centrality on the contact network (Amitai et al. 2004)
  betweenness  betweenness centrality on the contact network (del Sol et al. 2006)
  prs          perturbation-response scanning on an anisotropic network model of the
               C-alpha atoms (ANM, 15 A cutoff), computed with ProDy: the mean response
               of each residue to perturbing the active-site residues (Atilgan &
               Atilgan 2009; General et al. 2014)

Nothing in the pipeline is changed or rerun. For each protein the network, labels and
distance shells are rebuilt with the pipeline's own functions and checked against the
stored run (residue count, allosteric count and the plain classical walk's best AUC
beyond distance must all match), so every method is scored on the same residues.

    pip install prody            # only needed for this script
    python results/external_baselines.py

Writes results/external.csv and .md; make_results.py turns them into LaTeX macros.
"""
import csv, json, os, sys, tempfile
import numpy as np
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import labels as lb            # noqa: E402
import rin_builder as rb       # noqa: E402
import run_protein as rp       # noqa: E402
import walk_core as wc         # noqa: E402
from make_results import PROTEINS, RUNS     # noqa: E402

METHODS = ["closeness", "betweenness", "prs"]
N_PERM, SEED, ANM_CUTOFF = 10000, 0, 15.0


def prs_from_active(ca, active):
    """Mean normalised PRS response of every residue to perturbing the active site (ProDy)."""
    import prody
    prody.confProDy(verbosity="none")
    anm = prody.ANM()
    anm.buildHessian(ca, cutoff=ANM_CUTOFF)
    anm.calcModes(n_modes=None)                 # all non-zero modes
    prs, _, _ = prody.calcPerturbResponse(anm)  # prs[i, j]: response of j to perturbing i
    return prs[active].mean(axis=0)


def one(pid, tmp):
    with open(os.path.join(RUNS, pid, f"{pid}_result.json")) as f:
        stored = json.load(f)
    labels, _, _ = lb.resolve(pid, "auto")
    chains = set(lb.chains_of(labels))
    G = rp.structure_to_graph(pid, tmp, 0, chains, "cb", 8.0, 0, False)
    _, nodes, A, resnames = wc.network_arrays(G)
    N = len(nodes)
    rep = lb.match(labels, nodes, resnames)
    el = np.ones(N, dtype=bool)
    el[rep["active"]] = False
    pos = np.zeros(N, dtype=bool)
    pos[rep["allosteric"]] = True
    pos &= el
    hops = wc.hop_distance(G, nodes, [nodes[i] for i in rep["active"]])
    reach = el & (hops >= 0)
    shells = wc.distance_shells(hops, reach)
    adj = lambda sc: wc.roc_auc(wc.shell_percentile(sc, shells, N)[reach], pos[reach])

    # the same residues and shells as the stored run, or stop
    a = stored["allosteric"]
    cw = wc.classical_walk(A, list(rep["active"]), 30.0)
    cl = max(adj(cw(k)) for k in rp.CLASSICAL_RATES)
    ok = (N == stored["summary"]["residues"] and int(pos.sum()) == a["n_allosteric"]
          and abs(cl - a["adjusted"]["baselines"]["classical walk, best rate"]) < 1e-9)
    if not ok:
        sys.exit(f"{pid}: rebuilt network does not match the stored run; stopping.")

    # C-alpha coordinates for the ANM (the network node's C-beta if a residue has no C-alpha)
    ca_of = {rb.node_label(r): r["residue"]["CA"].get_coord() for r in
             rb.collect_residues(rb.load_structure(pid, tmpdir=tmp), 0, chains) if "CA" in r["residue"]}
    ca = np.array([ca_of.get(n, [G.nodes[n]["x"], G.nodes[n]["y"], G.nodes[n]["z"]]) for n in nodes], dtype=float)

    close, betw = nx.closeness_centrality(G), nx.betweenness_centrality(G)
    scores = {"closeness": np.array([close[n] for n in nodes]),
              "betweenness": np.array([betw[n] for n in nodes]),
              "prs": prs_from_active(ca, list(rep["active"]))}
    row = {"pdb": pid, "residues": N, "no_ca": sum(n not in ca_of for n in nodes)}
    for m in METHODS:
        s = scores[m]
        sig = wc.shell_permutation_test([s], shells, pos & reach, n_perm=N_PERM, seed=SEED)
        row[f"{m}_auc"] = wc.roc_auc(s[el], pos[el])
        row[f"{m}_beyond"] = adj(s)
        row[f"{m}_p"] = sig["p_value"]
    return row


def main():
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        for pid in PROTEINS:
            rows.append(one(pid, tmp))
            r = rows[-1]
            print(pid, "  ".join(f"{m} {r[m + '_beyond']:.3f} (p {r[m + '_p']:.3f})" for m in METHODS), flush=True)
    with open(os.path.join(HERE, "external.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    cols = ["pdb"] + [f"{m}_{k}" for m in METHODS for k in ("beyond", "p")]
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    lines += ["| " + " | ".join(r[c] if c == "pdb" else f"{r[c]:.3f}" for c in cols) + " |" for r in rows]
    with open(os.path.join(HERE, "external.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("Wrote results/external.csv and .md")


if __name__ == "__main__":
    main()
