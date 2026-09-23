#!/usr/bin/env python3
"""Does intermediate dephasing recover the KNOWN allosteric residues (per protein)?
Labels JSON: {"source": ["A:57"], "allosteric": ["A:196", ...]}. Real labels:
ALLO benchmark supplementary (Wu et al. 2022, Europe PMC). Build with make_labels.py.

    python score_allostery.py NET.graphml --labels labels.json
"""
import argparse, json, os, sys
import numpy as np, networkx as nx
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import walk_core as wc


def resolve_labels(G, nodes, labels_path, min_hops):
    """-> (source_indices, positive_mask, eligible_mask, report). Eligible =
    non-source residues (and, if min_hops>0, not within min_hops of a source)."""
    idx = {n: i for i, n in enumerate(nodes)}
    lab = json.load(open(labels_path))
    sources = [s for s in lab.get("source", []) if s in idx]
    allo = {s for s in lab.get("allosteric", []) if s in idx}
    if not sources or not allo:
        raise ValueError("no source/allosteric residues match the network node ids.")
    eligible = np.ones(len(nodes), dtype=bool)
    for s in sources:
        eligible[idx[s]] = False
    if min_hops > 0:
        for s in sources:
            for c in nx.single_source_shortest_path_length(G, s, cutoff=min_hops - 1):
                eligible[idx[c]] = False
    positive = np.array([n in allo for n in nodes], dtype=bool)
    if int((positive & eligible).sum()) == 0:
        raise ValueError("no allosteric residues left after the filter (try --min-hops 0).")
    report = f"{len(nodes)} residues | sources {len(sources)} | allosteric in scope {int((positive & eligible).sum())}"
    return [idx[s] for s in sources], positive, eligible, report


def main():
    ap = argparse.ArgumentParser(description="Allosteric recovery vs noise.")
    ap.add_argument("network")
    ap.add_argument("--labels", required=True)
    ap.add_argument("--site-energy", choices=["none", "hydropathy", "degree", "random"], default="hydropathy")
    ap.add_argument("--scale", type=float, default=3.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--gammas", default="0,0.05,0.1,0.2,0.4,0.8,1.5,3,6,12,25,50")
    ap.add_argument("--tmax", type=float, default=30.0)
    ap.add_argument("--ntime", type=int, default=200)
    ap.add_argument("--min-hops", type=int, default=0)
    ap.add_argument("--outdir", default="score_output")
    ap.add_argument("--prefix", default=None)
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    prefix = a.prefix or os.path.splitext(os.path.basename(a.network))[0]
    G, nodes, A, resnames = wc.load_network(a.network)
    src, positive, eligible, report = resolve_labels(G, nodes, a.labels, a.min_hops)
    print(f"{report} | site energies {a.site_energy} (scale {a.scale})")

    H = wc.build_hamiltonian(A, resnames, a.site_energy, a.scale, a.seed)
    gammas = [float(x) for x in a.gammas.split(",")]
    aucs = wc.sweep_auc(H, src, positive, eligible, gammas, np.linspace(0, a.tmax, a.ntime))
    for g, auc in zip(gammas, aucs):
        print(f"  gamma={g:6.2f}   AUC={auc:.3f}")

    with open(os.path.join(a.outdir, prefix + "_auc_vs_gamma.csv"), "w") as f:
        f.write("gamma,auc\n")
        f.writelines(f"{g},{auc:.6f}\n" for g, auc in zip(gammas, aucs))

    peak = int(np.nanargmax(aucs))
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.figure(figsize=(8, 5))
        plt.plot(gammas, aucs, "o-", lw=2)
        plt.axhline(0.5, color="gray", ls=":", alpha=0.7, label="random")
        plt.axvline(gammas[peak], color="crimson", ls="--", alpha=0.7, label=f"best g={gammas[peak]:g}")
        plt.xlabel("dephasing rate gamma (0 = quantum, large = classical)")
        plt.ylabel("ROC AUC for known allosteric residues")
        plt.title(f"{prefix}: allosteric recovery vs noise"); plt.legend(); plt.tight_layout()
        plt.savefig(os.path.join(a.outdir, prefix + "_auc_vs_gamma.png"), dpi=140); plt.close()
    except Exception as e:
        print(f"[warn] plot skipped: {e}")

    q0, qmax, best = aucs[0], aucs[-1], aucs[peak]
    print("\n" + (f"HUMP: best at intermediate gamma={gammas[peak]:g} (AUC {best:.3f}) "
                  f"vs quantum {q0:.3f}, classical {qmax:.3f}."
                  if 0 < peak < len(gammas) - 1 and best > q0 + 1e-6 and best > qmax + 1e-6
                  else f"No hump: best fully quantum (AUC {q0:.3f})." if peak == 0
                  else "No hump: best toward the classical end."))
    print("(One protein is one data point; use run_benchmark.py across many.)")


if __name__ == "__main__":
    main()
