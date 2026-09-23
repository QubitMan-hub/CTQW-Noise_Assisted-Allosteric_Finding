#!/usr/bin/env python3
"""Transport-vs-noise check for one protein: does signal from the source residue
peak at an intermediate dephasing rate (the environment-assisted hump)? Fast,
exact, classical reference the Qmod route must match.

    python classical_qsw.py NET.graphml --source A:151
"""
import argparse, os, sys
import numpy as np, networkx as nx
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import walk_core as wc


def main():
    ap = argparse.ArgumentParser(description="Transport vs noise on a residue network.")
    ap.add_argument("network")
    ap.add_argument("--source", required=True, help="Source residue id, e.g. A:151")
    ap.add_argument("--target", default=None, help="Target residue (default: farthest).")
    ap.add_argument("--site-energy", choices=["none", "hydropathy", "degree", "random"], default="hydropathy")
    ap.add_argument("--scale", type=float, default=3.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--gammas", default="0,0.05,0.1,0.2,0.4,0.8,1.5,3,6,12,25,50")
    ap.add_argument("--tmax", type=float, default=30.0)
    ap.add_argument("--ntime", type=int, default=200)
    ap.add_argument("--outdir", default="qsw_output")
    ap.add_argument("--prefix", default=None)
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    prefix = a.prefix or os.path.splitext(os.path.basename(a.network))[0] + "_qsw"
    G, nodes, A, resnames = wc.load_network(a.network)
    idx = {n: i for i, n in enumerate(nodes)}
    if a.source not in idx:
        sys.exit(f"[error] source '{a.source}' not in network, e.g. {nodes[:5]}")
    if a.target and a.target in idx:
        target = a.target
    else:
        target = max(nx.single_source_shortest_path_length(G, a.source).items(),
                     key=lambda kv: kv[1])[0]
    print(f"source {a.source} -> target {target} "
          f"(distance {nx.shortest_path_length(G, a.source, target)}), "
          f"site energies {a.site_energy} (scale {a.scale})")

    H = wc.build_hamiltonian(A, resnames, a.site_energy, a.scale, a.seed)
    gammas = [float(x) for x in a.gammas.split(",")]
    tlist = np.linspace(0, a.tmax, a.ntime)
    transport = []
    for g in gammas:
        v = wc.visiting_scores(wc.run_walk(H, idx[a.source], g, tlist), tlist)[idx[target]]
        transport.append(v)
        print(f"  gamma={g:6.2f}   transport={v:.4f}")

    with open(os.path.join(a.outdir, prefix + "_transport_vs_gamma.csv"), "w") as f:
        f.write("gamma,transport\n")
        f.writelines(f"{g},{v:.6f}\n" for g, v in zip(gammas, transport))

    peak = int(np.argmax(transport))
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.figure(figsize=(8, 5))
        plt.plot(gammas, transport, "o-", lw=2)
        plt.axvline(gammas[peak], color="crimson", ls="--", alpha=0.7, label=f"peak g={gammas[peak]:g}")
        plt.xlabel("dephasing rate gamma (0 = quantum, large = classical)")
        plt.ylabel(f"signal reaching {target}")
        plt.title(f"{prefix}: transport vs noise"); plt.legend(); plt.tight_layout()
        plt.savefig(os.path.join(a.outdir, prefix + "_transport_hump.png"), dpi=140); plt.close()
    except Exception as e:
        print(f"[warn] plot skipped: {e}")

    print("\n" + ("HUMP: transport peaks at intermediate noise." if 0 < peak < len(gammas) - 1
                  else "No hump: best at the quantum end." if peak == 0
                  else "No hump: best at the classical end."))


if __name__ == "__main__":
    main()
