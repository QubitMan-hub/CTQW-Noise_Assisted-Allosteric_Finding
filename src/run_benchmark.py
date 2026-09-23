#!/usr/bin/env python3
"""The headline experiment: run the allosteric-recovery-vs-noise sweep across
many labelled proteins and aggregate. Reports how often intermediate noise wins,
the mean AUC gained, and the mean AUC-vs-noise curve (the paper's figure).

Input: a folder with NAME.graphml + NAME.json (or NAME_labels.json) per protein.

    python run_benchmark.py --data-dir data/benchmark
"""
import argparse, glob, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import walk_core as wc
from score_allostery import resolve_labels


def find_pairs(data_dir):
    pairs = []
    for g in sorted(glob.glob(os.path.join(data_dir, "*.graphml"))):
        name = os.path.splitext(os.path.basename(g))[0]
        for cand in (name + ".json", name + "_labels.json"):
            lp = os.path.join(data_dir, cand)
            if os.path.isfile(lp):
                pairs.append((name, g, lp))
                break
    return pairs


def main():
    ap = argparse.ArgumentParser(description="Aggregate allosteric recovery across proteins.")
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--site-energy", choices=["none", "hydropathy", "degree", "random"], default="hydropathy")
    ap.add_argument("--scale", type=float, default=3.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--gammas", default="0,0.05,0.1,0.2,0.4,0.8,1.5,3,6,12,25,50")
    ap.add_argument("--tmax", type=float, default=30.0)
    ap.add_argument("--ntime", type=int, default=200)
    ap.add_argument("--min-hops", type=int, default=0)
    ap.add_argument("--outdir", default="benchmark_output")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    gammas = [float(x) for x in a.gammas.split(",")]
    tlist = np.linspace(0, a.tmax, a.ntime)
    pairs = find_pairs(a.data_dir)
    if not pairs:
        sys.exit(f"[error] no NAME.graphml + NAME.json pairs in {a.data_dir}")
    print(f"found {len(pairs)} labelled proteins\n")

    rows, curves = [], []
    for name, gpath, lpath in pairs:
        try:
            G, nodes, A, resnames = wc.load_network(gpath)
            src, positive, eligible, _ = resolve_labels(G, nodes, lpath, a.min_hops)
            aucs = wc.sweep_auc(wc.build_hamiltonian(A, resnames, a.site_energy, a.scale, a.seed),
                                src, positive, eligible, gammas, tlist)
        except Exception as e:
            print(f"  {name:12s}  skipped ({e})")
            continue
        peak = int(np.nanargmax(aucs))
        q0, qmax, best = aucs[0], aucs[-1], aucs[peak]
        hump = bool(0 < peak < len(gammas) - 1 and best > q0 + 1e-6 and best > qmax + 1e-6)
        rows.append((name, len(nodes), int(positive.sum()), q0, qmax, best,
                     gammas[peak], hump, best - max(q0, qmax)))
        curves.append(aucs)
        print(f"  {name:12s}  quantum={q0:.3f}  classical={qmax:.3f}  "
              f"best={best:.3f}@g{gammas[peak]:<5g} {'HUMP' if hump else '----'}  "
              f"gain={best - max(q0, qmax):+.3f}")
    if not rows:
        sys.exit("[error] no protein produced a usable result.")

    with open(os.path.join(a.outdir, "benchmark_summary.csv"), "w") as f:
        f.write("protein,residues,n_allosteric,auc_quantum,auc_classical,auc_best,"
                "best_gamma,intermediate_hump,auc_gain\n")
        f.writelines("{},{},{},{:.4f},{:.4f},{:.4f},{},{},{:.4f}\n".format(*r) for r in rows)

    curves = np.array(curves)
    mean, std = np.nanmean(curves, axis=0), np.nanstd(curves, axis=0)
    n, n_hump = len(rows), sum(r[7] for r in rows)
    agg_peak = int(np.nanargmax(mean))
    print("\n" + "=" * 58)
    print(f"proteins evaluated       : {n}")
    print(f"intermediate-noise wins  : {n_hump}/{n} ({100*n_hump/n:.0f}%)")
    print(f"mean AUC gain from noise  : {np.mean([r[8] for r in rows]):+.3f}")
    print(f"mean curve peaks at gamma : {gammas[agg_peak]:g} (AUC {mean[agg_peak]:.3f} "
          f"vs quantum {mean[0]:.3f}, classical {mean[-1]:.3f})")
    print("=" * 58)

    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.figure(figsize=(8, 5))
        plt.plot(gammas, mean, "o-", lw=2, color="navy", label="mean AUC")
        plt.fill_between(gammas, mean - std, mean + std, alpha=0.2, color="navy", label="+/- 1 s.d.")
        plt.axhline(0.5, color="gray", ls=":", alpha=0.7)
        plt.axvline(gammas[agg_peak], color="crimson", ls="--", alpha=0.7, label=f"mean peak g={gammas[agg_peak]:g}")
        plt.xlabel("dephasing rate gamma (0 = quantum, large = classical)")
        plt.ylabel("ROC AUC for known allosteric residues")
        plt.title(f"allosteric recovery vs noise, averaged over {n} proteins")
        plt.legend(); plt.tight_layout()
        plt.savefig(os.path.join(a.outdir, "benchmark_mean_curve.png"), dpi=150); plt.close()
    except Exception as e:
        print(f"[warn] plot skipped: {e}")
    print(f"\nsummary: {os.path.join(a.outdir, 'benchmark_summary.csv')}")
    print(f"figure:  {os.path.join(a.outdir, 'benchmark_mean_curve.png')}")


if __name__ == "__main__":
    main()
