#!/usr/bin/env python3
"""make_results.py - rebuild every table and figure in results/.

    python results/make_results.py            # rerun all walks (~25 min on 4 cores), then tables + figure
    python results/make_results.py --reuse    # rebuild tables + figure from results/runs/ only

Every protein studied is reported, including the ones where nothing is found.
Runs use the pipeline defaults (8 Å Cβ contacts, hydropathy energies at scale 3,
18 noise levels, 5 random-energy seeds, 10,000 label shuffles); the known sites
come from the bundled ALLO table. Only the small data files of each run are kept
(CSVs, parameters and result JSON); the network, structure and per-run figures
are left out to keep the repository light.
"""
import argparse, csv, glob, json, os, shutil, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import run_protein as rp      # noqa: E402
import sensitivity            # noqa: E402

PROTEINS = ["1T49", "3LSW", "1IWH", "3CSM"]      # every protein studied so far, positives and nulls
SENSITIVITY = ["1T49", "1IWH"]                   # cutoff x scale grid for the two with a noise hump
CUTOFFS, SCALES = [7.0, 8.0, 9.0], [1.0, 3.0, 5.0]
RUNS = os.path.join(HERE, "runs")
KEEP = ("*.csv", "*_parameters.json", "*_result.json")


def prune(folder):
    """Keep only the small data files of a run folder."""
    keep = {p for pat in KEEP for p in glob.glob(os.path.join(folder, pat))}
    for p in glob.glob(os.path.join(folder, "*")):
        if p not in keep:
            shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)


def run_all():
    for pid in PROTEINS:
        print(f"=== {pid}", flush=True)
        r = rp.analyze(pid, os.path.join(RUNS, pid), pid, log=lambda m: print("  " + m, flush=True))
        print(f"  done in {r['elapsed_s']} s", flush=True)
        prune(os.path.join(RUNS, pid))
    for pid in SENSITIVITY:
        print(f"=== sensitivity {pid}", flush=True)
        out = os.path.join(RUNS, "sensitivity", pid)
        shutil.rmtree(out, ignore_errors=True)
        sensitivity.grid(pid, CUTOFFS, SCALES, outdir=out, log=lambda m: print("  " + m, flush=True))
        for d in glob.glob(os.path.join(out, f"{pid}_c*")):     # per-setting run folders: numbers are in the CSV
            shutil.rmtree(d)


def protein_table():
    rows = []
    for pid in PROTEINS:
        with open(os.path.join(RUNS, pid, f"{pid}_result.json")) as f:
            r = json.load(f)
        a = r["allosteric"]
        d, st = a["adjusted"], a["adjusted"]["stats"]
        cl = d["baselines"]["classical walk, best rate"]
        rows.append({
            "pdb": pid, "protein": r["parameters"]["labels_used"]["protein"], "residues": r["summary"]["residues"],
            "active_residues": a["n_active"], "allosteric_residues": a["n_allosteric"],
            "raw_auc_best": a["stats"]["peak_value"], "proximity_auc": a["baselines"]["proximity to active site"],
            "beyond_distance_best": st["peak_value"], "best_gamma": st["peak_gamma"],
            "quantum_end": st["quantum_value"], "dephased_end": st["dephased_value"],
            "noise_gain": st["gain_abs"] if st["verdict"] == "hump" else 0.0,
            "p_value": d["significance"]["p_value"],
            "classical_best": cl, "quantum_minus_classical": st["peak_value"] - cl,
            "classical_p_value": a["classical"]["significance_adjusted"]["p_value"],
            "control_best_seed": max(d["control"]["best_seeds"]),
            "transport_verdict": r["transport"]["stats"]["verdict"],
            "transport_peak_gamma": r["transport"]["stats"]["peak_gamma"],
        })
    return rows


def write_table(rows, base, shown):
    with open(base + ".csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    cell = lambda v: f"{v:.3f}" if isinstance(v, float) else str(v)
    lines = ["| " + " | ".join(shown) + " |", "|" + "---|" * len(shown)]
    lines += ["| " + " | ".join(cell(r[k]) for k in shown) + " |" for r in rows]
    with open(base + ".md", "w") as f:
        f.write("\n".join(lines) + "\n")


def figure(path):
    """AUC beyond distance vs noise for every protein: quantum walk, random-energy band, classical walk."""
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fwd = lambda g: np.log10(1 + np.asarray(g) / 0.05)          # same x axis as the pipeline's figures
    inv = lambda u: 0.05 * (10 ** np.asarray(u) - 1)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "svg.fonttype": "none"})
    fig, axs = plt.subplots(2, 2, figsize=(12, 8.6), sharex=True, sharey=True)
    for ax, pid in zip(axs.flat, PROTEINS):
        with open(os.path.join(RUNS, pid, f"{pid}_result.json")) as f:
            r = json.load(f)
        g, d = np.array(r["gammas"]), r["allosteric"]["adjusted"]
        m, s = np.array(d["control"]["mean"]), np.array(d["control"]["sd"])
        ax.set_xscale("function", functions=(fwd, inv))
        ax.fill_between(g, m - s, m + s, color="#8a8a8a", alpha=0.2, lw=0)
        ax.plot(g, m, color="#8a8a8a", ls="--", lw=1.6, label="random site energies (mean ± sd, 5 seeds)")
        ax.plot(g, d["auc"], color="#111", lw=2.2, marker="o", ms=3.5, mec="white", mew=0.6,
                label="quantum walk, hydropathy site energies")
        ax.axhline(d["baselines"]["classical walk, best rate"], color="#444", lw=1.3, ls=(0, (1, 1.5)),
                   label="classical random walk (best hopping rate)")
        ax.axhline(0.5, color="#bbb", lw=1)
        ax.text(0.012, 0.506, "closeness alone", fontsize=8.5, color="#777")
        st = d["stats"]
        ax.set_title(f"{pid} ({r['summary']['residues']} residues): best {st['peak_value']:.3f} "
                     f"at γ = {st['peak_gamma']:g}, p = {d['significance']['p_value']:.3f}", fontsize=11)
        ax.set_ylim(0.3, 0.9)
        ax.set_xlim(0, 105)
        ax.set_xticks([0, 0.1, 1, 10, 100])
        ax.set_xticklabels(["0", "0.1", "1", "10", "100"])
        ax.minorticks_off()
        ax.grid(axis="y", color="#eee")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    for ax in axs[1]:
        ax.set_xlabel("dephasing rate γ")
    for ax in axs[:, 0]:
        ax.set_ylabel("ROC AUC among residues equally far\nfrom the active site")
    fig.legend(*axs[0, 0].get_legend_handles_labels(), loc="upper center", ncol=3, frameon=False, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path + ".png", dpi=300)
    fig.savefig(path + ".svg")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Rebuild the tables and figure in results/.")
    ap.add_argument("--reuse", action="store_true", help="Skip the walks; rebuild from results/runs/.")
    a = ap.parse_args()
    if not a.reuse:
        run_all()
    rows = protein_table()
    write_table(rows, os.path.join(HERE, "proteins"),
                ["pdb", "residues", "allosteric_residues", "beyond_distance_best", "best_gamma", "noise_gain",
                 "p_value", "classical_best", "control_best_seed"])
    for pid in SENSITIVITY:
        for ext in ("csv", "md"):
            shutil.copy(os.path.join(RUNS, "sensitivity", pid, f"{pid}_sensitivity.{ext}"),
                        os.path.join(HERE, f"sensitivity_{pid}.{ext}"))
    os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
    figure(os.path.join(HERE, "figures", "beyond_distance"))
    print("\nWrote results/proteins.csv/.md, results/sensitivity_*.csv/.md and results/figures/beyond_distance.png/.svg")


if __name__ == "__main__":
    main()
