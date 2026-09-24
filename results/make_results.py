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
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import run_protein as rp      # noqa: E402
import sensitivity            # noqa: E402

PROTEINS = ["1T49", "3LSW", "1IWH", "3CSM"]      # every protein studied so far, positives and nulls
SENSITIVITY = ["1T49", "1IWH"]                   # cutoff x scale grid for the two with a noise hump
CUTOFFS, SCALES = [7.0, 8.0, 9.0], [1.0, 3.0, 5.0]
RUNS = os.path.join(HERE, "runs")
KEEP = ("*.csv", "*_parameters.json", "*_result.json", "*_map.json")
MAP_KEYS = ("nodes", "edges", "gammas", "index", "rates", "signal", "classical", "aspect", "from")
FIG = os.path.join(HERE, "figures")


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
        with open(os.path.join(RUNS, pid, f"{pid}_map.json"), "w") as f:      # both walks at every gamma
            json.dump(rp.jsonable({k: r["map"][k] for k in MAP_KEYS}), f)
        prune(os.path.join(RUNS, pid))
    for pid in SENSITIVITY:
        print(f"=== sensitivity {pid}", flush=True)
        out = os.path.join(RUNS, "sensitivity", pid)
        shutil.rmtree(out, ignore_errors=True)
        sensitivity.grid(pid, CUTOFFS, SCALES, outdir=out, log=lambda m: print("  " + m, flush=True))
        for d in glob.glob(os.path.join(out, f"{pid}_c*")):     # per-setting run folders: numbers are in the CSV
            shutil.rmtree(d)


def load(pid, kind):
    with open(os.path.join(RUNS, pid, f"{pid}_{kind}.json")) as f:
        return json.load(f)


def difference_test(pid):
    """The pipeline's test of the quantum-minus-classical map (see run_protein.analyze)."""
    d = load(pid, "result")["allosteric"]["adjusted"]["difference"]
    return {"difference_gamma": d["gamma"], "difference_rate": d["classical_rate"], "difference_auc": d["auc"],
            "difference_p": d["p_value"], "allosteric_quantum_favoured": d["allosteric_quantum_favoured"]}


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
            **difference_test(pid),
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


def write_numbers(rows):
    """Every number the paper quotes, as LaTeX macros: \\res{1T49}{beyond}, \\sens{1T49}{nsig}, ...
    so the text can never drift from the runs."""
    f3, f4 = (lambda v: f"{v:.3f}"), (lambda v: f"{v:.4f}")
    lines = ["% written by results/make_results.py -- do not edit by hand"]
    put = lambda group, pid, key, val: lines.append(f"\\expandafter\\def\\csname {group}@{pid}@{key}\\endcsname{{{val}}}")
    for r in rows:
        res, pid = load(r["pdb"], "result"), r["pdb"]
        t = res["transport"]
        for key, val in {
                "protein": r["protein"], "residues": r["residues"], "nactive": r["active_residues"],
                "nallo": r["allosteric_residues"], "raw": f3(r["raw_auc_best"]), "prox": f3(r["proximity_auc"]),
                "beyond": f3(r["beyond_distance_best"]), "gamma": f"{r['best_gamma']:g}",
                "qend": f3(r["quantum_end"]), "dend": f3(r["dephased_end"]), "gain": f3(r["noise_gain"]),
                "p": f3(r["p_value"]), "classical": f3(r["classical_best"]), "clp": f3(r["classical_p_value"]),
                "margin": f"{r['quantum_minus_classical']:+.3f}", "seedbest": f3(r["control_best_seed"]),
                "tgamma": f"{r['transport_peak_gamma']:g}", "tgain": f"{100 * t['stats']['gain_rel']:.1f}",
                "tseeds": f"{t['control']['seed_verdicts'].count('hump')}", "tn": f"{t['control']['n_seeds']}",
                "dgamma": f"{r['difference_gamma']:g}", "drate": f"{r['difference_rate']:.3g}",
                "dauc": f3(r["difference_auc"]), "dp": f3(r["difference_p"]),
                "dfav": f"{round(100 * r['allosteric_quantum_favoured'])}",
                "shells": res["allosteric"]["adjusted"]["shells"], "elapsed": f"{res['elapsed_s']:.0f}",
                "clraw": f3(res["allosteric"]["baselines"]["classical walk, best rate"])}.items():
            put("res", pid, key, val)
        put("res", pid, "rateheavy", f4(res["quantum_vs_classical"]["classical_rates"][-1]))
    for pid in SENSITIVITY:
        with open(os.path.join(HERE, f"sensitivity_{pid}.csv")) as f:
            srows = list(csv.DictReader(f))
        ps = [float(x["p_value"]) for x in srows]
        margins = [float(x["quantum_minus_classical"]) for x in srows]
        best = [float(x["beyond_distance_best"]) for x in srows]
        for key, val in {"n": len(srows), "nsig": sum(p < 0.05 for p in ps), "nbeat": sum(mg > 0 for mg in margins),
                         "nhump": sum(x["noise_hump"] == "yes" for x in srows), "pmin": f3(min(ps)), "pmax": f3(max(ps)),
                         "bmin": f3(min(best)), "bmax": f3(max(best)), "mmin": f"{min(margins):+.3f}",
                         "mmax": f"{max(margins):+.3f}"}.items():
            put("sens", pid, key, val)
    with open(os.path.join(HERE, "numbers.tex"), "w") as f:
        f.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------- figures (greyscale, 300 dpi PNG + SVG)
def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5, "svg.fonttype": "none",
                         "axes.spines.top": False, "axes.spines.right": False})
    return plt


def _gamma_axis(ax):
    rp.gamma_axis(ax, 100.0)                  # the same axis as the pipeline's own figures
    ax.grid(axis="y", color="#eee")


def _save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=300)
    fig.savefig(os.path.join(FIG, name + ".svg"))


def fig_workflow():
    """The method in one picture: structure -> network -> walks over noise -> scores and tests."""
    plt = _plt()
    fig, ax = plt.subplots(figsize=(11, 2.9))
    ax.set_xlim(0, 11)
    ax.set_ylim(0, 2.9)
    ax.axis("off")
    boxes = [("Structure", "PDB id or file;\nknown sites from ALLO"),
             ("Contact network", "residues = nodes;\nCβ within 8 Å = edges"),
             ("Hamiltonian", "H = A + s·diag(ε)\nε: z-scored hydropathy"),
             ("Dephased walk", "ρ from the active site\n18 rates γ, 0 to 100"),
             ("Scores", "signal per residue\n= ∫ ρ_ii dt, t ≤ 30"),
             ("Tests", "AUC beyond distance, p;\nclassical walk; random ε")]
    w, gap = 1.62, 0.25
    for k, (title, body) in enumerate(boxes):
        x = 0.1 + k * (w + gap)
        ax.add_patch(plt.Rectangle((x, 0.55), w, 1.8, fc="#f4f4f4" if k % 2 else "#e9e9e9", ec="#222", lw=1.2))
        ax.text(x + w / 2, 2.0, title, ha="center", va="center", fontsize=11, weight="bold")
        ax.text(x + w / 2, 1.25, body, ha="center", va="center", fontsize=9, color="#333", linespacing=1.4)
        if k < len(boxes) - 1:
            ax.annotate("", xy=(x + w + gap - 0.02, 1.45), xytext=(x + w + 0.02, 1.45),
                        arrowprops=dict(arrowstyle="-|>", color="#222", lw=1.3))
    ax.text(5.5, 0.18, "Nulls on every run: 5 random-energy seeds; a classical random walk on the same contacts; "
            "10,000 label shuffles within distance shells", ha="center", fontsize=9, color="#555")
    fig.tight_layout()
    _save(fig, "workflow")
    plt.close(fig)


def fig_transport():
    """Transport to the distal residues versus noise, with the random-energy band."""
    plt = _plt()
    fig, axs = plt.subplots(1, 4, figsize=(13, 3.3), sharex=True)
    for ax, pid in zip(axs, PROTEINS):
        r = load(pid, "result")
        g, t = np.array(r["gammas"]), r["transport"]
        m, sd = np.array(t["control"]["mean"]), np.array(t["control"]["sd"])
        _gamma_axis(ax)
        ax.fill_between(g, m - sd, m + sd, color="#8a8a8a", alpha=0.2, lw=0)
        ax.plot(g, m, color="#8a8a8a", ls="--", lw=1.5, label="random site energies (mean ± sd)")
        ax.plot(g, t["distal_mean"], color="#111", lw=2, marker="o", ms=3, mec="white", mew=0.5,
                label="hydropathy site energies")
        st = t["stats"]
        ax.set_title(f"{pid}: peak at γ = {st['peak_gamma']:g}", fontsize=10.5)
        ax.set_xlabel("dephasing rate γ")
    axs[0].set_ylabel("mean signal on distal residues")
    fig.legend(*axs[0].get_legend_handles_labels(), loc="upper center", ncol=2, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    _save(fig, "transport")
    plt.close(fig)


def fig_beyond_distance():
    """AUC beyond distance vs noise for every protein: quantum walk, random-energy band, classical walk."""
    plt = _plt()
    fig, axs = plt.subplots(2, 2, figsize=(12, 8.4), sharex=True, sharey=True)
    for ax, pid in zip(axs.flat, PROTEINS):
        r = load(pid, "result")
        g, d = np.array(r["gammas"]), r["allosteric"]["adjusted"]
        m, sd = np.array(d["control"]["mean"]), np.array(d["control"]["sd"])
        _gamma_axis(ax)
        ax.fill_between(g, m - sd, m + sd, color="#8a8a8a", alpha=0.2, lw=0)
        ax.plot(g, m, color="#8a8a8a", ls="--", lw=1.6, label="random site energies (mean ± sd, 5 seeds)")
        ax.plot(g, d["auc"], color="#111", lw=2.2, marker="o", ms=3.5, mec="white", mew=0.6,
                label="quantum walk, hydropathy site energies")
        ax.axhline(d["baselines"]["classical walk, best rate"], color="#444", lw=1.3, ls=(0, (1, 1.5)),
                   label="classical random walk (best hopping rate)")
        ax.axhline(0.5, color="#bbb", lw=1)
        ax.text(0.012, 0.506, "closeness alone", fontsize=8.5, color="#777")
        st = d["stats"]
        ax.set_title(f"{pid} ({r['summary']['residues']} residues): best {st['peak_value']:.3f} "
                     f"at γ = {st['peak_gamma']:g}, p = {d['significance']['p_value']:.3f}", fontsize=10.5)
        ax.set_ylim(0.3, 0.9)
    for ax in axs[1]:
        ax.set_xlabel("dephasing rate γ")
    for ax in axs[:, 0]:
        ax.set_ylabel("ROC AUC among residues equally far\nfrom the active site")
    fig.legend(*axs[0, 0].get_legend_handles_labels(), loc="upper center", ncol=3, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    _save(fig, "beyond_distance")
    plt.close(fig)


def fig_sensitivity():
    """Cutoff x scale grids: best AUC beyond distance (shade), p and the margin over the classical walk."""
    plt = _plt()
    fig, axs = plt.subplots(1, len(SENSITIVITY), figsize=(5.2 * len(SENSITIVITY), 4.2))
    for ax, pid in zip(np.atleast_1d(axs), SENSITIVITY):
        with open(os.path.join(HERE, f"sensitivity_{pid}.csv")) as f:
            rows = list(csv.DictReader(f))
        cell = {(float(r["cutoff"]), float(r["scale"])): r for r in rows}
        grid = np.array([[float(cell[(c, s)]["beyond_distance_best"]) for s in SCALES] for c in CUTOFFS])
        ax.imshow(grid, cmap="Greys", vmin=0.45, vmax=0.95)
        for i, c in enumerate(CUTOFFS):
            for k, s in enumerate(SCALES):
                r = cell[(c, s)]
                v, p = float(r["beyond_distance_best"]), float(r["p_value"])
                margin = float(r["quantum_minus_classical"])
                ax.text(k, i, f"{v:.3f}\np = {p:.3f}\n{'+' if margin >= 0 else '−'}{abs(margin):.3f} vs cl.",
                        ha="center", va="center", fontsize=9, color="white" if v > 0.72 else "#111")
        ax.set_xticks(range(len(SCALES)), [f"{s:g}" for s in SCALES])
        ax.set_yticks(range(len(CUTOFFS)), [f"{c:g} Å" for c in CUTOFFS])
        ax.set_xlabel("site-energy scale s")
        ax.set_ylabel("contact cutoff")
        ax.set_title(pid)
        for side in ("left", "bottom"):
            ax.spines[side].set_visible(False)
        ax.tick_params(length=0)
    fig.tight_layout()
    _save(fig, "sensitivity")
    plt.close(fig)


def fig_map(pid="1T49"):
    """The residue map at the gamma where the AUC beyond distance peaks: (a) signal, (b) quantum minus classical."""
    plt = _plt()
    if not os.path.isfile(os.path.join(RUNS, pid, f"{pid}_map.json")):
        sys.exit(f"results/runs/{pid}/{pid}_map.json is missing: run without --reuse once to create it.")
    r, m = load(pid, "result"), load(pid, "map")
    j = r["allosteric"]["adjusted"]["stats"]["peak_index"]
    nodes = m["nodes"]
    xy = np.array([[v["x"], -v["y"]] for v in nodes])
    q, cl = np.array(m["signal"][j]), np.array(m["classical"][j])
    diff = q - cl
    src = np.array([v["role"] == "source" for v in nodes])
    allo = np.array([v["role"] == "allosteric" for v in nodes])
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.6))
    for ax, mode in zip(axs, ("signal", "diff")):
        for a, b in m["edges"]:
            ax.plot(*xy[[a, b]].T, color="#dcdcdc", lw=0.5, zorder=1)
        if mode == "signal":
            lo, hi = q[~src].min(), q[~src].max()
            ax.scatter(*xy[~src].T, c=(q[~src] - lo) / (hi - lo), cmap="Greys", vmin=-0.08, vmax=1, s=34,
                       edgecolors="white", linewidths=0.6, zorder=2)
            ax.set_title(f"(a) signal from the active site, γ = {m['gammas'][j]:g}", fontsize=11)
        else:
            t = np.abs(diff) / np.abs(diff[~src]).max()
            more, less = ~src & (diff >= 0), ~src & (diff < 0)
            grey = lambda v: plt.cm.Greys(0.15 + 0.85 * v)
            ax.scatter(*xy[more].T, c=[grey(v) for v in t[more]], s=34, edgecolors="white", linewidths=0.6, zorder=2)
            ax.scatter(*xy[less].T, facecolors="white", edgecolors=[grey(v) for v in t[less]], s=34,
                       linewidths=1.4, zorder=2)
            ax.set_title(f"(b) quantum minus classical (classical rate k = {m['rates'][j]:.3g})", fontsize=11)
        ax.scatter(*xy[src].T, c="#111", s=34, zorder=3)
        ax.scatter(*xy[src].T, facecolors="none", edgecolors="#111", s=150, linewidths=1.2, zorder=3)
        ax.scatter(*xy[allo].T, facecolors="none", edgecolors="#666", s=150, linewidths=1.2,
                   linestyles=(0, (2, 1.5)), zorder=3)
        ax.set_aspect("equal")
        ax.axis("off")
    from matplotlib.lines import Line2D
    keys = [Line2D([], [], marker="o", ls="", mfc="none", mec="#111", ms=11, label="active site (walk starts here)"),
            Line2D([], [], marker="o", ls="", mfc="none", mec="#666", ms=11, label="known allosteric residue"),
            Line2D([], [], marker="o", ls="", mfc="#333", mec="white", ms=7, label="(b) quantum puts more signal here"),
            Line2D([], [], marker="o", ls="", mfc="white", mec="#333", mew=1.4, ms=7, label="(b) classical puts more here")]
    fig.legend(handles=keys, loc="lower center", ncol=4, frameon=False, fontsize=9.5)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    _save(fig, f"map_{pid}")
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
                 "p_value", "classical_best", "control_best_seed", "difference_auc", "difference_p"])
    for pid in SENSITIVITY:
        for ext in ("csv", "md"):
            shutil.copy(os.path.join(RUNS, "sensitivity", pid, f"{pid}_sensitivity.{ext}"),
                        os.path.join(HERE, f"sensitivity_{pid}.{ext}"))
    write_numbers(rows)
    for make in (fig_workflow, fig_transport, fig_beyond_distance, fig_sensitivity, fig_map):
        make()
    print("\nWrote results/proteins.csv/.md, results/sensitivity_*.csv/.md and results/figures/")


if __name__ == "__main__":
    main()
