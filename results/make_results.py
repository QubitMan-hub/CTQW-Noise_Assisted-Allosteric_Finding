#!/usr/bin/env python3
"""make_results.py - rebuild every table and figure in results/.

    python results/make_results.py            # rerun all walks (about 1.5 h on 4 cores), then tables + figures
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

DEVELOPMENT = ["1T49", "3LSW", "1IWH", "3CSM"]   # picked during development, positives and nulls
VALIDATION = ["2RD5", "3HO6", "4B1F", "4BBG", "4PFK", "3PYY",    # prespecified: select_proteins.py,
              "3ZCW", "3HFR", "3M3F", "3H30"]                # in its fixed draw order
REPLICATION = ["2YHD", "4HO6", "2VD3", "3PXF", "1ZDS"]            # pre-registered: PREREGISTRATION.md
PROTEINS = DEVELOPMENT + VALIDATION + REPLICATION
HELDOUT = VALIDATION + REPLICATION
SAME_PROTEIN = {"3ZCW": "4BBG", "3HFR": "4B1F", "3M3F": "3LSW"}   # other structures of one protein
KEY_PROTEINS = ["1T49", "3PYY", "3H30", "2RD5"]    # the main-text figure: one of each kind of result
MARK = {"validation": "*", "replication": "†"}          # set marks in figure titles
SET_OF = {**{p: "development" for p in DEVELOPMENT}, **{p: "validation" for p in VALIDATION},
          **{p: "replication" for p in REPLICATION}}
SENSITIVITY = ["1T49", "1IWH"]                   # cutoff x scale grid for the two with a noise hump
CUTOFFS, SCALES = [7.0, 8.0, 9.0], [1.0, 3.0, 5.0]
RUNS = os.path.join(HERE, "runs")
KEEP = ("*.csv", "*_parameters.json", "*_result.json", "*_map.json")
MAP_KEYS = ("nodes", "edges", "gammas", "index", "rates", "signal", "classical", "aspect", "from")
FIG = os.path.join(HERE, "figures")
EXTERNAL = ("closeness", "betweenness", "prs")      # established predictors, see external_baselines.py


def prune(folder):
    """Keep only the small data files of a run folder."""
    keep = {p for pat in KEEP for p in glob.glob(os.path.join(folder, pat))}
    for p in glob.glob(os.path.join(folder, "*")):
        if p not in keep:
            shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)


def run_all(proteins=None, with_grid=True):
    for pid in proteins or PROTEINS:
        print(f"=== {pid}", flush=True)
        r = rp.analyze(pid, os.path.join(RUNS, pid), pid, log=lambda m: print("  " + m, flush=True))
        print(f"  done in {r['elapsed_s']} s", flush=True)
        with open(os.path.join(RUNS, pid, f"{pid}_map.json"), "w") as f:      # both walks at every gamma
            json.dump(rp.jsonable({k: r["map"][k] for k in MAP_KEYS}), f)
        prune(os.path.join(RUNS, pid))
    for pid in SENSITIVITY if with_grid else []:
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


def bootstrap_ci(pid, n_boot=2000, seed=0):
    """95% interval for the best AUC beyond distance, at the gamma where it peaks, from the
    stored run: residues are resampled with replacement within each distance shell (so the
    shells keep their sizes) and the within-shell score is recomputed each time. The gamma
    is held at the selected value, so the interval does not include the choice of gamma."""
    from scipy.stats import rankdata
    wc = rp.wc
    r, m = load(pid, "result"), load(pid, "map")
    ids = [v["id"] for v in m["nodes"]]
    with open(os.path.join(RUNS, pid, f"{pid}_quantum_vs_classical.csv")) as f:
        rows = {row["residue"]: row for row in csv.DictReader(f)}
    hops = np.array([int(rows[i]["hops_from_start"]) if i in rows else 0 for i in ids])
    reach = np.array([i in rows for i in ids]) & (hops >= 0)
    pos = np.array([i in rows and rows[i]["known_allosteric"] == "yes" for i in ids])
    shells = wc.distance_shells(hops, reach)
    score = np.array(m["signal"][r["allosteric"]["adjusted"]["stats"]["peak_index"]])
    point = wc.roc_auc(wc.shell_percentile(score, shells, len(ids))[reach], pos[reach])
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        q, lab = [], []
        for g in shells:
            pick = rng.choice(g, len(g))
            q.append((rankdata(score[pick]) - 0.5) / len(g))
            lab.append(pos[pick])
        q, lab = np.concatenate(q), np.concatenate(lab)
        if lab.any() and not lab.all():
            boots.append(wc.roc_auc(q, lab))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return point, float(lo), float(hi)


def protein_table():
    rows = []
    for pid in PROTEINS:
        with open(os.path.join(RUNS, pid, f"{pid}_result.json")) as f:
            r = json.load(f)
        a = r["allosteric"]
        d, st = a["adjusted"], a["adjusted"]["stats"]
        cl = d["baselines"]["classical walk, best rate"]
        ew = a["incoherent"]
        rows.append({
            "pdb": pid, "set": SET_OF[pid], "protein": r["parameters"]["labels_used"]["protein"], "residues": r["summary"]["residues"],
            "active_residues": a["n_active"], "allosteric_residues": a["n_allosteric"],
            "raw_auc_best": a["stats"]["peak_value"], "proximity_auc": a["baselines"]["proximity to active site"],
            "beyond_distance_best": st["peak_value"], "best_gamma": st["peak_gamma"],
            "quantum_end": st["quantum_value"], "dephased_end": st["dephased_value"],
            "noise_gain": st["gain_abs"] if st["verdict"] == "hump" else 0.0,
            "p_value": d["significance"]["p_value"],
            "classical_best": cl, "quantum_minus_classical": st["peak_value"] - cl,
            "classical_p_value": a["classical"]["significance_adjusted"]["p_value"],
            "energy_weighted_best": ew["significance_adjusted"]["observed_best"],
            "energy_weighted_gamma": ew["best_gamma_adjusted"],
            "energy_weighted_p_value": ew["significance_adjusted"]["p_value"],
            "quantum_minus_energy_weighted": st["peak_value"] - ew["significance_adjusted"]["observed_best"],
            "control_best_seed": max(d["control"]["best_seeds"]),
            "transport_verdict": r["transport"]["stats"]["verdict"],
            "transport_peak_gamma": r["transport"]["stats"]["peak_gamma"],
            **difference_test(pid),
        })
        point, lo, hi = bootstrap_ci(pid)
        if abs(point - st["peak_value"]) > 1e-6:           # the stored data must reproduce the stored score
            sys.exit(f"{pid}: bootstrap data give {point:.6f}, stored {st['peak_value']:.6f}")
        rows[-1].update(ci_low=lo, ci_high=hi)
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
    fp = lambda v: "\\ensuremath{<}0.001" if v < 0.001 else f"{v:.3f}"     # 10,000 shuffles: p is never 0
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
                "p": fp(r["p_value"]), "classical": f3(r["classical_best"]), "clp": fp(r["classical_p_value"]),
                "margin": f"{r['quantum_minus_classical']:+.3f}",
                "ew": f3(r["energy_weighted_best"]), "ewp": fp(r["energy_weighted_p_value"]),
                "ewgamma": f"{r['energy_weighted_gamma']:g}", "ewmargin": f"{r['quantum_minus_energy_weighted']:+.3f}", "seedbest": f3(r["control_best_seed"]),
                "tgamma": f"{r['transport_peak_gamma']:g}", "tgain": f"{100 * t['stats']['gain_rel']:.1f}",
                "tseeds": f"{t['control']['seed_verdicts'].count('hump')}", "tn": f"{t['control']['n_seeds']}",
                "dgamma": f"{r['difference_gamma']:g}", "drate": f"{r['difference_rate']:.3g}",
                "dauc": f3(r["difference_auc"]), "dp": fp(r["difference_p"]),
                "dfav": f"{round(100 * r['allosteric_quantum_favoured'])}",
                "shells": res["allosteric"]["adjusted"]["shells"], "elapsed": f"{res['elapsed_s']:.0f}",
                "clraw": f3(res["allosteric"]["baselines"]["classical walk, best rate"])}.items():
            put("res", pid, key, val)
        for key, val in {"p": r["p_value"], "clp": r["classical_p_value"], "ewp": r["energy_weighted_p_value"],
                         "dp": r["difference_p"]}.items():       # with the relation, for "$p\\res{..}{peq}$" in math mode
            put("res", pid, key + "eq", "<0.001" if val < 0.001 else f"={val:.3f}")
        put("res", pid, "rateheavy", f4(res["quantum_vs_classical"]["classical_rates"][-1]))
    for pid in SENSITIVITY:
        with open(os.path.join(HERE, f"sensitivity_{pid}.csv")) as f:
            srows = list(csv.DictReader(f))
        ps = [float(x["p_value"]) for x in srows]
        margins = [float(x["quantum_minus_classical"]) for x in srows]
        ewm = [float(x["quantum_minus_energy_weighted"]) for x in srows]
        best = [float(x["beyond_distance_best"]) for x in srows]
        for key, val in {"n": len(srows), "nsig": sum(p < 0.05 for p in ps), "nbeat": sum(mg > 0 for mg in margins),
                         "nhump": sum(x["noise_hump"] == "yes" for x in srows), "pmin": fp(min(ps)), "pmax": fp(max(ps)),
                         "bmin": f3(min(best)), "bmax": f3(max(best)), "mmin": f"{min(margins):+.3f}",
                         "mmax": f"{max(margins):+.3f}", "ewmin": f"{min(ewm):+.3f}", "ewmax": f"{max(ewm):+.3f}",
                         "nbeatew": sum(mg > 0.02 for mg in ewm)}.items():
            put("sens", pid, key, val)
    # the prespecified validation set as a whole: how many significant, and how often that happens by chance
    from scipy.stats import binom
    for name, members in (("dev", DEVELOPMENT), ("val", VALIDATION), ("rep", REPLICATION), ("held", HELDOUT),
                          ("all", PROTEINS)):
        sub = [r for r in rows if r["pdb"] in members]
        nsig = sum(r["p_value"] < 0.05 for r in sub)
        for key, val in {"n": len(sub), "nsig": nsig, "binp": fp(float(binom.sf(nsig - 1, len(sub), 0.05))),
                         "nbeatcl": sum(r["quantum_minus_classical"] > 0.02 for r in sub),
                         "nbeatew": sum(r["quantum_minus_energy_weighted"] > 0.02 for r in sub),
                         "nsigew": sum(r["energy_weighted_p_value"] < 0.05 for r in sub),
                         "nsigcl": sum(r["classical_p_value"] < 0.05 for r in sub),
                         "nhump": sum(r["noise_gain"] > 0.01 for r in sub),
                         "nseed": sum(r["beyond_distance_best"] > r["control_best_seed"] for r in sub),
                         "mean": f3(float(np.mean([r["beyond_distance_best"] for r in sub]))),
                         "meanew": f3(float(np.mean([r["energy_weighted_best"] for r in sub]))),
                         "meancl": f3(float(np.mean([r["classical_best"] for r in sub]))),
                         "bonf": fp(min(1.0, len(sub) * min(r["p_value"] for r in sub))),
                         "thump": sum(r["transport_verdict"] == "hump" for r in sub),
                         "tseedhump": sum(load(r["pdb"], "result")["transport"]["control"]["seed_verdicts"].count("hump")
                                          for r in sub),
                         "tseedn": sum(load(r["pdb"], "result")["transport"]["control"]["n_seeds"] for r in sub)}.items():
            put("agg", name, key, val)
    first = [r for r in rows if r["pdb"] in VALIDATION[:6]]              # the first validation stage
    put("agg", "val6", "nsig", sum(r["p_value"] < 0.05 for r in first))
    put("agg", "val6", "binp", fp(float(binom.sf(sum(r["p_value"] < 0.05 for r in first) - 1, len(first), 0.05))))
    for r in rows:
        put("res", r["pdb"], "cilo", f3(r["ci_low"]))
        put("res", r["pdb"], "cihi", f3(r["ci_high"]))
    # what survives every control: the quantum walk is significant, and the energy-weighted walk beats the
    # plain classical walk (by > 0.02) and every random-energy seed, and is significant itself
    survive = lambda r: (r["p_value"] < 0.05 and r["energy_weighted_p_value"] < 0.05
                         and r["energy_weighted_best"] > r["classical_best"] + 0.02
                         and r["energy_weighted_best"] > r["control_best_seed"])
    # beats both classical walks and every random seed with the quantum walk itself (interference would show here)
    quantum_only = lambda r: (r["p_value"] < 0.05 and r["quantum_minus_energy_weighted"] > 0.02
                              and r["quantum_minus_classical"] > 0.02 and r["beyond_distance_best"] > r["control_best_seed"])
    for name, members in (("val", VALIDATION), ("rep", REPLICATION), ("held", HELDOUT), ("all", PROTEINS)):
        sub = [r for r in rows if r["pdb"] in members]
        ps = sorted(r["p_value"] for r in sub)                   # Benjamini-Hochberg at 5%
        put("agg", name, "nbh", max([k + 1 for k, p in enumerate(ps) if p <= 0.05 * (k + 1) / len(ps)], default=0))
        put("agg", name, "nsurvive", sum(survive(r) for r in sub))
        put("agg", name, "survivors", ", ".join(r["pdb"] for r in sub if survive(r)) or "none")
        put("agg", name, "nquantum", sum(quantum_only(r) for r in sub))
        # distinct proteins: several structures of one protein count once (significant if any is); among
        # held-out sets, the structure that repeats a development protein is left out
        groups = {}
        for r in sub:
            if name != "all" and SAME_PROTEIN.get(r["pdb"]) in DEVELOPMENT:
                continue
            key = SAME_PROTEIN.get(r["pdb"], r["pdb"])
            groups[key] = groups.get(key, False) or r["p_value"] < 0.05
        nd, sd = len(groups), sum(groups.values())
        put("agg", name, "ndistinct", nd)
        put("agg", name, "nsigdistinct", sd)
        put("agg", name, "binpdistinct", fp(float(binom.sf(sd - 1, nd, 0.05))) if nd else "--")
    # the established predictors of external_baselines.py, when it has been run
    ext = os.path.join(HERE, "external.csv")
    if os.path.isfile(ext):
        with open(ext) as f:
            erows = {r["pdb"]: r for r in csv.DictReader(f)}
        for pid, r in erows.items():
            for m in EXTERNAL:
                put("ext", pid, m, f3(float(r[f"{m}_beyond"])))
                put("ext", pid, m + "p", fp(float(r[f"{m}_p"])))
        for name, members in (("val", VALIDATION), ("rep", REPLICATION), ("held", HELDOUT), ("all", PROTEINS)):
            sub = [erows[p] for p in members if p in erows]
            for m in EXTERNAL:
                put("agg", name, "nsig" + m, sum(float(r[f"{m}_p"]) < 0.05 for r in sub))
                put("agg", name, "mean" + m, f3(float(np.mean([float(r[f"{m}_beyond"]) for r in sub]))))
    # one setting at a time for the significant proteins (sensitivity_more.py), when it has been run
    more = os.path.join(HERE, "sensitivity_more.csv")
    if os.path.isfile(more):
        with open(more) as f:
            mrows = list(csv.DictReader(f))
        sig = lambda r: r["p_value"] != "" and float(r["p_value"]) < 0.05     # blank: no allosteric test
        beats_ew = lambda r: (r["beyond_distance_best"] != "" and r["energy_weighted_best"] != ""
                              and float(r["beyond_distance_best"]) > float(r["energy_weighted_best"]) + 0.02)
        put("agg", "more", "n", len(mrows))
        put("agg", "more", "nprot", len({r["pdb"] for r in mrows}))
        put("agg", "more", "nsig", sum(sig(r) for r in mrows))
        put("agg", "more", "nbeatew", sum(beats_ew(r) for r in mrows))
        for pid in sorted({r["pdb"] for r in mrows}):
            sub = [r for r in mrows if r["pdb"] == pid]
            put("more", pid, "nsig", f"{sum(sig(r) for r in sub)} of {len(sub)}")
        for what in ("T", "cutoff", "scale"):
            sub = [r for r in mrows if r["setting"] == what]
            put("agg", "more", "nsig" + what, f"{sum(sig(r) for r in sub)} of {len(sub)}")
    with open(os.path.join(HERE, "numbers.tex"), "w") as f:
        f.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------- figures (greyscale, 300 dpi PNG + SVG)
def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5, "svg.fonttype": "none",
                         "svg.hashsalt": "qubitman",          # same ids every run, so unchanged figures stay unchanged
                         "axes.spines.top": False, "axes.spines.right": False})
    return plt


def _gamma_axis(ax):
    rp.gamma_axis(ax, 100.0)                  # the same axis as the pipeline's own figures
    ax.grid(axis="y", color="#eee")


def _save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=300, metadata={"Software": None})
    fig.savefig(os.path.join(FIG, name + ".svg"), metadata={"Date": None, "Creator": None})


def fig_workflow():
    """The method in one picture: structure -> network -> walks over noise -> scores and tests."""
    plt = _plt()
    fig, ax = plt.subplots(figsize=(12, 2.9))
    ax.set_xlim(0, 12.3)
    ax.set_ylim(-0.15, 2.9)
    ax.axis("off")
    boxes = [("Structure", "PDB id or file;\nknown sites from ALLO"),
             ("Contact network", "residues = nodes;\nCβ within 8 Å = edges"),
             ("Hamiltonian", "H = A + s·diag(ε)\nε: z-scored hydropathy"),
             ("Dephased walk", "ρ from the active site\n18 rates γ, 0 to 100"),
             ("Scores", "signal per residue\n= ∫ ρ_ii dt, t ≤ 30"),
             ("Tests", "AUC beyond distance,\np; classical walks;\nrandom energies")]
    w, gap = 1.75, 0.3
    for k, (title, body) in enumerate(boxes):
        x = 0.1 + k * (w + gap)
        ax.add_patch(plt.Rectangle((x, 0.55), w, 1.8, fc="#f4f4f4" if k % 2 else "#e9e9e9", ec="#222", lw=1.2))
        ax.text(x + w / 2, 2.0, title, ha="center", va="center", fontsize=11, weight="bold")
        ax.text(x + w / 2, 1.2, body, ha="center", va="center", fontsize=9, color="#333", linespacing=1.4)
        if k < len(boxes) - 1:
            ax.annotate("", xy=(x + w + gap - 0.02, 1.45), xytext=(x + w + 0.02, 1.45),
                        arrowprops=dict(arrowstyle="-|>", color="#222", lw=1.3))
    ax.text(6.15, 0.05, "Nulls on every run: 5 random-energy seeds; 10,000 label shuffles within distance shells;\n"
            "two classical walks on the same contacts (plain, and with the energy-dependent rates of the noisy walk)", ha="center", fontsize=9, color="#555")
    fig.tight_layout()
    _save(fig, "workflow")
    plt.close(fig)


def fig_transport():
    """Transport to the distal residues versus noise, with the random-energy band."""
    plt = _plt()
    rows = -(-len(PROTEINS) // 5)
    fig, axs = plt.subplots(rows, 5, figsize=(15, 3.1 * rows + 0.4), sharex=True, squeeze=False)
    for ax in axs.flat[len(PROTEINS):]:
        ax.set_visible(False)
    for ax, pid in zip(axs.flat, PROTEINS):
        r = load(pid, "result")
        g, t = np.array(r["gammas"]), r["transport"]
        m, sd = np.array(t["control"]["mean"]), np.array(t["control"]["sd"])
        _gamma_axis(ax)
        ax.fill_between(g, m - sd, m + sd, color="#8a8a8a", alpha=0.2, lw=0)
        ax.plot(g, m, color="#8a8a8a", ls="--", lw=1.5, label="random site energies (mean ± sd)")
        ax.plot(g, t["distal_mean"], color="#111", lw=2, marker="o", ms=3, mec="white", mew=0.5,
                label="hydropathy site energies")
        st = t["stats"]
        ax.set_title(f"{pid}{MARK.get(SET_OF[pid], '')}: peak at γ = {st['peak_gamma']:g}", fontsize=10.5)
    for ax in axs.flat[len(PROTEINS) - 5:len(PROTEINS)]:
        ax.set_xlabel("dephasing rate γ")
        ax.xaxis.set_tick_params(labelbottom=True)
    for ax in axs[:, 0]:
        ax.set_ylabel("mean signal on distal residues")
    fig.legend(*axs[0, 0].get_legend_handles_labels(), loc="upper center", ncol=2, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.25 / rows))
    _save(fig, "transport")
    plt.close(fig)


def fig_beyond_distance(proteins=None, name="beyond_distance_all", cols=3):
    """AUC beyond distance vs noise: quantum walk, random-energy band, both classical walks (all proteins,
    or a chosen few for a larger figure)."""
    proteins = proteins or PROTEINS
    plt = _plt()
    rows = -(-len(proteins) // cols)
    fig, axs = plt.subplots(rows, cols, figsize=(6 * cols, 3.6 * rows + 0.5), sharex=True, sharey=True, squeeze=False)
    for ax in axs.flat[len(proteins):]:
        ax.set_visible(False)
    for ax, pid in zip(axs.flat, proteins):
        r = load(pid, "result")
        g, d = np.array(r["gammas"]), r["allosteric"]["adjusted"]
        m, sd = np.array(d["control"]["mean"]), np.array(d["control"]["sd"])
        _gamma_axis(ax)
        ax.fill_between(g, m - sd, m + sd, color="#8a8a8a", alpha=0.2, lw=0)
        ax.plot(g, m, color="#8a8a8a", ls="--", lw=1.6, label="random site energies (mean ± sd, 5 seeds)")
        ax.plot(g, d["auc"], color="#111", lw=2.2, marker="o", ms=3.5, mec="white", mew=0.6,
                label="quantum walk, hydropathy site energies")
        ax.plot(g[1:], r["allosteric"]["incoherent"]["auc_adjusted"][1:], color="#555", lw=1.5, ls="-.",
                label="energy-weighted classical walk (same γ, no interference)")
        ax.axhline(d["baselines"]["classical walk, best rate"], color="#444", lw=1.3, ls=(0, (1, 1.5)),
                   label="plain classical random walk (best hopping rate)")
        ax.axhline(0.5, color="#bbb", lw=1)
        ax.text(0.012, 0.506, "closeness alone", fontsize=8.5, color="#777")
        st = d["stats"]
        pv = d["significance"]["p_value"]
        ax.set_title(f"{pid}{MARK.get(SET_OF[pid], '')} ({r['summary']['residues']} residues): best {st['peak_value']:.3f} "
                     f"at γ = {st['peak_gamma']:g}, p {'< 0.001' if pv < 0.001 else f'= {pv:.3f}'}", fontsize=10.5)
        ax.set_ylim(0.2, 1.0)
    for c in range(cols):                                   # the lowest panel in each column gets the x axis
        shown = [a for a in axs[:, c] if a.get_visible()]
        if shown:
            shown[-1].set_xlabel("dephasing rate γ")
            shown[-1].xaxis.set_tick_params(labelbottom=True)
    for ax in axs[:, 0]:
        ax.set_ylabel("AUC beyond distance")
    fig.legend(*axs[0, 0].get_legend_handles_labels(), loc="upper center", ncol=2, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.18 / rows))
    _save(fig, name)
    plt.close(fig)


def fig_methods():
    """Every method's AUC beyond distance on every protein; filled markers are significant (p < 0.05)."""
    ext = os.path.join(HERE, "external.csv")
    if not os.path.isfile(ext):
        return
    plt = _plt()
    with open(ext) as f:
        erows = {r["pdb"]: r for r in csv.DictReader(f)}
    with open(os.path.join(HERE, "proteins.csv")) as f:
        prows = {r["pdb"]: r for r in csv.DictReader(f)}
    methods = [("quantum walk", "o", lambda p: (prows[p]["beyond_distance_best"], prows[p]["p_value"])),
               ("energy-weighted classical", "s", lambda p: (prows[p]["energy_weighted_best"], prows[p]["energy_weighted_p_value"])),
               ("plain classical", "D", lambda p: (prows[p]["classical_best"], prows[p]["classical_p_value"])),
               ("closeness", "^", lambda p: (erows[p]["closeness_beyond"], erows[p]["closeness_p"])),
               ("betweenness", "v", lambda p: (erows[p]["betweenness_beyond"], erows[p]["betweenness_p"])),
               ("perturbation response", "P", lambda p: (erows[p]["prs_beyond"], erows[p]["prs_p"]))]
    fig, ax = plt.subplots(figsize=(10, 0.42 * len(PROTEINS) + 1.6))
    ys = np.arange(len(PROTEINS))[::-1]
    for k, (label, mk, get) in enumerate(methods):
        off = (k - (len(methods) - 1) / 2) * 0.11
        for y, pid in zip(ys, PROTEINS):
            v, pv = (float(x) for x in get(pid))
            ax.scatter(v, y + off, marker=mk, s=34, edgecolors="#111", linewidths=0.9,
                       facecolors="#111" if pv < 0.05 else "white", zorder=3)
    ax.axvline(0.5, color="#bbb", lw=1, zorder=1)
    for y in ys[:-1]:
        ax.axhline(y - 0.5, color="#eee", lw=0.8, zorder=0)
    ax.set_yticks(ys, [f"{p}{MARK.get(SET_OF[p], '')}" for p in PROTEINS])
    ax.set_ylim(-0.6, len(PROTEINS) - 0.4)
    ax.set_xlabel("AUC beyond distance (0.5 = distance alone)")
    from matplotlib.lines import Line2D
    keys = [Line2D([], [], marker=mk, ls="", mfc="white", mec="#111", ms=7, label=label) for label, mk, _ in methods]
    keys.append(Line2D([], [], marker="o", ls="", mfc="#111", mec="#111", ms=7, label="filled: p < 0.05"))
    fig.legend(handles=keys, loc="upper center", ncol=4, frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.75 / (0.42 * len(PROTEINS) + 1.6)))
    _save(fig, "methods")
    plt.close(fig)


def fig_sensitivity():
    """Cutoff x scale grids: best AUC beyond distance (shade), p and the margin over the classical walk."""
    plt = _plt()
    fig, axs = plt.subplots(1, len(SENSITIVITY), figsize=(5.6 * len(SENSITIVITY), 4.8))
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
                sg = lambda x: f"{'+' if round(x, 3) >= 0 else '−'}{abs(x):.3f}"
                ax.text(k, i, f"{v:.3f}\n{'p < 0.001' if p < 0.001 else f'p = {p:.3f}'}\n{sg(float(r['quantum_minus_classical']))} vs cl."
                        f"\n{sg(float(r['quantum_minus_energy_weighted']))} vs e-w",
                        ha="center", va="center", fontsize=8.5, color="white" if (v - 0.45) / 0.5 > 0.62 else "#111")
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
        ax.scatter(*xy[src].T, c="#111" if mode == "signal" else "#b5b5b5", s=34, zorder=3)   # (b): not compared
        ax.scatter(*xy[src].T, facecolors="none", edgecolors="#111", s=150, linewidths=1.2, zorder=3)
        ax.scatter(*xy[allo].T, facecolors="none", edgecolors="#666", s=150, linewidths=1.2,
                   linestyles=(0, (2, 1.5)), zorder=3)
        ax.set_aspect("equal")
        ax.axis("off")
    from matplotlib.lines import Line2D
    keys = [Line2D([], [], marker="o", ls="", mfc="none", mec="#111", ms=11, label="active site (walk starts; grey in (b))"),
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
                ["pdb", "set", "residues", "allosteric_residues", "beyond_distance_best", "best_gamma", "noise_gain",
                 "p_value", "classical_best", "energy_weighted_best", "control_best_seed", "difference_auc", "difference_p"])
    for pid in SENSITIVITY:
        for ext in ("csv", "md"):
            shutil.copy(os.path.join(RUNS, "sensitivity", pid, f"{pid}_sensitivity.{ext}"),
                        os.path.join(HERE, f"sensitivity_{pid}.{ext}"))
    write_numbers(rows)
    for make in (fig_workflow, fig_transport, fig_beyond_distance, fig_sensitivity, fig_map):
        make()
    fig_beyond_distance(KEY_PROTEINS, "beyond_distance", cols=2)
    fig_methods()
    print("\nWrote results/proteins.csv/.md, results/sensitivity_*.csv/.md and results/figures/")


if __name__ == "__main__":
    main()
