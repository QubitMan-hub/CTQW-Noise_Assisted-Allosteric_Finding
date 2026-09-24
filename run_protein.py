#!/usr/bin/env python3
"""run_protein.py - one command, one protein.

    python run_protein.py STRUCTURE            (a PDB/mmCIF file or a 4-character PDB id)

It builds the residue network and sweeps the dephasing rate gamma, writing to
the output folder:

    NAME.graphml             the residue network
    NAME_hump.png/.svg/.csv  transport vs noise: mean signal reaching the distal
                             residues (the far half of the network from the start)
    NAME_allosteric.*        if active/allosteric labels are known: how well the
                             walk from the active site ranks the known allosteric
                             residues (ROC AUC) at each gamma, with baselines
    NAME_ranking.csv         every residue ranked by the signal it receives
    NAME_parameters.json     every setting, the labels used and a code version
    NAME_result.json         all numbers behind the figures

Any PDB id or structure file works. The walk starts at the active site when one
is known: looked up in the bundled ALLO benchmark table (118 proteins) by PDB id,
or given with --active. Known allosteric residues (from ALLO or --allosteric)
add the ROC AUC test; without them you still get the ranked candidate residues.
With no active site at all, the walk starts at --source (default: the
most-connected residue) and gives the transport curve and ranking.
The null model runs by default: the same sweep with random site energies over
several seeds, drawn as a band (--no-control skips it for a quick look).
Sweeps run in parallel over the CPU cores.
"""
import argparse, csv, hashlib, json, os, subprocess, sys, time
import numpy as np, networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rin_builder as rb
import walk_core as wc
import labels as lb

PIPELINE_FILES = ("run_protein.py", "walk_core.py", "rin_builder.py", "labels.py",
                  os.path.join("data", "allo_labels.csv"))

DEFAULTS = dict(
    source=None, chains=None, model=0, method="cb", cutoff=None, min_seq_sep=0, weighted=False,
    site_energy="hydropathy", scale=3.0, seed=0, gammas=list(wc.DEFAULT_GAMMAS),
    tmax=30.0, ntime=200, distal_fraction=0.5,
    control=True, control_seeds=5, permutations=10000,
    labels="auto", active=None, allosteric=None, site=None,
    workers=None,
)

CLASSICAL_RATES = [float(f"{k:.4g}") for k in np.logspace(-2, 2, 17)]   # classical hopping rates tried


class InputError(Exception):
    """A problem with the input, reported as a short message (no traceback)."""


# ---------------------------------------------------------------- building blocks
def structure_to_graph(inp, outdir, model, chains, method, cutoff, min_seq_sep, weighted):
    residues = rb.collect_residues(rb.load_structure(inp, tmpdir=outdir), model, chains)
    if not residues:
        raise ValueError("no amino-acid residues found.")
    parts = (rb.build_by_heavy_atoms(residues, cutoff, min_seq_sep) if method == "heavy"
             else rb.build_by_representative(residues, method, cutoff, min_seq_sep))
    return rb.assemble_graph(*parts, weighted)


def code_version():
    """Git commit (+ whether pipeline files differ from it) and a hash of the pipeline files."""
    h = hashlib.sha256()
    for name in PIPELINE_FILES:
        p = os.path.join(HERE, name)
        if os.path.isfile(p):
            with open(p, "rb") as f:
                h.update(name.encode() + b"\0" + f.read())
    out = {"pipeline_sha256": h.hexdigest()[:16], "git_commit": None, "git_dirty": None}
    try:
        r = subprocess.run(["git", "-C", HERE, "rev-parse", "HEAD"], capture_output=True, text=True, timeout=5)
        if r.returncode == 0:
            out["git_commit"] = r.stdout.strip()
            d = subprocess.run(["git", "-C", HERE, "status", "--porcelain", "--", *PIPELINE_FILES],
                               capture_output=True, text=True, timeout=5)
            out["git_dirty"] = bool(d.stdout.strip())
    except (OSError, subprocess.SubprocessError):
        pass
    return out


def gamma_axis(ax, gmax):
    """The x axis every figure shares: log(1 + gamma/0.05), so gamma = 0 sits at the left
    edge and the low-noise points stay readable; ticks at 0, 0.1, 1, 10, 100."""
    from matplotlib.ticker import FixedLocator, NullLocator
    ax.set_xscale("function", functions=(lambda g: np.log10(1 + np.asarray(g) / 0.05),
                                         lambda u: 0.05 * (10 ** np.asarray(u) - 1)))
    ax.set_xlim(0, gmax * 1.05)
    ticks = [t for t in (0, 0.1, 1, 10, 100) if t <= gmax * 1.05]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xticklabels([f"{t:g}" for t in ticks])


def save_figure(base, gammas, main, main_label, ylabel, control=None, control_label=None,
                baselines=(), chance=None, peak_gamma=None):
    """Greyscale publication figure (300 dpi PNG + SVG) on the shared gamma axis."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 13, "axes.linewidth": 0.8,
                         "svg.fonttype": "none"})
    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    g = np.asarray(gammas, dtype=float)
    lo = [float(np.nanmin(main))]
    ax.plot(g, main, color="#161616", lw=2.0, marker="o", ms=3.6, mfc="#161616", mec="white",
            mew=0.8, label=main_label, zorder=3)
    if control is not None:
        m, s = np.asarray(control[0]), np.asarray(control[1])
        ax.fill_between(g, m - s, m + s, color="#8a8a8a", alpha=0.18, lw=0, zorder=1)
        ax.plot(g, m, color="#8a8a8a", lw=1.8, ls=(0, (5, 3)), label=control_label, zorder=2)
        lo.append(float(np.nanmin(m - s)))
    for value, label in baselines:
        ax.axhline(value, color="#6b6b6b", lw=1, ls=(0, (1, 2)), zorder=1)
    if chance is not None:
        ax.axhline(chance, color="#bdbdbd", lw=0.8, zorder=0)
    if peak_gamma is not None:
        ax.axvline(peak_gamma, color="#8a8a8a", lw=1, ls=(0, (3, 3)), zorder=0)
    top = max([float(np.nanmax(main))] + [float(np.nanmax(np.asarray(control[0]) + np.asarray(control[1])))
                                          if control is not None else 0] + [v for v, _ in baselines])
    bottom = 0.0 if chance is None else min(chance, min(lo), min([v for v, _ in baselines] or [1])) - 0.05
    ax.set_ylim(max(0.0, bottom), top * 1.08 if chance is None else min(1.0, top + 0.05))
    y0, y1 = ax.get_ylim()
    last = None                                   # labels right of the plot, nudged apart when lines are close
    for value, label in sorted(baselines, reverse=True):
        y = value if last is None else min(value, last - 0.055 * (y1 - y0))
        ax.text(g[-1], y, f" {label}", va="center", ha="left", fontsize=10, color="#3a3a3a")
        last = y
    gamma_axis(ax, g[-1])
    ax.set_xlabel("dephasing rate γ")
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", color="#e5e5e5", lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#6b6b6b")
    ax.tick_params(colors="#3a3a3a", length=3)
    if control is not None:
        ax.legend(frameon=False, fontsize=10.5, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=1,
                  handlelength=3.2, borderaxespad=0.3)
    fig.tight_layout()
    fig.savefig(base + ".png", dpi=300)
    fig.savefig(base + ".svg")
    plt.close(fig)


def jsonable(obj):
    """Plain-JSON copy: numpy -> Python, NaN/inf -> None (strict JSON for browsers)."""
    if isinstance(obj, dict):
        return {str(k): jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return jsonable(obj.tolist())
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        return float(obj) if np.isfinite(obj) else None
    return obj


def _control(gammas, curves, main_gain):
    """Random-energy seeds as a mean +- sd band, and where the real hump's gain falls among theirs."""
    c = np.asarray(curves, dtype=float)
    m = c.mean(axis=0)
    cg = np.array([wc.hump_stats(gammas, x)["gain_abs"] for x in c])
    sd = float(cg.std())
    return {"mean": m.tolist(), "sd": c.std(axis=0).tolist(), "stats": wc.hump_stats(gammas, m),
            "control_gains": cg.tolist(), "exceeds_seeds": int((main_gain > cg).sum()), "n_seeds": len(cg),
            "z": float((main_gain - cg.mean()) / sd) if sd > 0 else None}


# ---------------------------------------------------------------- the analysis
def analyze(inp, outdir, prefix, opts=None, log=print):
    """Run the whole pipeline on one structure; write files; return every number."""
    o = dict(DEFAULTS, **(opts or {}))
    t0 = time.time()
    os.makedirs(outdir, exist_ok=True)
    cutoff = o["cutoff"] if o["cutoff"] is not None else (5.0 if o["method"] == "heavy" else 8.0)
    gammas = [float(g) for g in o["gammas"]]
    if gammas[0] != 0 or any(b <= a for a, b in zip(gammas, gammas[1:])):
        raise InputError("gammas must start at 0 and increase.")
    tlist = np.linspace(0, o["tmax"], o["ntime"])
    notes = []

    # 1) labels first, because they decide which chains the network uses
    try:
        labels, label_note, allo_entries = lb.resolve(inp, o["labels"], o["active"], o["allosteric"],
                                                       o["site"], (o["chains"] or "").split(",")[0].strip() or None)
    except (ValueError, OSError) as e:
        raise InputError(f"labels: {e}")
    chains = {c.strip() for c in o["chains"].split(",") if c.strip()} if o["chains"] else None
    if labels and chains is None:
        chains = set(lb.chains_of(labels))
        notes.append(f"network restricted to chain(s) {','.join(sorted(chains))}, the chains the labels "
                     f"are numbered against (pass --chains to override)")
    log(f"[labels] {label_note}")

    # 2) structure -> network
    try:
        G = structure_to_graph(inp, outdir, o["model"], chains, o["method"], cutoff,
                               o["min_seq_sep"], o["weighted"])
    except Exception as e:
        msg = (str(e).strip().splitlines() or [type(e).__name__])[0][:160]
        if "no amino-acid residues" in msg:
            raise InputError("no amino-acid residues found" + (f" in chain(s) {','.join(sorted(chains))}." if chains else "."))
        if isinstance(e, ValueError):                  # the structure was read; say what is wrong with it
            raise InputError(msg)
        if not os.path.isfile(inp):
            raise InputError(f"could not fetch PDB id {inp}; check the id and your internet connection.")
        raise InputError(f"could not read the structure file: {msg}")
    if G.number_of_nodes() < 3:
        raise InputError(f"only {G.number_of_nodes()} residue(s) found; the walk needs a real network.")
    graphml = os.path.join(outdir, prefix + ".graphml")
    nx.write_graphml(G, graphml)
    _, nodes, A, resnames = wc.network_arrays(G)
    idx = {n: i for i, n in enumerate(nodes)}
    N = len(nodes)
    degree = dict(G.degree())
    chains_used = sorted({n.split(":")[0] for n in nodes})
    log(f"[network] {N} residues, {G.number_of_edges()} contacts, chains {chains_used}")

    if o["source"] and o["source"] not in idx and o["source"].upper() in idx:
        o["source"] = o["source"].upper()                  # chain ids are case-sensitive; accept a:151 for A:151
    if o["source"] and o["source"] not in idx:
        raise InputError(f"start residue {o['source']} is not in the network; ids look like {', '.join(nodes[:3])}.")

    # 3) labels onto the network
    allo = None
    if labels:
        rep = lb.match(labels, nodes, resnames)
        for k in ("active_missing", "allosteric_missing"):
            if rep[k]:
                notes.append(f"{len(rep[k])} {k.split('_')[0]} residue(s) not in the network: {', '.join(rep[k])}")
        if rep["name_mismatch"]:
            notes.append("residue names differ from the labels (possible numbering offset): "
                         + "; ".join(rep["name_mismatch"]))
        eligible = np.ones(N, dtype=bool)
        eligible[rep["active"]] = False
        positive = np.zeros(N, dtype=bool)
        positive[rep["allosteric"]] = True
        positive &= eligible
        if not rep["active"]:
            notes.append("none of the active-site residues are in the network; walking from one residue instead")
        elif positive.any() and (eligible & ~positive).any():
            allo = {"rep": rep, "eligible": eligible, "positive": positive}
        elif labels["allosteric"]:
            notes.append("allosteric test skipped: none of the known allosteric residues are in the network")

    # 4) one walk per noise level: from the active site when it is known (feeding
    #    the transport curve, the ranking and the allosteric test), else from one residue
    active = rep["active"] if labels else []
    if active:
        walk_src = list(active)
        walk_from = "active site"
    else:
        source_id = o["source"] or max(degree, key=degree.get)
        walk_src = [idx[source_id]]
        walk_from = source_id
    hops = wc.hop_distance(G, nodes, [nodes[i] for i in walk_src])        # contacts from the start
    distal, distal_hops = wc.distal_from_hops(hops, o["distal_fraction"])
    if not distal.any():
        raise InputError("the network is too small: no residues lie far from the walk's start.")

    # 5) every walk this run needs, in one parallel batch
    hams = {"main": wc.build_hamiltonian(A, resnames, o["site_energy"], o["scale"], o["seed"])}
    control_keys = []
    if o["control"]:
        if o["site_energy"] == "random":
            notes.append("control skipped: the main run already uses random site energies")
        else:
            for k in range(int(o["control_seeds"])):
                hams[f"rand{k}"] = wc.build_hamiltonian(A, resnames, "random", o["scale"], o["seed"] + k)
                control_keys.append(f"rand{k}")
    jobs = [(key, walk_src, g) for key in hams for g in gammas]
    log(f"[sweep] {len(jobs)} walks over {len(gammas)} dephasing rates"
        f"{' with ' + str(len(control_keys)) + ' random-energy seeds' if control_keys else ''}")
    flat = wc.sweep_many(hams, jobs, tlist, o["workers"])
    scores = {key: np.array(flat[i * len(gammas):(i + 1) * len(gammas)]) for i, key in enumerate(hams)}

    cwalk = wc.classical_walk(A, walk_src, o["tmax"])     # classical baseline: one eigendecomposition, any rate

    # 6) transport metrics
    S = scores["main"]
    transport = S[:, distal].mean(axis=1)
    t_stats = wc.hump_stats(gammas, transport)
    t_ctrl = None
    if control_keys:
        curves = [scores[k][:, distal].mean(axis=1) for k in control_keys]
        t_ctrl = dict(_control(gammas, curves, t_stats["gain_abs"]),
                      seed_verdicts=[wc.hump_stats(gammas, c)["verdict"] for c in curves])

    # quantum vs classical at every noise level: the classical hopping rate is set so both walks send
    # the same mean signal to the distal residues, then they are compared residue by residue
    rates = [wc.matched_rate(cwalk, distal, transport[j]) for j in range(len(gammas))]
    CM = np.array([cwalk(k) for k in rates])
    off = [f"{g:g}" for g, c, t in zip(gammas, CM, transport) if abs(c[distal].mean() - t) > 1e-3 * t]
    if off:                                        # outside the rates the matching searches (1e-3 to 1e3)
        notes.append(f"the classical walk could not match the quantum walk's spread at gamma {', '.join(off)}; "
                     "its difference map there compares walks that spread differently")

    # 7) allosteric metrics
    a_res = None
    if allo:
        el, pos = allo["eligible"], allo["positive"]
        reach = el & (hops >= 0)                     # candidates connected to the active site
        if not (pos & reach).any() or not (~pos & reach).any():
            notes.append("allosteric test skipped: no known allosteric residue is connected to the active site "
                         "in this network (try other chains or a larger cutoff)")
            allo = None
    if allo:
        auc = np.array([wc.roc_auc(s[el], pos[el]) for s in S])
        a_stats = wc.hump_stats(gammas, auc)
        dist = hops                                  # the walk starts at the active site here
        prox = np.where(dist >= 0, -dist, -N).astype(float)
        deg = np.array([degree[n] for n in nodes], dtype=float)
        # distance-adjusted: each residue compared only with residues equally far from the
        # active site, so the test measures what the walk adds beyond plain proximity
        shells = wc.distance_shells(dist, reach)
        adj_auc = lambda sc: wc.roc_auc(wc.shell_percentile(sc, shells, N)[reach], pos[reach])
        auc_adj = np.array([adj_auc(s) for s in S])
        adj_stats = wc.hump_stats(gammas, auc_adj)
        baselines = {"contact degree": wc.roc_auc(deg[el], pos[el]),
                     "proximity to active site": wc.roc_auc(prox[el], pos[el])}
        adj_baselines = {"contact degree, same distance": adj_auc(deg)}
        # classical baseline: a random walk on the same contacts from the same start,
        # over a wide range of hopping rates; its best rate is the number to beat
        C = np.array([cwalk(k) for k in CLASSICAL_RATES])
        c_raw = np.array([wc.roc_auc(c[el], pos[el]) for c in C])
        c_adj = np.array([adj_auc(c) for c in C])
        baselines["classical walk, best rate"] = float(np.nanmax(c_raw))
        adj_baselines["classical walk, best rate"] = float(np.nanmax(c_adj))
        # energy-weighted classical walk: the hopping rates the dephased walk reduces to (same
        # contacts, energies and gamma, no interference); its best gamma is the stricter baseline
        IC = wc.incoherent_scores(hams["main"], walk_src, gammas, o["tmax"])
        live = [j for j, g in enumerate(gammas) if g > 0]
        e_raw, e_adj = np.full(len(gammas), np.nan), np.full(len(gammas), np.nan)
        e_raw[live] = [wc.roc_auc(IC[j][el], pos[el]) for j in live]
        e_adj[live] = [adj_auc(IC[j]) for j in live]
        baselines["energy-weighted classical walk, best γ"] = float(np.nanmax(e_raw))
        adj_baselines["energy-weighted classical walk, best γ"] = float(np.nanmax(e_adj))
        # significance: shuffle the allosteric labels within each distance shell
        sig = wc.shell_permutation_test(S, shells, pos & reach, n_perm=o["permutations"], seed=o["seed"])
        c_sig = wc.shell_permutation_test(C, shells, pos & reach, n_perm=o["permutations"], seed=o["seed"])
        e_sig = wc.shell_permutation_test(IC[live], shells, pos & reach, n_perm=o["permutations"], seed=o["seed"])
        classical = {"rates": list(CLASSICAL_RATES), "auc": c_raw.tolist(), "auc_adjusted": c_adj.tolist(),
                     "best_rate_raw": float(CLASSICAL_RATES[int(np.nanargmax(c_raw))]),
                     "best_rate_adjusted": float(CLASSICAL_RATES[int(np.nanargmax(c_adj))]),
                     "significance_adjusted": c_sig}
        incoherent = {"auc": e_raw.tolist(), "auc_adjusted": e_adj.tolist(),
                      "best_gamma_adjusted": float(gammas[int(np.nanargmax(e_adj))]),
                      "significance_adjusted": e_sig}
        a_ctrl = adj_ctrl = None
        if control_keys:
            a_ctrl = _control(gammas, [[wc.roc_auc(s[el], pos[el]) for s in scores[k]] for k in control_keys],
                              a_stats["gain_abs"])
            curves = [[adj_auc(s) for s in scores[k]] for k in control_keys]
            adj_ctrl = dict(_control(gammas, curves, adj_stats["gain_abs"]),
                            best_seeds=[float(np.nanmax(c)) for c in curves])
        # does quantum-minus-classical signal itself pick out the allosteric residues, at the gamma
        # where the AUC beyond distance peaks? (same shells and permutation test)
        j = adj_stats["peak_index"]
        dq = S[j] - CM[j]
        d_sig = wc.shell_permutation_test([dq], shells, pos & reach, n_perm=o["permutations"], seed=o["seed"])
        difference = {"gamma": gammas[j], "classical_rate": rates[j], "auc": adj_auc(dq),
                      "p_value": d_sig["p_value"] if d_sig else None,
                      "allosteric_quantum_favoured": float((dq[pos] > 0).mean())}
        rep = allo["rep"]
        a_res = {"auc": auc.tolist(), "stats": a_stats, "baselines": baselines, "control": a_ctrl,
                 "adjusted": {"auc": auc_adj.tolist(), "stats": adj_stats, "baselines": adj_baselines,
                              "control": adj_ctrl, "shells": len(shells), "significance": sig,
                              "difference": difference},
                 "classical": classical, "incoherent": incoherent,
                 "n_active": len(rep["active"]), "n_allosteric": int(pos.sum()),
                 "n_candidates": int(el.sum()),
                 "active": [nodes[i] for i in rep["active"]], "allosteric": [nodes[i] for i in np.where(pos)[0]],
                 "missing": {"active": rep["active_missing"], "allosteric": rep["allosteric_missing"]},
                 "name_mismatch": rep["name_mismatch"]}

    # 8) signal map + ranking where the AUC beyond distance peaks (labels) or at the transport peak
    peak = a_res["adjusted"]["stats"]["peak_index"] if a_res else t_stats["peak_index"]
    map_scores, map_gamma, map_from = S[peak], gammas[peak], walk_from
    sources_idx = set(walk_src)
    known = set(np.where(allo["positive"])[0]) if allo else set()
    order = [i for i in np.argsort(-map_scores, kind="stable") if i not in sources_idx]
    ranking = [{"rank": r + 1, "id": nodes[i], "resname": resnames[i], "score": float(map_scores[i]),
                "degree": int(degree[nodes[i]]), "known_allosteric": i in known} for r, i in enumerate(order)]

    k_match, map_classical = rates[peak], CM[peak]
    map_diff = map_scores - map_classical
    qvc = {"gamma": map_gamma, "classical_rate": k_match, "distal_quantum": float(transport[peak]),
           "distal_classical": float(map_classical[distal].mean()), "classical_rates": rates}

    # 9) files
    files = {"graphml": prefix + ".graphml", "ranking_csv": prefix + "_ranking.csv",
             "parameters": prefix + "_parameters.json", "result": prefix + "_result.json",
             "qvc_csv": prefix + "_quantum_vs_classical.csv"}
    # (file tag, csv column, values, stats, control, y label, baselines, chance line)
    curves = [("hump", "transport_distal_mean", transport, t_stats, t_ctrl, "mean signal reaching distal residues", {}, None)]
    if a_res:
        d = a_res["adjusted"]
        curves += [("allosteric", "auc", a_res["auc"], a_stats, a_res["control"], "ROC AUC, known allosteric residues",
                    a_res["baselines"], 0.5),
                   ("adjusted", "auc", d["auc"], d["stats"], d["control"],
                    "ROC AUC among residues equally far from the active site", d["baselines"], 0.5)]
    for tag, column, values, st, ctrl, ylabel, bl, chance in curves:
        files.update({f"{tag}_{ext}": f"{prefix}_{tag}.{ext}" for ext in ("csv", "png", "svg")})
        with open(os.path.join(outdir, files[f"{tag}_csv"]), "w") as f:
            f.write(f"gamma,{column}" + (",control_mean,control_sd" if ctrl else "") + "\n")
            for i, g in enumerate(gammas):
                f.write(f"{g},{values[i]:.6f}" + (f",{ctrl['mean'][i]:.6f},{ctrl['sd'][i]:.6f}" if ctrl else "") + "\n")
        save_figure(os.path.join(outdir, f"{prefix}_{tag}"), gammas, values, f"{o['site_energy']} site energies", ylabel,
                    control=(ctrl["mean"], ctrl["sd"]) if ctrl else None,
                    control_label=f"random site energies ({len(control_keys)} seeds, mean ± sd)",
                    baselines=[(v, k) for k, v in bl.items()], chance=chance, peak_gamma=st["peak_gamma"])
    if a_res:
        files["energy_weighted_csv"] = prefix + "_energy_weighted.csv"
        with open(os.path.join(outdir, files["energy_weighted_csv"]), "w") as f:
            f.write("gamma,auc,auc_beyond_distance\n")
            for g, r_, d_ in zip(gammas, a_res["incoherent"]["auc"], a_res["incoherent"]["auc_adjusted"]):
                f.write(f"{g},{r_:.6f},{d_:.6f}\n")
    with open(os.path.join(outdir, files["ranking_csv"]), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "residue", "resname", f"signal_from_{map_from.replace(' ', '_')}_at_gamma_{map_gamma:g}",
                    "contacts", "known_allosteric"])
        w.writerows([[r["rank"], r["id"], r["resname"], f"{r['score']:.6f}", r["degree"],
                      "yes" if r["known_allosteric"] else ""] for r in ranking])
    with open(os.path.join(outdir, files["qvc_csv"]), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "residue", "resname", "hops_from_start", f"quantum_at_gamma_{map_gamma:g}",
                    f"classical_at_rate_{k_match:.4g}", "quantum_minus_classical", "known_allosteric"])
        w.writerows([[r + 1, nodes[i], resnames[i], int(hops[i]), f"{map_scores[i]:.6f}", f"{map_classical[i]:.6f}",
                      f"{map_diff[i]:.6f}", "yes" if i in known else ""]
                     for r, i in enumerate(i for i in np.argsort(-map_diff, kind="stable") if i not in sources_idx)])

    params = {k: o[k] for k in DEFAULTS if k != "workers"}
    params.update(input=inp if not os.path.isfile(inp) else os.path.basename(inp), cutoff=cutoff, gammas=gammas,
                  chains_used=chains_used,
                  walk_from=walk_from, distal_min_hops=distal_hops, labels_used=(
                      {k: v for k, v in labels.items() if k not in ("active", "allosteric")} if labels else None),
                  code_version=code_version(), created_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    if os.path.isfile(inp):
        with open(inp, "rb") as f:
            params["input_sha256"] = hashlib.sha256(f.read()).hexdigest()
    with open(os.path.join(outdir, files["parameters"]), "w") as f:
        json.dump(jsonable(params), f, indent=2)

    xy = _layout(G, nodes)
    result = {
        "name": prefix, "elapsed_s": round(time.time() - t0, 1), "notes": notes, "label_note": label_note,
        "summary": {"residues": N, "contacts": G.number_of_edges(),
                    "chains": chains_used, "walk_from": walk_from,
                    "n_start": len(walk_src),
                    "distal_count": int(distal.sum()), "distal_min_hops": distal_hops,
                    "site_energy": o["site_energy"], "scale": o["scale"]},
        "gammas": gammas,
        "transport": {"distal_mean": transport.tolist(),
                      "stats": t_stats, "control": t_ctrl},
        "allosteric": a_res,
        "allo_entries": allo_entries,
        "map": {"nodes": [{"id": n, "resname": resnames[i], "x": float(xy[i, 0]), "y": float(xy[i, 1]),
                           "score": float(map_scores[i]), "classical": float(map_classical[i]),
                           "diff": float(map_diff[i]), "degree": int(degree[n]),
                           "role": ("source" if i in sources_idx else "allosteric" if i in known else "")}
                          for i, n in enumerate(nodes)],
                "edges": [[idx[u], idx[v]] for u, v in G.edges()], "gamma": map_gamma, "from": map_from,
                # both walks at every gamma, so the page can slide through the noise levels
                "index": peak, "gammas": gammas, "rates": rates,
                "signal": np.round(S, 7).tolist(), "classical": np.round(CM, 7).tolist(),
                "aspect": float(max(xy[:, 1].max(), 1e-3) / max(xy[:, 0].max(), 1e-3))},
        "quantum_vs_classical": qvc, "ranking": ranking, "files": files, "parameters": params,
    }
    with open(os.path.join(outdir, files["result"]), "w") as f:
        json.dump(jsonable({k: v for k, v in result.items() if k != "map"}), f, indent=2)
    return result


def _layout(G, nodes):
    """Residue positions projected onto their two principal axes (spring layout
    if coordinates are missing), scaled to [0, 1] with the aspect ratio kept."""
    try:
        xyz = np.array([[float(G.nodes[n][k]) for k in ("x", "y", "z")] for n in nodes])
        missing = ~xyz.any(axis=1)                  # the network stores (0, 0, 0) when a residue has no C-alpha
        if missing.all():
            raise ValueError("no coordinates")
        idx = {n: i for i, n in enumerate(nodes)}
        for i in np.where(missing)[0]:             # place it among its contacts instead of at the origin
            near = [idx[m] for m in G.neighbors(nodes[i]) if not missing[idx[m]]]
            xyz[i] = xyz[near].mean(axis=0) if near else xyz[~missing].mean(axis=0)
        xyz -= xyz.mean(axis=0)
        _, _, vt = np.linalg.svd(xyz, full_matrices=False)
        xy = xyz @ vt[:2].T
    except (KeyError, ValueError, TypeError):
        pos = nx.spring_layout(G, seed=1)
        xy = np.array([pos[n] for n in nodes])
    xy = xy - xy.min(axis=0)
    return xy / (float(xy.max()) or 1.0)


# ---------------------------------------------------------------- command line
def _verdict_line(stats, what):
    if stats["verdict"] == "hump":
        lo, hi = stats["helps_range"]
        return (f"HUMP in {what}: peak at gamma={stats['peak_gamma']:g} ({stats['peak_value']:.4f}), "
                f"+{stats['gain_abs']:.4f} ({100 * stats['gain_rel']:.1f}%) over the better end; "
                f"noise helps for gamma {lo:g}-{hi:g}")
    return f"No hump in {what}: best {'fully quantum (gamma=0)' if stats['verdict'] == 'quantum' else 'at the most dephased end'}"


def main():
    ap = argparse.ArgumentParser(description="One protein: structure -> network, noise sweep, allosteric test.")
    ap.add_argument("input", help="PDB/mmCIF file, or a 4-character PDB id.")
    ap.add_argument("--source", default=None, help="Start residue when no active site is known, e.g. A:151 (default: most connected).")
    ap.add_argument("--chains", default=None, help="Chains to include (default: the labels' chains, else all).")
    ap.add_argument("--model", type=int, default=0, help="Model index for NMR files.")
    ap.add_argument("--method", choices=["cb", "ca", "heavy"], default="cb")
    ap.add_argument("--cutoff", type=float, default=None)
    ap.add_argument("--min-seq-sep", type=int, default=0)
    ap.add_argument("--weighted", action="store_true")
    ap.add_argument("--site-energy", choices=["none", "hydropathy", "degree", "random"], default="hydropathy")
    ap.add_argument("--scale", type=float, default=3.0, help="Site-energy disorder strength.")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--gammas", default=None, help="Comma-separated dephasing rates starting at 0 (default: 0 + 4 per decade from 0.01 to 100).")
    ap.add_argument("--tmax", type=float, default=30.0)
    ap.add_argument("--ntime", type=int, default=200)
    ap.add_argument("--control", action=argparse.BooleanOptionalAction, default=True,
                    help="Random site energies (null model) over several seeds; on by default, --no-control "
                         "skips it (about 4 times faster, but a hump without it is not evidence).")
    ap.add_argument("--control-seeds", type=int, default=5, help="Random-energy seeds for the control.")
    ap.add_argument("--permutations", type=int, default=10000, help="Label shuffles for the p-value.")
    ap.add_argument("--labels", choices=["auto", "none"], default="auto", help="auto: ALLO table by PDB id; none: skip.")
    ap.add_argument("--active", default=None, help="Active-site residues, e.g. A:57,A:102 (the walk starts here).")
    ap.add_argument("--allosteric", default=None, help="Known allosteric residues, e.g. A:196,A:203 (adds the AUC test).")
    ap.add_argument("--site", default=None, help="Which ALLO entry when a PDB has several (e.g. 2 for 1CE8_2).")
    ap.add_argument("--workers", type=int, default=None, help="Parallel processes (default: CPU cores, max 8).")
    ap.add_argument("--outdir", default="output")
    ap.add_argument("--prefix", default=None)
    a = ap.parse_args()

    prefix = a.prefix or (os.path.splitext(os.path.basename(a.input))[0] or "protein")
    opts = {k: getattr(a, k) for k in DEFAULTS if k != "gammas" and hasattr(a, k)}
    if a.gammas:
        opts["gammas"] = [float(x) for x in a.gammas.split(",")]
    try:
        r = analyze(a.input, a.outdir, prefix, opts)
    except InputError as e:
        sys.exit(f"\n[error] {e}")

    s, t = r["summary"], r["transport"]
    print(f"\n{r['name']}: {s['residues']} residues, walk from {s['walk_from']}, "
          f"{s['distal_count']} distal residues (>= {s['distal_min_hops']} contacts away)")
    print("  " + _verdict_line(t["stats"], "transport"))
    if t["control"]:
        c = t["control"]
        print(f"  random-energy control: {c['stats']['verdict']} (mean curve), hump in "
              f"{c['seed_verdicts'].count('hump')}/{len(c['seed_verdicts'])} seeds; real gain exceeds "
              f"{c['exceeds_seeds']}/{c['n_seeds']} seed gains")
    if r["allosteric"]:
        al = r["allosteric"]
        st = al["stats"]
        print(f"  allosteric test ({al['n_active']} active-site sources, {al['n_allosteric']} known allosteric "
              f"of {al['n_candidates']} candidates): AUC {st['quantum_value']:.3f} quantum, "
              f"{st['peak_value']:.3f} best at gamma={st['peak_gamma']:g}, {st['dephased_value']:.3f} dephased")
        print("  " + _verdict_line(st, "allosteric AUC"))
        print("  baselines: " + ", ".join(f"{k} {v:.3f}" for k, v in al["baselines"].items()))
        ad = al["adjusted"]
        print(f"  beyond distance (same-distance shells, proximity = 0.500): AUC {ad['stats']['quantum_value']:.3f} "
              f"quantum, {ad['stats']['peak_value']:.3f} best at gamma={ad['stats']['peak_gamma']:g}; baseline "
              + ", ".join(f"{k} {v:.3f}" for k, v in ad["baselines"].items())
              + (f"; control best {max(ad['control']['best_seeds']):.3f} over {len(ad['control']['best_seeds'])} seeds"
                 if ad["control"] else ""))
        sg, cs = ad["significance"], al["classical"]["significance_adjusted"]
        if sg:
            print(f"  significance (labels shuffled within distance shells, best over all gammas): "
                  f"p = {sg['p_value']:.4f}; classical walk best {cs['observed_best']:.3f}, p = {cs['p_value']:.4f}")
            es = al["incoherent"]["significance_adjusted"]
            if es:
                print(f"  energy-weighted classical walk (no interference): best {es['observed_best']:.3f} "
                      f"at gamma {al['incoherent']['best_gamma_adjusted']:g}, p = {es['p_value']:.4f}")
    for n in r["notes"]:
        print(f"  [note] {n}")
    print(f"\nDone in {r['elapsed_s']} s. Outputs in {a.outdir}/")


if __name__ == "__main__":
    main()
