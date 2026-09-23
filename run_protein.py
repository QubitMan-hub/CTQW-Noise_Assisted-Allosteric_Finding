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
    NAME.qmod                the Classiq circuit (noise-free walk) + Trotter check
    NAME_parameters.json     every setting, the labels used and a code version
    NAME_result.json         all numbers behind the figures

Labels are looked up automatically in the bundled ALLO benchmark table (118
proteins) by PDB id, or given with --active/--allosteric. With labels, the walk
starts at the active site and that one walk gives both curves; without them it
starts at --source (default: the most-connected residue).
--control adds the null model: the same sweep with random site energies over
several seeds, drawn as a band. Sweeps run in parallel over the CPU cores.
"""
import argparse, hashlib, json, os, subprocess, sys, time
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
    control=False, control_seeds=5,
    labels="auto", active=None, allosteric=None, site=None,
    qmod=True, qmod_time=1.0, qmod_order=2, qmod_repetitions=4, qmod_precision=12,
    max_residues=200, force=False, workers=None,
)

EXEC_SCAFFOLD = '''\
# Run this .qmod on Classiq after classiq.authenticate().
# Noise-free (quantum end):
#   qprog = synthesize(create_model(main)); execute(qprog).result()
# Environment-assisted (the hump): run on SIMULATOR_DENSITY_MATRIX with a
# dephasing channel whose rate is your gamma dial. Confirm the exact
# ClassiqSimulatorNoiseSpecification fields on the live "custom noise models"
# docs page, then sweep the rate and check it matches the hump graph here.
# Coefficients are written with 12 decimals; see NAME_result.json for the
# Trotter fidelity of this circuit against the exact evolution.
'''


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


def write_qmod(terms, n_qubits, source_index, evolution_time, order, repetitions, out_path,
               decimal_precision=12):
    """Emit the Qmod (offline, no login). Source is amplitude-encoded with qubit
    0 as the least significant bit, matching Classiq's convention. Coefficients
    are written with 12 decimals (Classiq's default of 4 perturbs H by ~1e-3)."""
    from pathlib import Path
    from classiq import (Output, QArray, QBit, qfunc, allocate, X, suzuki_trotter,
                         create_model, Pauli, PauliTerm, write_qmod as classiq_write_qmod)
    pmap = {"I": Pauli.I, "X": Pauli.X, "Y": Pauli.Y, "Z": Pauli.Z}
    ham = [PauliTerm(pauli=[pmap[c] for c in s], coefficient=v) for s, v in terms]
    set_bits = [q for q in range(n_qubits) if (source_index >> q) & 1]

    @qfunc
    def main(q: Output[QArray[QBit]]):
        allocate(n_qubits, q)
        for b in set_bits:
            X(q[b])
        suzuki_trotter(ham, evolution_coefficient=evolution_time,
                       order=order, repetitions=repetitions, qbv=q)

    base = os.path.splitext(out_path)[0]
    classiq_write_qmod(create_model(main), os.path.basename(base),
                       directory=Path(os.path.dirname(base) or "."), decimal_precision=decimal_precision)
    return base + ".qmod", len(ham)


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


def save_figure(base, gammas, main, main_label, ylabel, control=None, control_label=None,
                baselines=(), chance=None, peak_gamma=None):
    """Greyscale publication figure (300 dpi PNG + SVG). x is log(1 + gamma/0.05)
    so gamma = 0 sits at the left edge and the low-noise points stay readable."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, NullLocator

    fwd = lambda g: np.log10(1 + np.asarray(g) / 0.05)
    inv = lambda u: 0.05 * (10 ** np.asarray(u) - 1)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 13, "axes.linewidth": 0.8,
                         "svg.fonttype": "none"})
    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    ax.set_xscale("function", functions=(fwd, inv))
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
        ax.text(g[-1], value, f" {label}", va="center", ha="left", fontsize=10, color="#3a3a3a")
    if chance is not None:
        ax.axhline(chance, color="#bdbdbd", lw=0.8, zorder=0)
    if peak_gamma is not None:
        ax.axvline(peak_gamma, color="#8a8a8a", lw=1, ls=(0, (3, 3)), zorder=0)
    ax.set_xlim(0, g[-1] * 1.05)
    top = max([float(np.nanmax(main))] + [float(np.nanmax(np.asarray(control[0]) + np.asarray(control[1])))
                                          if control is not None else 0] + [v for v, _ in baselines])
    bottom = 0.0 if chance is None else min(chance, min(lo), min([v for v, _ in baselines] or [1])) - 0.05
    ax.set_ylim(max(0.0, bottom), top * 1.08 if chance is None else min(1.0, top + 0.05))
    ticks = [t for t in (0, 0.1, 1, 10, 100) if t <= g[-1] * 1.05]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xticklabels([f"{t:g}" for t in ticks])
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


def _band(curves):
    c = np.asarray(curves, dtype=float)
    return c.mean(axis=0), c.std(axis=0)


def _gain_rank(main_gain, control_gains):
    """Where the real hump's gain falls among the random-energy seeds."""
    cg = np.asarray(control_gains, dtype=float)
    sd = float(cg.std())
    return {"control_gains": [float(x) for x in cg], "exceeds_seeds": int((main_gain > cg).sum()),
            "n_seeds": int(len(cg)), "z": (float((main_gain - cg.mean()) / sd) if sd > 0 else None)}


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
                                                       o["site"], (o["chains"] or "").split(",")[0] or None)
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
        if not os.path.isfile(inp):
            raise InputError(f"could not fetch PDB id {inp}; check the id and your internet connection.")
        raise InputError(f"could not read the structure file: {msg}")
    if G.number_of_nodes() < 3:
        raise InputError(f"only {G.number_of_nodes()} residue(s) found; the walk needs a real network.")
    graphml = os.path.join(outdir, prefix + ".graphml")
    nx.write_graphml(G, graphml)
    _, nodes, A, resnames = wc.load_network(graphml)
    idx = {n: i for i, n in enumerate(nodes)}
    N = len(nodes)
    degree = dict(G.degree())
    log(f"[network] {N} residues, {G.number_of_edges()} contacts, chains {sorted({n.split(':')[0] for n in nodes})}")

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
        if rep["active"] and positive.any() and (eligible & ~positive).any():
            allo = {"rep": rep, "eligible": eligible, "positive": positive}
        else:
            notes.append("allosteric test skipped: no usable active-site or allosteric residues in this network")

    # 4) one walk per noise level: from the active site when labels exist (it feeds
    #    both the transport curve and the allosteric test), else from one residue
    if allo:
        walk_src = list(allo["rep"]["active"])
        source_id = o["source"] or max((nodes[i] for i in walk_src), key=degree.get)
        walk_from = "active site"
    else:
        source_id = o["source"] or max(degree, key=degree.get)
        walk_src = [idx[source_id]]
        walk_from = source_id
    distal, distal_hops = wc.distal_mask(G, nodes, [nodes[i] for i in walk_src], o["distal_fraction"])
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
    jobs, where = [], {}
    def add(tag, key, src):
        where[(tag, key)] = len(jobs)
        jobs.extend((key, src, g) for g in gammas)
    for key in ["main"] + control_keys:
        add("walk", key, walk_src)
    log(f"[sweep] {len(jobs)} walks over {len(gammas)} dephasing rates"
        f"{' with ' + str(len(control_keys)) + ' random-energy seeds' if control_keys else ''}")
    flat = wc.sweep_many(hams, jobs, tlist, o["workers"])
    block = lambda tag, key: np.array(flat[where[(tag, key)]: where[(tag, key)] + len(gammas)])

    # 6) transport metrics
    S = block("walk", "main")
    transport = S[:, distal].mean(axis=1)
    t_stats = wc.hump_stats(gammas, transport)
    t_ctrl = None
    if control_keys:
        curves = [block("walk", k)[:, distal].mean(axis=1) for k in control_keys]
        m, s = _band(curves)
        t_ctrl = {"mean": m.tolist(), "sd": s.tolist(), "stats": wc.hump_stats(gammas, m),
                  "seed_verdicts": [wc.hump_stats(gammas, c)["verdict"] for c in curves],
                  **_gain_rank(t_stats["gain_abs"], [wc.hump_stats(gammas, c)["gain_abs"] for c in curves])}

    # 7) allosteric metrics
    a_res = None
    if allo:
        el, pos = allo["eligible"], allo["positive"]
        auc = np.array([wc.roc_auc(s[el], pos[el]) for s in S])
        a_stats = wc.hump_stats(gammas, auc)
        dist_act = {}
        for a in allo["rep"]["active"]:
            for k, d in nx.single_source_shortest_path_length(G, nodes[a]).items():
                dist_act[k] = min(d, dist_act.get(k, d))
        prox = np.array([-dist_act.get(n, N) for n in nodes], dtype=float)
        deg = np.array([degree[n] for n in nodes], dtype=float)
        baselines = {"contact degree": wc.roc_auc(deg[el], pos[el]),
                     "proximity to active site": wc.roc_auc(prox[el], pos[el])}
        a_ctrl = None
        if control_keys:
            curves = [np.array([wc.roc_auc(s[el], pos[el]) for s in block("walk", k)]) for k in control_keys]
            m, s = _band(curves)
            a_ctrl = {"mean": m.tolist(), "sd": s.tolist(), "stats": wc.hump_stats(gammas, m),
                      **_gain_rank(a_stats["gain_abs"], [wc.hump_stats(gammas, c)["gain_abs"] for c in curves])}
        rep = allo["rep"]
        a_res = {"auc": auc.tolist(), "stats": a_stats, "baselines": baselines, "control": a_ctrl,
                 "n_active": len(rep["active"]), "n_allosteric": int(pos.sum()),
                 "n_candidates": int(el.sum()),
                 "active": [nodes[i] for i in rep["active"]], "allosteric": [nodes[i] for i in np.where(pos)[0]],
                 "missing": {"active": rep["active_missing"], "allosteric": rep["allosteric_missing"]},
                 "name_mismatch": rep["name_mismatch"]}

    # 6) signal map + ranking at the best-AUC gamma (labels) or the transport peak
    peak = a_res["stats"]["peak_index"] if a_res else t_stats["peak_index"]
    map_scores, map_gamma, map_from = S[peak], gammas[peak], ("active site" if allo else "source")
    sources_idx = set(walk_src)
    known = set(np.where(allo["positive"])[0]) if allo else set()
    order = [i for i in np.argsort(-map_scores, kind="stable") if i not in sources_idx]
    ranking = [{"rank": r + 1, "id": nodes[i], "resname": resnames[i], "score": float(map_scores[i]),
                "degree": int(degree[nodes[i]]), "known_allosteric": i in known} for r, i in enumerate(order)]

    # 7) Qmod + Trotter check
    qmod = {"status": None, "file": None}
    H = hams["main"]
    if not o["qmod"]:
        qmod["status"] = "skipped_by_user"
    elif N > o["max_residues"] and not o["force"]:
        qmod["status"] = "skipped_too_large"
    else:
        try:
            import classiq  # noqa: F401
        except ImportError:
            qmod["status"] = "unavailable"
        if qmod["status"] is None:
            Hp, nq = wc.pad_to_power_of_two(H)
            terms, n = wc.pauli_decompose(Hp)
            err = float(np.max(np.abs(wc.pauli_reconstruct(terms, n) - Hp)))
            trot = wc.trotter_check(terms, nq, Hp, idx[source_id], o["qmod_time"], o["qmod_order"],
                                    o["qmod_repetitions"], N)
            try:
                qfile, n_terms = write_qmod(terms, nq, idx[source_id], o["qmod_time"], o["qmod_order"],
                                            o["qmod_repetitions"], os.path.join(outdir, prefix + ".qmod"),
                                            o["qmod_precision"])
                with open(os.path.join(outdir, prefix + "_execute.txt"), "w") as f:
                    f.write(EXEC_SCAFFOLD)
                qmod.update(status="written", file=os.path.basename(qfile), pauli_terms=n_terms, qubits=nq,
                            decomposition_error=err, size_kb=round(os.path.getsize(qfile) / 1024),
                            trotter=trot, source_basis_index=idx[source_id])
            except Exception as e:
                qmod.update(status="failed", error=(str(e).splitlines() or [type(e).__name__])[0][:160])

    # 8) files
    files = {"graphml": prefix + ".graphml", "hump_csv": prefix + "_hump.csv",
             "hump_png": prefix + "_hump.png", "hump_svg": prefix + "_hump.svg",
             "ranking_csv": prefix + "_ranking.csv", "parameters": prefix + "_parameters.json",
             "result": prefix + "_result.json"}
    with open(os.path.join(outdir, files["hump_csv"]), "w") as f:
        f.write("gamma,transport_distal_mean" + (",control_mean,control_sd" if t_ctrl else "") + "\n")
        for i, g in enumerate(gammas):
            f.write(f"{g},{transport[i]:.6f}"
                    + (f",{t_ctrl['mean'][i]:.6f},{t_ctrl['sd'][i]:.6f}" if t_ctrl else "") + "\n")
    save_figure(os.path.join(outdir, prefix + "_hump"), gammas, transport, f"{o['site_energy']} site energies",
                "mean signal reaching distal residues",
                control=(t_ctrl["mean"], t_ctrl["sd"]) if t_ctrl else None,
                control_label=f"random site energies ({len(control_keys)} seeds, mean ± sd)",
                peak_gamma=t_stats["peak_gamma"])
    if a_res:
        files.update(allosteric_csv=prefix + "_allosteric.csv", allosteric_png=prefix + "_allosteric.png",
                     allosteric_svg=prefix + "_allosteric.svg")
        ac = a_res["control"]
        with open(os.path.join(outdir, files["allosteric_csv"]), "w") as f:
            f.write("gamma,auc" + (",control_mean,control_sd" if ac else "") + "\n")
            for i, g in enumerate(gammas):
                f.write(f"{g},{a_res['auc'][i]:.6f}" + (f",{ac['mean'][i]:.6f},{ac['sd'][i]:.6f}" if ac else "") + "\n")
        save_figure(os.path.join(outdir, prefix + "_allosteric"), gammas, a_res["auc"],
                    f"{o['site_energy']} site energies", "ROC AUC, known allosteric residues",
                    control=(ac["mean"], ac["sd"]) if ac else None,
                    control_label=f"random site energies ({len(control_keys)} seeds, mean ± sd)",
                    baselines=[(v, k) for k, v in a_res["baselines"].items()], chance=0.5,
                    peak_gamma=a_stats["peak_gamma"])
    import csv
    with open(os.path.join(outdir, files["ranking_csv"]), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "residue", "resname", f"signal_from_{map_from.replace(' ', '_')}_at_gamma_{map_gamma:g}",
                    "contacts", "known_allosteric"])
        w.writerows([[r["rank"], r["id"], r["resname"], f"{r['score']:.6f}", r["degree"],
                      "yes" if r["known_allosteric"] else ""] for r in ranking])
    if qmod["file"]:
        files["qmod"] = qmod["file"]

    params = {k: o[k] for k in DEFAULTS if k != "workers"}
    params.update(input=inp if not os.path.isfile(inp) else os.path.basename(inp), cutoff=cutoff, gammas=gammas,
                  chains_used=sorted({n.split(':')[0] for n in nodes}), qmod_source=source_id,
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
        "summary": {"residues": N, "contacts": G.number_of_edges(), "qubits": max(1, int(np.ceil(np.log2(N)))),
                    "chains": sorted({n.split(':')[0] for n in nodes}), "walk_from": walk_from,
                    "qmod_source": source_id,
                    "distal_count": int(distal.sum()), "distal_min_hops": distal_hops,
                    "site_energy": o["site_energy"], "scale": o["scale"], "max_residues": o["max_residues"]},
        "gammas": gammas,
        "transport": {"distal_mean": transport.tolist(),
                      "stats": t_stats, "control": t_ctrl},
        "allosteric": a_res,
        "qmod": qmod, "allo_entries": allo_entries,
        "map": {"nodes": [{"id": n, "resname": resnames[i], "x": float(xy[i, 0]), "y": float(xy[i, 1]),
                           "score": float(map_scores[i]), "degree": int(degree[n]),
                           "role": ("source" if i in sources_idx else "allosteric" if i in known else "")}
                          for i, n in enumerate(nodes)],
                "edges": [[idx[u], idx[v]] for u, v in G.edges()], "gamma": map_gamma, "from": map_from,
                "aspect": float(max(xy[:, 1].max(), 1e-3) / max(xy[:, 0].max(), 1e-3))},
        "ranking": ranking, "files": files, "parameters": params,
    }
    with open(os.path.join(outdir, files["result"]), "w") as f:
        json.dump(jsonable({k: v for k, v in result.items() if k != "map"}), f, indent=2)
    return result


def _layout(G, nodes):
    """Residue positions projected onto their two principal axes (spring layout
    if coordinates are missing), scaled to [0, 1] with the aspect ratio kept."""
    try:
        xyz = np.array([[float(G.nodes[n][k]) for k in ("x", "y", "z")] for n in nodes])
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
    ap = argparse.ArgumentParser(description="One protein: structure -> network, noise sweep, allosteric test, Qmod.")
    ap.add_argument("input", help="PDB/mmCIF file, or a 4-character PDB id.")
    ap.add_argument("--source", default=None, help="Start residue without labels, e.g. A:151 (default: most connected); also the Qmod start.")
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
    ap.add_argument("--control", action="store_true", help="Also run random site energies (null model) over several seeds.")
    ap.add_argument("--control-seeds", type=int, default=5)
    ap.add_argument("--labels", choices=["auto", "none"], default="auto", help="auto: ALLO table by PDB id; none: skip.")
    ap.add_argument("--active", default=None, help="Active-site residues, e.g. A:57,A:102.")
    ap.add_argument("--allosteric", default=None, help="Known allosteric residues, e.g. A:196,A:203.")
    ap.add_argument("--site", default=None, help="Which ALLO entry when a PDB has several (e.g. 2 for 1CE8_2).")
    ap.add_argument("--time", type=float, default=1.0, help="Evolution time in the circuit.")
    ap.add_argument("--order", type=int, default=2)
    ap.add_argument("--repetitions", type=int, default=4)
    ap.add_argument("--max-residues", type=int, default=200, help="Skip the Qmod above this size unless --force.")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--no-qmod", action="store_true", help="Skip the circuit.")
    ap.add_argument("--workers", type=int, default=None, help="Parallel processes (default: CPU cores, max 8).")
    ap.add_argument("--outdir", default="output")
    ap.add_argument("--prefix", default=None)
    a = ap.parse_args()

    prefix = a.prefix or (os.path.splitext(os.path.basename(a.input))[0] or "protein")
    opts = dict(source=a.source, chains=a.chains, model=a.model, method=a.method, cutoff=a.cutoff,
                min_seq_sep=a.min_seq_sep, weighted=a.weighted, site_energy=a.site_energy, scale=a.scale,
                seed=a.seed, tmax=a.tmax, ntime=a.ntime, control=a.control, control_seeds=a.control_seeds,
                labels=a.labels, active=a.active, allosteric=a.allosteric, site=a.site, qmod=not a.no_qmod,
                qmod_time=a.time, qmod_order=a.order, qmod_repetitions=a.repetitions,
                max_residues=a.max_residues, force=a.force, workers=a.workers)
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
    q = r["qmod"]
    if q["status"] == "written":
        print(f"  qmod: {q['pauli_terms']} Pauli terms, Trotter fidelity {q['trotter']['fidelity']:.6f} "
              f"(order {a.order}, {a.repetitions} steps)")
    else:
        print(f"  qmod: {q['status'].replace('_', ' ')}")
    for n in r["notes"]:
        print(f"  [note] {n}")
    print(f"\nDone in {r['elapsed_s']} s. Outputs in {a.outdir}/")


if __name__ == "__main__":
    main()
