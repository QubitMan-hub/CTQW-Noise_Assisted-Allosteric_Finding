#!/usr/bin/env python3
"""Local web app for the protein walk pipeline.

    python web/app.py            then open http://127.0.0.1:8000

Wraps the existing pipeline without editing it: the network, the transport
sweep and the Qmod come from run_protein.py / walk_core.py / rin_builder.py,
called in the same order and with the same defaults as `run_protein.py`.
Each run writes its files under web/runs/<id>/ and is served from there.
"""
import csv, hashlib, json, os, re, shutil, subprocess, sys, time, traceback, uuid

import numpy as np
import networkx as nx
from flask import Flask, abort, jsonify, request, send_from_directory
from werkzeug.utils import secure_filename

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import run_protein as rp          # noqa: E402  (pipeline entry point, imported not edited)
import walk_core as wc            # noqa: E402

RUNS = os.path.join(HERE, "runs")
KEEP_RUNS = 30
PIPELINE_FILES = ("run_protein.py", "walk_core.py", "rin_builder.py")

# run_protein.py defaults, restated so every run records them explicitly
DEFAULTS = dict(model=0, method="cb", cutoff=None, min_seq_sep=0, weighted=False,
                site_energy="hydropathy", scale=3.0, seed=0,
                gammas=[0, 0.05, 0.1, 0.2, 0.4, 0.8, 1.5, 3, 6, 12, 25, 50],
                tmax=30.0, ntime=200, qmod_time=1.0, qmod_order=2,
                qmod_repetitions=4, max_residues=200)

app = Flask(__name__, static_folder=os.path.join(HERE, "static"), static_url_path="/static")
app.config["MAX_CONTENT_LENGTH"] = 64 * 1024 * 1024


class UserError(Exception):
    """A problem with the input, shown to the user as a short message."""


# ---------------------------------------------------------------- versioning
def code_version():
    """Git commit (+dirty flag) when available, and always a hash of the pipeline files."""
    h = hashlib.sha256()
    for name in PIPELINE_FILES:
        with open(os.path.join(ROOT, name), "rb") as f:
            h.update(name.encode() + b"\0" + f.read())
    out = {"pipeline_sha256": h.hexdigest()[:16], "git_commit": None, "git_dirty": None}
    try:
        commit = subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True,
                                text=True, timeout=5)
        if commit.returncode == 0:
            out["git_commit"] = commit.stdout.strip()
            dirty = subprocess.run(["git", "-C", ROOT, "status", "--porcelain", "--", *PIPELINE_FILES],
                                   capture_output=True, text=True, timeout=5)
            out["git_dirty"] = bool(dirty.stdout.strip())
    except (OSError, subprocess.SubprocessError):
        pass
    return out


# ---------------------------------------------------------------- the run
def sweep(H, source_i, gammas, tlist):
    """Visiting scores for every residue at every gamma (the walk run_protein uses)."""
    return np.array([wc.visiting_scores(wc.run_walk(H, source_i, g, tlist), tlist) for g in gammas])


def verdict(peak, n):
    return ("hump" if 0 < peak < n - 1 else "quantum" if peak == 0 else "classical")


def layout_2d(G, nodes):
    """Project residue coordinates (CB/CA) onto their two principal axes, so the
    map reads like the protein; fall back to a spring layout without coordinates."""
    try:
        xyz = np.array([[float(G.nodes[n][k]) for k in ("x", "y", "z")] for n in nodes])
        xyz -= xyz.mean(axis=0)
        _, _, vt = np.linalg.svd(xyz, full_matrices=False)
        xy = xyz @ vt[:2].T
    except (KeyError, ValueError, TypeError):
        pos = nx.spring_layout(G, seed=1)
        xy = np.array([pos[n] for n in nodes])
    xy = xy - xy.min(axis=0)
    span = float(xy.max()) or 1.0                    # one scale for both axes keeps the shape
    return xy / span


def paper_figure(run_dir, prefix, gammas, curves, peak_gamma):
    """High-resolution greyscale hump chart (300 dpi PNG + SVG) for a manuscript."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, NullLocator

    fwd = lambda g: np.log10(1 + np.asarray(g) / 0.05)          # same axis as the web chart
    inv = lambda u: 0.05 * (10 ** np.asarray(u) - 1)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 13, "axes.linewidth": 0.8,
                         "svg.fonttype": "none"})
    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    ax.set_xscale("function", functions=(fwd, inv))
    main = curves[0]
    ax.fill_between(gammas, 0, main["values"], color="#161616", alpha=0.07, lw=0)
    for i, c in enumerate(curves):
        solid = i == 0                               # main run solid, control dashed
        ax.plot(gammas, c["values"], color="#161616" if solid else "#8a8a8a", lw=2.0,
                ls="-" if solid else (0, (5, 3)), marker="o", ms=5,
                mfc="#161616" if solid else "#8a8a8a", mec="white", mew=1.2,
                label=c["label"], zorder=3 if solid else 2)
    ax.axvline(peak_gamma, color="#8a8a8a", lw=1, ls=(0, (3, 3)), zorder=1)
    ax.set_xlim(0, gammas[-1] * 1.05)
    ax.set_ylim(0, max(max(c["values"]) for c in curves) * 1.12)
    ticks = [t for t in (0, 0.1, 0.5, 1, 5, 10, 50) if t <= gammas[-1] * 1.05]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.set_xlabel("dephasing rate γ")
    ax.set_ylabel("signal reaching target")
    ax.grid(axis="y", color="#e5e5e5", lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#6b6b6b")
    ax.tick_params(colors="#3a3a3a", length=3)
    if len(curves) > 1:                              # above the plot, clear of the curves
        ax.legend(frameon=False, fontsize=11, loc="lower left", bbox_to_anchor=(0, 1.0),
                  ncol=2, handlelength=3.2, columnspacing=2.0, borderaxespad=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(run_dir, prefix + "_figure.png"), dpi=300)
    fig.savefig(os.path.join(run_dir, prefix + "_figure.svg"))
    plt.close(fig)


def run_pipeline(run_id, run_dir, inp, input_record, prefix, opts):
    t0 = time.time()
    method, cutoff = DEFAULTS["method"], DEFAULTS["cutoff"]
    cutoff = cutoff if cutoff is not None else (5.0 if method == "heavy" else 8.0)
    chains = {c.strip() for c in opts["chains"].split(",")} if opts["chains"] else None

    # 1) structure -> network  (run_protein.structure_to_graph, then write + reload as it does)
    try:
        G = rp.structure_to_graph(inp, run_dir, DEFAULTS["model"], chains, method, cutoff,
                                  DEFAULTS["min_seq_sep"], DEFAULTS["weighted"])
    except Exception as e:
        msg = str(e).strip().splitlines()[0][:160] if str(e).strip() else type(e).__name__
        if "no amino-acid residues" in msg:
            raise UserError("No amino-acid residues found"
                            + (f" in chain(s) {opts['chains']}." if chains else "."))
        if input_record["kind"] == "pdb_id":
            raise UserError(f"Could not fetch PDB id {input_record['value']}. "
                            "Check the id and your internet connection.")
        raise UserError(f"Could not read the structure file: {msg}")
    if G.number_of_nodes() < 3:
        raise UserError(f"Only {G.number_of_nodes()} residue(s) found; the walk needs a real network.")
    graphml = os.path.join(run_dir, prefix + ".graphml")
    nx.write_graphml(G, graphml)
    _, nodes, A, resnames = wc.load_network(graphml)
    idx = {n: i for i, n in enumerate(nodes)}
    N = len(nodes)
    n_qubits = max(1, int(np.ceil(np.log2(N))))

    # source and target, exactly as run_protein picks them
    if opts["source"] and opts["source"] not in idx:
        raise UserError(f"Start residue {opts['source']} is not in the network. "
                        f"Residue ids look like {', '.join(nodes[:3])}.")
    degree = dict(G.degree())
    source_id = opts["source"] or max(degree, key=degree.get)
    target_id = max(nx.single_source_shortest_path_length(G, source_id).items(),
                    key=lambda kv: kv[1])[0]

    # 2) transport vs noise; keep every residue's score, not just the target's
    gammas = [float(g) for g in DEFAULTS["gammas"]]
    tlist = np.linspace(0, DEFAULTS["tmax"], DEFAULTS["ntime"])
    H = wc.build_hamiltonian(A, resnames, opts["site_energy"], opts["scale"], DEFAULTS["seed"])
    scores = sweep(H, idx[source_id], gammas, tlist)
    transport = scores[:, idx[target_id]]
    peak = int(np.argmax(transport))
    curves = [{"key": opts["site_energy"], "label": f"{opts['site_energy']} site energies",
               "values": [float(v) for v in transport], "peak_index": peak,
               "verdict": verdict(peak, len(gammas))}]

    with open(os.path.join(run_dir, prefix + "_hump.csv"), "w") as f:     # same file run_protein writes
        f.write("gamma,transport\n")
        f.writelines(f"{g},{v:.6f}\n" for g, v in zip(gammas, transport))

    # 3) optional null-model control: same network, random site energies
    if opts["control"] and opts["site_energy"] != "random":
        Hr = wc.build_hamiltonian(A, resnames, "random", opts["scale"], DEFAULTS["seed"])
        tr = sweep(Hr, idx[source_id], gammas, tlist)[:, idx[target_id]]
        pr = int(np.argmax(tr))
        curves.append({"key": "random", "label": "random site energies (control)",
                       "values": [float(v) for v in tr], "peak_index": pr,
                       "verdict": verdict(pr, len(gammas))})
        with open(os.path.join(run_dir, prefix + "_control_hump.csv"), "w") as f:
            f.write("gamma,transport\n")
            f.writelines(f"{g},{v:.6f}\n" for g, v in zip(gammas, tr))
    # 4) Qmod (same guard as run_protein; classiq is optional)
    qmod = {"status": None, "file": None}
    if opts["no_qmod"]:
        qmod["status"] = "skipped_by_user"
    elif N > DEFAULTS["max_residues"]:
        qmod["status"] = "skipped_too_large"
    else:
        try:
            import classiq  # noqa: F401
        except ImportError:
            qmod["status"] = "unavailable"
        else:
            try:
                Hp, nq = wc.pad_to_power_of_two(H)
                terms, n = wc.pauli_decompose(Hp)
                err = float(np.max(np.abs(wc.pauli_reconstruct(terms, n) - Hp)))
                qfile, n_terms = rp.write_qmod(terms, nq, idx[source_id], DEFAULTS["qmod_time"],
                                               DEFAULTS["qmod_order"], DEFAULTS["qmod_repetitions"],
                                               os.path.join(run_dir, prefix + ".qmod"))
                with open(os.path.join(run_dir, prefix + "_execute.txt"), "w") as f:
                    f.write(rp.EXEC_SCAFFOLD)
                qmod.update(status="written", file=os.path.basename(qfile), pauli_terms=n_terms,
                            decomposition_error=err, size_kb=round(os.path.getsize(qfile) / 1024))
            except Exception as e:                       # graph still stands without the circuit
                qmod.update(status="failed", error=str(e).splitlines()[0][:160] if str(e) else type(e).__name__)

    # 5) residue signal map + ranking at the peak gamma of the main curve
    peak_scores = scores[peak]
    xy = layout_2d(G, nodes)
    map_nodes = [{"id": n, "resname": resnames[i], "x": float(xy[i, 0]), "y": float(xy[i, 1]),
                  "score": float(peak_scores[i]), "degree": int(degree[n])}
                 for i, n in enumerate(nodes)]
    map_edges = [[idx[u], idx[v]] for u, v in G.edges()]
    order = [i for i in np.argsort(-peak_scores, kind="stable") if nodes[i] != source_id]
    ranking = [{"rank": r + 1, "id": nodes[i], "resname": resnames[i],
                "score": float(peak_scores[i]), "degree": int(degree[nodes[i]])}
               for r, i in enumerate(order)]
    with open(os.path.join(run_dir, prefix + "_ranking.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "residue", "resname", f"signal_at_gamma_{gammas[peak]:g}", "contacts"])
        w.writerows([[r["rank"], r["id"], r["resname"], f"{r['score']:.6f}", r["degree"]] for r in ranking])

    # 6) reproducibility: parameters + paper figure
    params = {
        "run_id": run_id, "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "input": input_record, "source": source_id,
        "source_mode": "user" if opts["source"] else "auto: most connected",
        "target": target_id, "target_mode": "auto: farthest", "chains": opts["chains"] or "all",
        "model": DEFAULTS["model"], "method": method, "cutoff": cutoff,
        "min_seq_sep": DEFAULTS["min_seq_sep"], "weighted": DEFAULTS["weighted"],
        "site_energy": opts["site_energy"], "scale": opts["scale"], "seed": DEFAULTS["seed"],
        "gammas": gammas, "tmax": DEFAULTS["tmax"], "ntime": DEFAULTS["ntime"],
        "control_random_site_energies": len(curves) > 1,
        "qmod": {"requested": not opts["no_qmod"], "status": qmod["status"],
                 "evolution_time": DEFAULTS["qmod_time"], "order": DEFAULTS["qmod_order"],
                 "repetitions": DEFAULTS["qmod_repetitions"], "max_residues": DEFAULTS["max_residues"]},
        "code_version": code_version(),
    }
    with open(os.path.join(run_dir, prefix + "_parameters.json"), "w") as f:
        json.dump(params, f, indent=2)
    paper_figure(run_dir, prefix, gammas, curves, gammas[peak])

    files = {"graphml": prefix + ".graphml", "hump_csv": prefix + "_hump.csv",
             "ranking_csv": prefix + "_ranking.csv", "parameters": prefix + "_parameters.json",
             "figure_png": prefix + "_figure.png", "figure_svg": prefix + "_figure.svg"}
    if qmod["file"]:
        files["qmod"] = qmod["file"]
    if len(curves) > 1:
        files["control_csv"] = prefix + "_control_hump.csv"

    return {
        "run_id": run_id, "name": prefix, "elapsed_s": round(time.time() - t0, 1),
        "summary": {"residues": N, "contacts": G.number_of_edges(), "qubits": n_qubits,
                    "chains": sorted({n.split(":")[0] for n in nodes}), "source": source_id,
                    "source_auto": not opts["source"], "target": target_id,
                    "verdict": curves[0]["verdict"], "peak_gamma": gammas[peak],
                    "peak_signal": float(transport[peak]), "site_energy": opts["site_energy"],
                    "scale": opts["scale"], "max_residues": DEFAULTS["max_residues"]},
        "gammas": gammas, "curves": curves, "qmod": qmod,
        "map": {"nodes": map_nodes, "edges": map_edges, "source": idx[source_id],
                "aspect": float(max(xy[:, 1].max(), 1e-3) / max(xy[:, 0].max(), 1e-3)),
                "target": idx[target_id], "gamma": gammas[peak]},
        "ranking": ranking[:15], "ranking_total": len(ranking),
        "files": {k: f"/runs/{run_id}/{v}" for k, v in files.items()},
        "parameters": params,
    }


# ---------------------------------------------------------------- HTTP
def prune_runs():
    if not os.path.isdir(RUNS):
        return
    dirs = sorted((os.path.join(RUNS, d) for d in os.listdir(RUNS)
                   if re.fullmatch(r"[0-9a-f]{12}", d)), key=os.path.getmtime, reverse=True)
    for d in dirs[KEEP_RUNS:]:
        shutil.rmtree(d, ignore_errors=True)


def parse_options(form):
    source = form.get("source", "").strip().upper()
    if source and not re.fullmatch(r"[A-Z0-9]+:-?\d+[A-Z]?", source):
        raise UserError("Start residue should look like A:151 (chain, colon, number).")
    chains = form.get("chains", "").strip()
    if chains and not re.fullmatch(r"[A-Za-z0-9]+(\s*,\s*[A-Za-z0-9]+)*", chains):
        raise UserError("Chains should be letters separated by commas, for example A or A,B.")
    site_energy = form.get("site_energy", "hydropathy")
    if site_energy not in ("hydropathy", "random"):
        raise UserError("Site energy must be hydropathy or random.")
    try:
        scale = float(form.get("scale", "") or DEFAULTS["scale"])
    except ValueError:
        raise UserError("Disorder scale must be a number.")
    if not (0 <= scale <= 50) or not np.isfinite(scale):
        raise UserError("Disorder scale must be between 0 and 50.")
    flag = lambda k: form.get(k, "") in ("1", "true", "on")
    return {"source": source, "chains": chains, "site_energy": site_energy, "scale": scale,
            "no_qmod": flag("no_qmod"), "control": flag("control")}


@app.post("/api/run")
def api_run():
    run_id = uuid.uuid4().hex[:12]
    run_dir = os.path.join(RUNS, run_id)
    try:
        opts = parse_options(request.form)
        upload = request.files.get("structure")
        has_file = bool(upload and upload.filename)
        pdb_id = request.form.get("pdb_id", "").strip()
        if has_file == bool(pdb_id):
            raise UserError("Give either a structure file or a PDB id, not both."
                            if has_file else "Choose a structure file or type a PDB id.")
        os.makedirs(run_dir)
        if has_file:
            name = secure_filename(upload.filename) or "structure.pdb"
            stem, ext = os.path.splitext(name)
            if ext.lower() not in (".pdb", ".ent", ".cif", ".mmcif"):
                raise UserError("The file should be a .pdb or .cif structure.")
            inp = os.path.join(run_dir, name)
            upload.save(inp)
            with open(inp, "rb") as f:
                data = f.read()
            if not data.strip():
                raise UserError("The uploaded file is empty.")
            prefix = re.sub(r"[^A-Za-z0-9_.-]", "_", stem)[:40] or "protein"
            record = {"kind": "file", "value": upload.filename,
                      "sha256": hashlib.sha256(data).hexdigest()}
        else:
            if not re.fullmatch(r"[0-9][A-Za-z0-9]{3}", pdb_id):
                raise UserError("A PDB id is 4 characters starting with a digit, for example 1A8O.")
            inp = prefix = pdb_id.upper()
            record = {"kind": "pdb_id", "value": prefix}
        result = run_pipeline(run_id, run_dir, inp, record, prefix, opts)
        prune_runs()
        return jsonify(result)
    except UserError as e:
        shutil.rmtree(run_dir, ignore_errors=True)
        return jsonify({"error": str(e)}), 400
    except Exception:
        traceback.print_exc()                        # full detail stays in the terminal
        shutil.rmtree(run_dir, ignore_errors=True)
        return jsonify({"error": "Something went wrong while running the pipeline. "
                                 "The terminal running web/app.py has the details."}), 500


@app.errorhandler(413)
def too_large(_e):
    return jsonify({"error": "That file is larger than 64 MB."}), 413


@app.get("/runs/<run_id>/<path:name>")
def run_file(run_id, name):
    if not re.fullmatch(r"[0-9a-f]{12}", run_id) or "/" in name or name.startswith("."):
        abort(404)
    return send_from_directory(os.path.join(RUNS, run_id), name, as_attachment=True)


@app.get("/")
def index():
    return send_from_directory(app.static_folder, "index.html")


if __name__ == "__main__":
    port = int(os.environ.get("PORT", "8000"))
    os.makedirs(RUNS, exist_ok=True)
    print(f"\n  Protein walk  ->  http://127.0.0.1:{port}\n")
    app.run(host="127.0.0.1", port=port, debug=False, threaded=True)
