#!/usr/bin/env python3
"""Local web app for the protein walk pipeline.

    python web/app.py            then open http://127.0.0.1:8000

A thin wrapper: every number comes from run_protein.analyze(), the same function
the command line uses, so the page and `python run_protein.py` agree exactly.
Each run writes its files under web/runs/<id>/ and is served from there.
"""
import json, os, re, shutil, sys, traceback, uuid

import numpy as np
from flask import Flask, abort, jsonify, request, send_from_directory
from werkzeug.utils import secure_filename

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import run_protein as rp          # noqa: E402

RUNS = os.path.join(HERE, "runs")
KEEP_RUNS = 30

app = Flask(__name__, static_folder=os.path.join(HERE, "static"), static_url_path="/static")
app.config["MAX_CONTENT_LENGTH"] = 64 * 1024 * 1024


class UserError(Exception):
    """A problem with the input, shown to the user as a short message."""


def prune_runs():
    if not os.path.isdir(RUNS):
        return
    dirs = sorted((os.path.join(RUNS, d) for d in os.listdir(RUNS) if re.fullmatch(r"[0-9a-f]{12}", d)),
                  key=os.path.getmtime, reverse=True)
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
        scale = float(form.get("scale", "") or 3.0)
    except ValueError:
        raise UserError("Disorder scale must be a number.")
    if not (0 <= scale <= 50) or not np.isfinite(scale):
        raise UserError("Disorder scale must be between 0 and 50.")
    flag = lambda k: form.get(k, "") in ("1", "true", "on")
    mode = form.get("labels", "auto")
    if mode not in ("auto", "none", "custom"):
        raise UserError("Labels must be auto, none or custom.")
    active = allosteric = None
    if mode == "custom":
        active, allosteric = form.get("active", "").strip(), form.get("allosteric", "").strip()
        if not active or not allosteric:
            raise UserError("Custom labels need both active-site and allosteric residues, e.g. A:57,A:102.")
    return {"source": source or None, "chains": chains or None, "site_energy": site_energy, "scale": scale,
            "qmod": not flag("no_qmod"), "control": flag("control"),
            "labels": "none" if mode == "none" else "auto", "active": active, "allosteric": allosteric}


def to_json(result, run_id):
    """The analysis result, with file names turned into download URLs."""
    out = dict(result)
    out["run_id"] = run_id
    out["files"] = {k: f"/runs/{run_id}/{v}" for k, v in result["files"].items()}
    out["ranking_total"] = len(result["ranking"])
    out["ranking"] = result["ranking"][:15]
    return out


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
            if os.path.getsize(inp) == 0:
                raise UserError("The uploaded file is empty.")
            prefix = re.sub(r"[^A-Za-z0-9_.-]", "_", stem)[:40] or "protein"
        else:
            if not re.fullmatch(r"[0-9][A-Za-z0-9]{3}", pdb_id):
                raise UserError("A PDB id is 4 characters starting with a digit, for example 1A8O.")
            inp = prefix = pdb_id.upper()
        try:
            result = rp.analyze(inp, run_dir, prefix, opts, log=lambda _m: None)
        except rp.InputError as e:
            msg = str(e)
            raise UserError(msg[:1].upper() + msg[1:])
        prune_runs()
        return app.response_class(json.dumps(rp.jsonable(to_json(result, run_id))),
                                  mimetype="application/json")
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
