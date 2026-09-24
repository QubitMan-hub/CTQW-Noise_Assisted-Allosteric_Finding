#!/usr/bin/env python3
"""Local web app for the protein walk pipeline.

    python web/app.py            then open http://127.0.0.1:8000

A thin wrapper: every number comes from run_protein.analyze(), the same function
the command line uses, so the page and `python run_protein.py` agree exactly.
Each run writes its files under web/runs/<id>/ and is served from there.
"""
import json, os, re, shutil, sys, threading, time, traceback, uuid

import numpy as np
from flask import Flask, abort, jsonify, request, send_from_directory
from werkzeug.utils import secure_filename

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import run_protein as rp          # noqa: E402

RUNS = os.path.join(HERE, "runs")
KEEP_RUNS = 30
WEB_RESULT = "web_result.json"
# one summary line per finished walk, kept for good (full results are pruned after
# KEEP_RUNS runs, the table row stays). Lives next to the runs, never in git.
TABLE = os.path.join(RUNS, "table.jsonl")
RUN_ID = re.compile(r"[0-9a-f]{12}")                 # the folder name of one walk under web/runs/
PDB_ID = re.compile(r"[0-9][A-Za-z0-9]{3}")
PDB_UNREACHABLE = {"error": "Could not reach the Protein Data Bank. Check your internet connection."}
# one walk at a time: each walk already uses every CPU core, and two at once
# would compete for cores and memory (easy to trigger with a second tab)
RUN_LOCK = threading.Lock()
RUN_STARTED = {"t": None}

app = Flask(__name__, static_folder=os.path.join(HERE, "static"), static_url_path="/static")
app.config["MAX_CONTENT_LENGTH"] = 64 * 1024 * 1024


class UserError(Exception):
    """A problem with the input, shown to the user as a short message."""


def prune_runs():
    if not os.path.isdir(RUNS):
        return
    dirs = sorted((os.path.join(RUNS, d) for d in os.listdir(RUNS) if RUN_ID.fullmatch(d)),
                  key=os.path.getmtime, reverse=True)
    for d in dirs[KEEP_RUNS:]:
        shutil.rmtree(d, ignore_errors=True)


def parse_options(form):
    source = form.get("source", "").strip()             # chain ids are case-sensitive: kept as typed
    if source and not re.fullmatch(r"[A-Za-z0-9]+:-?\d+[A-Za-z]?", source):
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
        active, allosteric = form.get("active", "").strip(), form.get("allosteric", "").strip() or None
        if not active:
            raise UserError("Give the active-site residues, e.g. A:57,A:102 (allosteric residues are optional).")
    return {"source": source or None, "chains": chains or None, "site_energy": site_energy, "scale": scale,
            "control": flag("control"),
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
    if not RUN_LOCK.acquire(blocking=False):
        busy = int(time.time() - RUN_STARTED["t"]) if RUN_STARTED["t"] else 0
        return jsonify({"error": f"A walk is already running (started {busy} s ago). "
                                 "Wait for it to finish, then try again."}), 409
    RUN_STARTED["t"] = time.time()
    try:
        return _run()
    finally:
        RUN_STARTED["t"] = None
        RUN_LOCK.release()


def _run():
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
            if not PDB_ID.fullmatch(pdb_id):
                raise UserError("A PDB id is 4 characters starting with a digit, for example 1A8O.")
            inp = prefix = pdb_id.upper()
        try:
            result = rp.analyze(inp, run_dir, prefix, opts, log=lambda _m: None)
        except rp.InputError as e:
            msg = str(e)
            raise UserError(msg[:1].upper() + msg[1:])
        out = rp.jsonable(to_json(result, run_id))
        out["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        body = json.dumps(out)
        with open(os.path.join(run_dir, WEB_RESULT), "w") as f:     # lets the page reopen it later
            f.write(body)
        row = table_row(out, run_id, source="file" if has_file else "pdb")
        if not has_file:                                # the entry's title, if the page already looked it up
            row["title"] = (_INFO_CACHE.get(inp) or {}).get("title")
        add_table_row(row)
        prune_runs()
        return app.response_class(body, mimetype="application/json")
    except UserError as e:
        shutil.rmtree(run_dir, ignore_errors=True)
        return jsonify({"error": str(e)}), 400
    except Exception:
        traceback.print_exc()                        # full detail stays in the terminal
        shutil.rmtree(run_dir, ignore_errors=True)
        return jsonify({"error": "Something went wrong while running the pipeline. "
                                 "The terminal running web/app.py has the details."}), 500


# ---------------------------------------------------------------- PDB entry info + random picks
_INFO_CACHE = {}
_ID_CHARS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"


def _first(d, *keys):
    """First non-empty value among mmCIF keys (lists give their first element)."""
    for k in keys:
        v = d.get(k)
        if isinstance(v, list):
            v = next((x for x in v if x not in ("?", ".", "")), None)
        if v not in (None, "?", ".", ""):
            return v
    return None


def _as_list(v):
    return v if isinstance(v, list) else ([] if v is None else [v])


def _nice(text):
    """'CRYSTAL STRUCTURE OF GDP-BOUND HUMAN KRAS' -> 'Crystal structure of GDP-bound human KRAS'."""
    if not text:
        return text
    text = " ".join(str(text).split())
    if text.upper() != text:
        return text
    words = []
    for w in text.split(" "):
        keep = any(ch.isdigit() for ch in w) or (len(w) <= 4 and w.isalpha() and w not in (
            "OF", "THE", "AND", "IN", "WITH", "FROM", "TO", "BY", "FOR", "ON", "AT", "AN", "A", "AS", "ITS"))
        words.append(w if keep else w.lower())
    out = " ".join(words)
    return out[:1].upper() + out[1:]


def pdb_info(pdb_id):
    """Small summary of a PDB entry from its RCSB mmCIF header (no coordinates)."""
    pdb_id = pdb_id.upper()
    if pdb_id in _INFO_CACHE:
        return _INFO_CACHE[pdb_id]
    import io, urllib.error, urllib.request
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    try:
        with urllib.request.urlopen(f"https://files.rcsb.org/header/{pdb_id}.cif", timeout=10) as r:
            text = r.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None
        raise
    d = MMCIF2Dict(io.StringIO(text))
    ent_ids = _as_list(d.get("_entity.id"))
    ent_desc = dict(zip(ent_ids, _as_list(d.get("_entity.pdbx_description"))))
    chains, molecules = [], []
    for eid, ptype, strands, seq in zip(_as_list(d.get("_entity_poly.entity_id")), _as_list(d.get("_entity_poly.type")),
                                        _as_list(d.get("_entity_poly.pdbx_strand_id")),
                                        _as_list(d.get("_entity_poly.pdbx_seq_one_letter_code_can"))):
        if "polypeptide" not in ptype:
            continue
        n = len("".join(str(seq).split()))
        ids = [c.strip() for c in str(strands).split(",") if c.strip()]
        chains += [{"chain": c, "residues": n} for c in ids]
        name = _nice(ent_desc.get(eid)) or "protein"
        molecules.append({"name": name, "chains": ids, "residues": n})
    ligands = sorted({c for c in _as_list(d.get("_pdbx_entity_nonpoly.comp_id")) if c not in ("HOH", "DOD")})
    lig_names = dict(zip(_as_list(d.get("_pdbx_entity_nonpoly.comp_id")), _as_list(d.get("_pdbx_entity_nonpoly.name"))))
    res = _first(d, "_refine.ls_d_res_high", "_reflns.d_resolution_high", "_em_3d_reconstruction.resolution")
    date = _first(d, "_pdbx_database_status.recvd_initial_deposition_date")
    info = {
        "id": pdb_id,
        "title": _nice(_first(d, "_struct.title")),
        "classification": _nice(_first(d, "_struct_keywords.pdbx_keywords")),
        "organism": _nice(_first(d, "_entity_src_gen.pdbx_gene_src_scientific_name",
                                 "_entity_src_nat.pdbx_organism_scientific", "_pdbx_entity_src_syn.organism_scientific")),
        "method": _nice(_first(d, "_exptl.method")),
        "resolution": float(res) if res else None,
        "year": int(date[:4]) if date and date[:4].isdigit() else None,
        "chains": chains, "molecules": molecules,
        "ligands": [{"id": c, "name": _nice(lig_names.get(c))} for c in ligands][:8],
    }
    known = [r for r in rp.lb.load_table() if r["pdb"].upper() == pdb_id]
    info["known_site"] = ({"entries": [r["entry"] for r in known], "protein": known[0]["protein"],
                           "ligand": known[0]["allosteric_ligand"].split()[0],
                           "n_allosteric": len(known[0]["allosteric_site"].split()),
                           "chains": sorted({t.split(":")[0] for t in (known[0]["active_site"] + " "
                                                                         + known[0]["allosteric_site"]).split()})}
                          if known else None)
    # suggest one protein chain: crystal files often hold several copies
    if info["known_site"]:
        info["suggested_chains"] = ",".join(info["known_site"]["chains"])
        n = sum(c["residues"] for c in chains if c["chain"] in info["known_site"]["chains"])
    else:
        first = chains[0] if chains else None
        info["suggested_chains"] = first["chain"] if first and len(chains) > 1 else ""
        n = first["residues"] if first and len(chains) > 1 else sum(c["residues"] for c in chains)
    info["run_residues"] = n
    # measured: 169 residues ~26 s on 4 cores; cost grows ~ with the square of the size
    info["estimate_s"] = int(round(max(5, 26 * (n / 169) ** 2))) if n else None
    if len(_INFO_CACHE) >= 500:                        # keep the cache small
        _INFO_CACHE.pop(next(iter(_INFO_CACHE)))
    _INFO_CACHE[pdb_id] = info
    return info


@app.get("/api/info/<pdb_id>")
def api_info(pdb_id):
    if not PDB_ID.fullmatch(pdb_id):
        return jsonify({"error": "A PDB id is 4 characters starting with a digit."}), 400
    try:
        info = pdb_info(pdb_id)
    except Exception:
        return jsonify(PDB_UNREACHABLE), 502
    if info is None:
        return jsonify({"error": f"{pdb_id.upper()} is not a PDB entry."}), 404
    return jsonify(info)


@app.get("/api/random")
def api_random():
    import random
    rng = random.SystemRandom()
    known = request.args.get("kind") == "allosteric"
    table = rp.lb.load_table()
    for _ in range(15):
        pid = (rng.choice(table)["pdb"] if known else
               rng.choice("123456789") + "".join(rng.choice(_ID_CHARS) for _ in range(3)))
        try:
            info = pdb_info(pid)
        except Exception:
            return jsonify(PDB_UNREACHABLE), 502
        if info is None or not info["chains"]:
            continue                                   # not an entry, or no protein in it
        if known or 40 <= info["run_residues"] <= 350:
            return jsonify(info)
    return jsonify({"error": "No luck this time. Roll again."}), 503


# ---------------------------------------------------------------- stored walks and the table
def _run_dirs():
    if not os.path.isdir(RUNS):
        return []
    dirs = [d for d in os.listdir(RUNS) if RUN_ID.fullmatch(d)
            and os.path.isfile(os.path.join(RUNS, d, WEB_RESULT))]
    return sorted(dirs, key=lambda d: os.path.getmtime(os.path.join(RUNS, d, WEB_RESULT)), reverse=True)


def table_row(r, run_id, source):
    """The numbers the table shows for one walk (everything else stays in its result file)."""
    t, a, s, p = r["transport"]["stats"], r.get("allosteric"), r["summary"], r.get("parameters", {})
    d = (a or {}).get("adjusted")
    sig = (d or {}).get("significance")
    row = {"run_id": run_id, "name": r["name"], "source": source, "finished_utc": r.get("finished_utc"),
           "residues": s["residues"], "contacts": s["contacts"], "chains": s["chains"], "walk_from": s["walk_from"],
           "site_energy": s["site_energy"], "scale": s["scale"], "cutoff": p.get("cutoff"),
           "control": bool(r["transport"].get("control")), "elapsed_s": r.get("elapsed_s"),
           "transport_verdict": t["verdict"], "transport_gamma": t["peak_gamma"],
           "labelled": bool(a), "raw_best": a["stats"]["peak_value"] if a else None}
    if d:
        row.update(beyond=d["stats"]["peak_value"], beyond_gamma=d["stats"]["peak_gamma"],
                   beyond_verdict=d["stats"]["verdict"], beyond_gain=d["stats"]["gain_abs"],
                   p_value=sig["p_value"] if sig else None,
                   classical=d["baselines"].get("classical walk, best rate"),
                   energy_weighted=d["baselines"].get("energy-weighted classical walk, best γ"),
                   control_best=max(d["control"]["best_seeds"]) if d.get("control") else None)
    return row


def add_table_row(row):
    os.makedirs(RUNS, exist_ok=True)
    with open(TABLE, "a") as f:
        f.write(json.dumps(rp.jsonable(row)) + "\n")


def read_table():
    rows = []
    try:
        with open(TABLE) as f:
            for line in f:
                try:
                    rows.append(json.loads(line))
                except ValueError:
                    continue                           # a half-written line is skipped, not fatal
    except OSError:
        pass
    return rows


def backfill_table():
    """Walks finished before the table existed get their row once, at start-up."""
    known = {r.get("run_id") for r in read_table()}
    for d in reversed(_run_dirs()):
        if d in known:
            continue
        try:
            with open(os.path.join(RUNS, d, WEB_RESULT)) as f:
                r = json.load(f)
            src = "pdb" if re.fullmatch(r"[0-9][A-Z0-9]{3}", r.get("name", "")) else "file"
            add_table_row(table_row(r, d, src))
        except (OSError, ValueError, KeyError):
            continue


@app.get("/api/table")
def api_table():
    """Every walk this machine has run, newest first; 'stored' says if its full result can still be opened."""
    stored = set(_run_dirs())
    rows = [dict(r, stored=r.get("run_id") in stored) for r in reversed(read_table())]
    return jsonify(rows)


@app.get("/api/runs/<run_id>")
def api_run_result(run_id):
    if not RUN_ID.fullmatch(run_id):
        abort(404)
    path = os.path.join(RUNS, run_id, WEB_RESULT)
    if not os.path.isfile(path):
        return jsonify({"error": f"That walk's full result is no longer stored (only the last {KEEP_RUNS} are kept). "
                                 "Its numbers stay in the table; run it again to see the charts."}), 404
    with open(path) as f:
        return app.response_class(f.read(), mimetype="application/json")


@app.errorhandler(413)
def too_large(_e):
    return jsonify({"error": "That file is larger than 64 MB."}), 413


@app.get("/runs/<run_id>/<path:name>")
def run_file(run_id, name):
    if not RUN_ID.fullmatch(run_id) or "/" in name or name.startswith("."):
        abort(404)
    return send_from_directory(os.path.join(RUNS, run_id), name, as_attachment=True)


@app.get("/")
def index():
    return send_from_directory(app.static_folder, "index.html")


if __name__ == "__main__":
    port = int(os.environ.get("PORT", "8000"))
    os.makedirs(RUNS, exist_ok=True)
    backfill_table()
    print(f"\n  QubitMan  ->  http://127.0.0.1:{port}\n")
    app.run(host="127.0.0.1", port=port, debug=False, threaded=True)
