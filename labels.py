#!/usr/bin/env python3
"""labels.py - known active-site and allosteric residues for a protein.

Labels come from, in order of precedence:
  1. residues you give (--active / --allosteric, or a JSON file via --labels),
  2. the bundled ALLO benchmark table (data/allo_labels.csv), looked up by the
     structure's PDB id. That table is Supplementary Table S2 of
     Wu, Stromich & Yaliraki, "Prediction of allosteric sites and signaling:
     insights from benchmarking datasets", Patterns 3(1), 100408 (2022),
     doi:10.1016/j.patter.2021.100408 (CC BY 4.0), curated from the Allosteric
     Database (ASD). Cite it if you use it.

Residue ids use the network's convention, CHAIN:RESSEQ[ICODE] in author
numbering, e.g. A:57 or A:57B.
"""
import csv, json, os, re

HERE = os.path.dirname(os.path.abspath(__file__))
TABLE_PATH = os.path.join(HERE, "data", "allo_labels.csv")
ALLO_CITATION = ("Wu N, Stromich L, Yaliraki SN. Prediction of allosteric sites and signaling: "
                 "insights from benchmarking datasets. Patterns 3(1), 100408 (2022), "
                 "doi:10.1016/j.patter.2021.100408, Table S2.")
_ID = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
_RES = re.compile(r"^([A-Za-z0-9]+):(-?\d+[A-Za-z]?)(?::([A-Za-z]{3}))?$")


def load_table(path=TABLE_PATH):
    if not os.path.isfile(path):
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def structure_id(inp):
    """4-character PDB id for a PDB id string or a structure file, else None.
    Files are read for their own id (PDB HEADER / mmCIF data_ block) before
    falling back to the file name (1a3w.pdb, pdb1a3w.ent, 1A3W.cif)."""
    if not os.path.isfile(inp):
        tok = inp.strip()
        return tok.upper() if _ID.match(tok) else None
    try:
        with open(inp, "r", errors="replace") as f:
            for i, line in enumerate(f):
                if line.startswith("HEADER") and len(line) >= 66 and _ID.match(line[62:66]):
                    return line[62:66].upper()
                if line.startswith("data_") and _ID.match(line[5:].strip()):
                    return line[5:].strip().upper()
                if line.startswith("_entry.id") and _ID.match(line.split()[-1]):
                    return line.split()[-1].upper()
                if i > 60:
                    break
    except OSError:
        pass
    stem = os.path.splitext(os.path.basename(inp))[0]
    m = re.fullmatch(r"(?:pdb)?([0-9][A-Za-z0-9]{3})", stem, re.I)
    return m.group(1).upper() if m else None


def parse_residues(text, default_chain=None):
    """'A:57,A:102' or '57 102' (with a default chain) -> [(id, resname or None), ...]."""
    out = []
    for tok in re.split(r"[,\s]+", (text or "").strip()):
        if not tok:
            continue
        if ":" not in tok:
            if not default_chain:
                raise ValueError(f"residue '{tok}' has no chain; write it as CHAIN:NUMBER, e.g. A:{tok}.")
            tok = f"{default_chain}:{tok}"
        m = _RES.match(tok)
        if not m:
            raise ValueError(f"cannot read residue '{tok}'; expected CHAIN:NUMBER, e.g. A:57.")
        out.append((f"{m.group(1)}:{m.group(2).upper()}", m.group(3).upper() if m.group(3) else None))
    return out


def from_table_row(row):
    return {"origin": "ALLO", "entry": row["entry"], "pdb": row["pdb"], "protein": row["protein"],
            "allosteric_ligand": row["allosteric_ligand"], "citation": ALLO_CITATION,
            "active": parse_residues(row["active_site"]),
            "allosteric": parse_residues(row["allosteric_site"])}


def from_json(path):
    d = json.load(open(path))
    act = d.get("active", d.get("source", []))
    allo = d.get("allosteric", [])
    join = lambda v: v if isinstance(v, str) else ",".join(v)
    return {"origin": "file", "entry": os.path.basename(path), "active": parse_residues(join(act)),
            "allosteric": parse_residues(join(allo))}


def resolve(inp, labels_arg="auto", active=None, allosteric=None, site=None, default_chain=None):
    """Work out which labels apply. Returns (labels or None, note, other ALLO entries)."""
    if active or allosteric:
        if not (active and allosteric):
            raise ValueError("give both --active and --allosteric, or neither.")
        return ({"origin": "user", "entry": "command line", "active": parse_residues(active, default_chain),
                 "allosteric": parse_residues(allosteric, default_chain)}, "labels from the command line", [])
    if labels_arg in (None, "", "none"):
        return None, "labels off", []
    if labels_arg != "auto":
        return from_json(labels_arg), f"labels from {os.path.basename(labels_arg)}", []
    pid = structure_id(inp)
    rows = [r for r in load_table() if r["pdb"].upper() == pid] if pid else []
    if not rows:
        return None, (f"PDB id {pid} is not in the ALLO table" if pid else "no PDB id found for this structure"), []
    entries = [r["entry"] for r in rows]
    if site is not None:
        pick = [r for r in rows if r["entry"] == site or r["entry"].endswith(f"_{site}")]
        if not pick:
            raise ValueError(f"ALLO has no site '{site}' for {pid}; choose one of {', '.join(entries)}.")
        row = pick[0]
    else:
        row = rows[0]
    note = f"ALLO entry {row['entry']} ({row['protein']})"
    if len(rows) > 1:
        note += f"; other sites for this PDB: {', '.join(e for e in entries if e != row['entry'])} (use --site)"
    return from_table_row(row), note, entries


def chains_of(labels):
    return sorted({rid.split(":")[0] for rid, _ in labels["active"] + labels["allosteric"]})


def match(labels, nodes, resnames):
    """Map label residues onto the network and check residue names where known.
    Nothing is dropped silently: every miss and every name clash is reported."""
    idx = {n: i for i, n in enumerate(nodes)}
    rep = {"active": [], "allosteric": [], "active_missing": [], "allosteric_missing": [],
           "name_mismatch": []}
    for key in ("active", "allosteric"):
        for rid, rname in labels[key]:
            if rid not in idx:
                rep[f"{key}_missing"].append(rid)
                continue
            if rname and resnames[idx[rid]] and resnames[idx[rid]].upper() != rname:
                rep["name_mismatch"].append(f"{rid} labels={rname} structure={resnames[idx[rid]]}")
            if idx[rid] not in rep[key]:
                rep[key].append(idx[rid])
    return rep
