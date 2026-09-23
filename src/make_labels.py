#!/usr/bin/env python3
"""Curated residue numbers -> a labels JSON matching a network's node ids.

Network node ids are CHAIN:RESSEQ (author/PDB numbering), e.g. "A:196", and
curated databases use the same numbering, so mapping is usually just joining
them. This checks each id exists and warns about mismatches (wrong chain letter
or numbering) before they silently corrupt scoring.

    python make_labels.py NET.graphml --source A:57,A:102 --allosteric A:196,A:203 --out labels.json
    python make_labels.py NET.graphml --chain A --source 57,102 --allosteric 196,203 --out labels.json
"""
import argparse, json, os, sys
import networkx as nx


def parse_residues(spec, default_chain):
    out = []
    for tok in (t.strip() for t in spec.split(",")):
        if not tok:
            continue
        if ":" in tok:
            out.append(tok)
        elif default_chain:
            out.append(f"{default_chain}:{tok}")
        else:
            raise ValueError(f"residue '{tok}' has no chain and no --chain given.")
    return out


def main():
    ap = argparse.ArgumentParser(description="Build a labels JSON from residue numbers.")
    ap.add_argument("network")
    ap.add_argument("--source", required=True, help="e.g. 'A:57,A:102' or '57,102'")
    ap.add_argument("--allosteric", required=True)
    ap.add_argument("--chain", default=None, help="default chain for bare numbers")
    ap.add_argument("--out", required=True)
    ap.add_argument("--strict", action="store_true", help="fail on any mismatch")
    a = ap.parse_args()

    node_set = set(nx.read_graphml(a.network).nodes())
    src_req = parse_residues(a.source, a.chain)
    allo_req = parse_residues(a.allosteric, a.chain)
    src = [r for r in src_req if r in node_set]
    allo = [r for r in allo_req if r in node_set]
    src_miss = [r for r in src_req if r not in node_set]
    allo_miss = [r for r in allo_req if r not in node_set]

    if src_miss:
        print(f"[warn] source not in network: {src_miss}")
    if allo_miss:
        print(f"[warn] allosteric not in network: {allo_miss}")
    if (src_miss or allo_miss) and a.strict:
        sys.exit("[error] --strict: some residues did not match. Check chain letters/numbering.")
    if not src or not allo:
        sys.exit("[error] no source or no allosteric residues matched the network.")

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    json.dump({"source": src, "allosteric": allo,
               "_unmatched": {"source": src_miss, "allosteric": allo_miss}},
              open(a.out, "w"), indent=2)
    print(f"wrote {a.out}  (source {len(src)}/{len(src_req)}, allosteric {len(allo)}/{len(allo_req)})")


if __name__ == "__main__":
    main()
