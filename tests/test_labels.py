#!/usr/bin/env python3
"""Checks for the label lookup. Run: python tests/test_labels.py (or pytest)."""
import os, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import labels as lb


def test_table_loads():
    rows = lb.load_table()
    assert len(rows) == 118
    assert all(r["active_site"] and r["allosteric_site"] for r in rows)


def test_lookup_by_pdb_id_and_site():
    lab, note, entries = lb.resolve("1iwh")
    assert lab["entry"] == "1IWH" and lab["origin"] == "ALLO" and "Hemoglobin" in note
    lab2, _, entries2 = lb.resolve("1CE8", site="2")
    assert lab2["entry"] == "1CE8_2" and set(entries2) == {"1CE8_1", "1CE8_2"}
    assert lb.resolve("1A8O")[0] is None                   # not in the table
    assert lb.resolve("1IWH", labels_arg="none")[0] is None


def test_structure_id_from_file_header_and_name():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "whatever.pdb")
        with open(p, "w") as f:
            f.write("HEADER    OXYGEN TRANSPORT                        05-SEP-02   1IWH              \n")
        assert lb.structure_id(p) == "1IWH"
        q = os.path.join(d, "pdb2qpy.ent")
        open(q, "w").close()
        assert lb.structure_id(q) == "2QPY"
        c = os.path.join(d, "x.cif")
        with open(c, "w") as f:
            f.write("data_3MW9\n#\n")
        assert lb.structure_id(c) == "3MW9"


def test_parse_and_match():
    assert lb.parse_residues("A:57, A:102") == [("A:57", None), ("A:102", None)]
    assert lb.parse_residues("57 58", default_chain="B") == [("B:57", None), ("B:58", None)]
    assert lb.parse_residues("A:14:asp") == [("A:14", "ASP")]
    labels = {"active": [("A:1", None)], "allosteric": [("A:2", "GLY"), ("A:9", None)]}
    rep = lb.match(labels, ["A:1", "A:2", "A:3"], ["ALA", "SER", "GLY"])
    assert rep["active"] == [0] and rep["allosteric"] == [1] and rep["allosteric_missing"] == ["A:9"]
    assert rep["name_mismatch"] == ["A:2 labels=GLY structure=SER"]


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print(f"ok  {name}")
    print("all checks passed")
