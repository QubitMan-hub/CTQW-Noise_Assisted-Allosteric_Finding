#!/usr/bin/env python3
"""sensitivity_more.py - one setting at a time, for every protein with a significant result.

The full cutoff x scale grid (sensitivity_*.csv) covers two proteins. This script
checks the others that matter, the proteins significant at p < 0.05 in
proteins.csv, by changing one setting at a time away from the defaults:
time window T = 15 or 60 (default 30), contact cutoff 7 or 9 A (default 8) and
energy scale 1 or 5 (default 3). No random-energy control, so it takes about a
minute per setting. The distal cutoff is not varied: it only defines the
transport curve, not the AUC beyond distance.

    python results/sensitivity_more.py        # after make_results.py; writes sensitivity_more.csv/.md
"""
import csv, os, shutil, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import sensitivity             # noqa: E402

SETTINGS = [("T", 15.0), ("T", 60.0), ("cutoff", 7.0), ("cutoff", 9.0), ("scale", 1.0), ("scale", 5.0)]
COLUMNS = ["pdb", "setting", "value", "beyond_distance_best", "best_gamma", "p_value", "classical_best",
           "energy_weighted_best"]


def fmt_row(r):
    if r.get("p_value") is None:
        return "no allosteric test in this setting"
    return f"{r['beyond_distance_best']:.3f} p {r['p_value']:.3f}"


def main():
    with open(os.path.join(HERE, "proteins.csv")) as f:
        sig = [r["pdb"] for r in csv.DictReader(f) if float(r["p_value"]) < 0.05]
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        for pid in sig:
            for what, val in SETTINGS:
                kw = {"cutoff": 8.0, "scale": 3.0, "tmax": 30.0}
                kw[{"T": "tmax"}.get(what, what)] = val
                out = os.path.join(tmp, f"{pid}_{what}{val:g}")
                r = sensitivity.one(pid, kw["cutoff"], kw["scale"], False, 5, out, pid, tmax=kw["tmax"])
                shutil.rmtree(out, ignore_errors=True)
                rows.append({"pdb": pid, "setting": what, "value": val, **{k: r.get(k) for k in COLUMNS[3:]}})
                print(pid, what, val, fmt_row(r), flush=True)
    with open(os.path.join(HERE, "sensitivity_more.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)
    fmt = lambda v: f"{v:.3f}" if isinstance(v, float) else ("" if v is None else str(v))
    lines = ["| " + " | ".join(COLUMNS[:7]) + " |", "|" + "---|" * 7]
    lines += ["| " + " | ".join(fmt(r[k]) for k in COLUMNS[:7]) + " |" for r in rows]
    with open(os.path.join(HERE, "sensitivity_more.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("Wrote results/sensitivity_more.csv and .md")


if __name__ == "__main__":
    main()
