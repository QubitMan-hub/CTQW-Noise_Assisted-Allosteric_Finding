#!/usr/bin/env python3
"""sensitivity.py - does the result survive other reasonable settings?

Reruns one protein over a grid of contact cutoffs and site-energy scales and
tabulates the key numbers, so you can show the result was not tuned:

    python sensitivity.py 1T49                         # cutoff 7,8,9 x scale 1,3,5
    python sensitivity.py 1T49 --cutoffs 7,8 --scales 3 --control

Writes output/sensitivity/NAME_sensitivity.csv (and .md). Each row is a full
run_protein analysis; the random-energy control is off unless --control is
given, which makes the grid about 4 times slower.
"""
import argparse, csv, os, sys, time

import run_protein as rp

COLUMNS = ["cutoff", "scale", "residues", "contacts", "beyond_distance_best", "best_gamma", "noise_hump",
           "gain_over_ends", "p_value", "classical_best", "quantum_minus_classical",
           "energy_weighted_best", "quantum_minus_energy_weighted", "raw_best", "proximity",
           "transport_verdict", "control_best_seed", "seconds"]


def one(inp, cutoff, scale, control, seeds, outdir, prefix, tmax=30.0):
    t0 = time.time()
    r = rp.analyze(inp, outdir, prefix, {"cutoff": cutoff, "scale": scale, "control": control,
                                          "control_seeds": seeds, "tmax": tmax}, log=lambda *a: None)
    a = r["allosteric"]
    row = {"cutoff": cutoff, "scale": scale, "residues": r["summary"]["residues"],
           "contacts": r["summary"]["contacts"], "transport_verdict": r["transport"]["stats"]["verdict"],
           "seconds": round(time.time() - t0, 1)}
    if a:
        d, st = a["adjusted"], a["adjusted"]["stats"]
        cb = d["baselines"]["classical walk, best rate"]
        ew = d["baselines"]["energy-weighted classical walk, best γ"]
        row.update(beyond_distance_best=st["peak_value"], best_gamma=st["peak_gamma"],
                   noise_hump="yes" if st["verdict"] == "hump" else "no",
                   gain_over_ends=st["gain_abs"] if st["verdict"] == "hump" else 0.0,
                   p_value=d["significance"]["p_value"] if d["significance"] else None,
                   classical_best=cb, quantum_minus_classical=st["peak_value"] - cb,
                   energy_weighted_best=ew, quantum_minus_energy_weighted=st["peak_value"] - ew,
                   raw_best=a["stats"]["peak_value"], proximity=a["baselines"]["proximity to active site"],
                   control_best_seed=max(d["control"]["best_seeds"]) if d["control"] else None)
    return row


def fmt(v):
    return f"{v:.3f}" if isinstance(v, float) else ("" if v is None else str(v))


def grid(inp, cutoffs, scales, control=False, seeds=5, outdir=os.path.join("output", "sensitivity"), log=print):
    """Run every cutoff x scale setting; write NAME_sensitivity.csv/.md; return the rows."""
    name = os.path.splitext(os.path.basename(inp))[0]
    os.makedirs(outdir, exist_ok=True)
    rows = []
    for c in cutoffs:
        for s in scales:
            prefix = f"{name}_c{c:g}_s{s:g}"
            rows.append(one(inp, c, s, control, seeds, os.path.join(outdir, prefix), prefix))
            r = rows[-1]
            log(f"cutoff {c:g} scale {s:g}: beyond distance {fmt(r.get('beyond_distance_best'))} "
                f"at gamma {fmt(r.get('best_gamma'))}, p {fmt(r.get('p_value'))}, "
                f"classical {fmt(r.get('classical_best'))} ({r['seconds']} s)")
    base = os.path.join(outdir, f"{name}_sensitivity")
    with open(base + ".csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows([{k: r.get(k) for k in COLUMNS} for r in rows])
    shown = ["cutoff", "scale", "beyond_distance_best", "best_gamma", "noise_hump", "p_value",
             "classical_best", "quantum_minus_classical"]
    lines = ["| " + " | ".join(shown) + " |", "|" + "---|" * len(shown)]
    lines += ["| " + " | ".join(fmt(r.get(k)) for k in shown) + " |" for r in rows]
    with open(base + ".md", "w") as f:
        f.write("\n".join(lines) + "\n")
    log("\n" + "\n".join(lines) + f"\n\nWrote {base}.csv and .md")
    return rows


def main():
    ap = argparse.ArgumentParser(description="Rerun one protein over several cutoffs and site-energy scales.")
    ap.add_argument("input", help="PDB/mmCIF file or 4-character PDB id (known sites looked up as usual).")
    ap.add_argument("--cutoffs", default="7,8,9", help="Contact cutoffs in Å (default 7,8,9).")
    ap.add_argument("--scales", default="1,3,5", help="Site-energy scales (default 1,3,5).")
    ap.add_argument("--control", action="store_true", help="Also run the random-energy control for every setting.")
    ap.add_argument("--control-seeds", type=int, default=5)
    ap.add_argument("--outdir", default=os.path.join("output", "sensitivity"))
    a = ap.parse_args()
    try:
        grid(a.input, [float(x) for x in a.cutoffs.split(",")], [float(x) for x in a.scales.split(",")],
             a.control, a.control_seeds, a.outdir, log=lambda m: print(m, flush=True))
    except rp.InputError as e:
        sys.exit(f"[error] {e}")


if __name__ == "__main__":
    main()
