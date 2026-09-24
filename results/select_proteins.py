#!/usr/bin/env python3
"""select_proteins.py - the prespecified validation set, chosen before any of it is run.

The four proteins in PROTEINS of make_results.py were picked during development.
The validation set is drawn from the ALLO table by a rule fixed in advance:

  1. one entry per PDB id (the first in the table);
  2. every active-site and allosteric residue lies on one chain;
  3. at least 5 allosteric residues;
  4. that chain has 100 to 350 amino-acid residues in the deposited structure
     (the size range of the development set, and a run time of minutes);
  5. not one of the four development proteins.

The eligible entries, sorted by PDB id, are shuffled with seed 2026 and the first
six are the validation set. Every one of them is reported whatever it shows; if one
cannot be analysed at all, the next in the shuffled order replaces it and that is
reported too. Writes results/selection.csv.

    python results/select_proteins.py
"""
import csv, os, random, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import labels as lb            # noqa: E402
import rin_builder as rb       # noqa: E402

DEVELOPMENT = {"1T49", "3LSW", "1IWH", "3CSM"}
MIN_ALLO, SIZE, SEED, DRAW = 5, (100, 350), 2026, 6


def chain_size(pid, chain, tmpdir):
    s = rb.load_structure(pid, tmpdir)
    return len(rb.collect_residues(s, 0, {chain}))


def main():
    seen, eligible = set(), []
    with tempfile.TemporaryDirectory() as tmp:
        for row in lb.load_table():
            pid = row["pdb"].upper()
            if pid in seen:
                continue
            seen.add(pid)
            active, allo = row["active_site"].split(), row["allosteric_site"].split()
            chains = {r.split(":")[0] for r in active + allo}
            if pid in DEVELOPMENT or len(chains) != 1 or len(allo) < MIN_ALLO or not active:
                continue
            chain = chains.pop()
            try:
                n = chain_size(pid, chain, tmp)
            except Exception as e:                       # noqa: BLE001  unreadable entries are not eligible
                print(f"{pid}: could not read ({e})", flush=True)
                continue
            print(f"{pid} chain {chain}: {n} residues", flush=True)
            if SIZE[0] <= n <= SIZE[1]:
                eligible.append({"pdb": pid, "protein": row["protein"], "chain": chain, "residues": n,
                                 "n_allosteric": len(allo)})
    eligible.sort(key=lambda r: r["pdb"])
    order = random.Random(SEED).sample(eligible, len(eligible))
    for i, r in enumerate(order):
        r["draw_order"], r["selected"] = i + 1, "yes" if i < DRAW else ""
    with open(os.path.join(HERE, "selection.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["draw_order", "pdb", "protein", "chain", "residues", "n_allosteric", "selected"])
        w.writeheader()
        w.writerows(order)
    print(f"{len(eligible)} eligible; drawn: {', '.join(r['pdb'] for r in order[:DRAW])}")


if __name__ == "__main__":
    main()
