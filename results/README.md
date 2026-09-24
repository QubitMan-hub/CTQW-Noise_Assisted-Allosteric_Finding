# Results

Every number and figure here is produced by one script from the pipeline in this
repository, so they can be rebuilt and checked:

    python results/make_results.py            # rerun all walks (~25 min on 4 cores), then tables + figure
    python results/make_results.py --reuse    # rebuild tables + figure from results/runs/ only

The pipeline defaults are used throughout: 8 Å Cβ contacts, hydropathy site
energies at scale 3, 18 dephasing rates (0 and 0.01 to 100), 5 random-energy
seeds and 10,000 label shuffles for p. The active and allosteric residues come
from the bundled ALLO table (Wu, Strömich & Yaliraki, Patterns 2022, Table S2).
Every protein studied so far is listed, including the ones where nothing is found.

## What is here

| File | What it holds |
|---|---|
| `proteins.csv` / `.md` | one row per protein: best distance-adjusted AUC, the γ where it peaks, the noise gain, p, the classical walk, the best random-energy seed |
| `sensitivity_1T49.*`, `sensitivity_1IWH.*` | the same numbers over contact cutoff 7/8/9 Å × energy scale 1/3/5 (no control) |
| `figures/*.png` / `.svg` | 300 dpi figures: `workflow`, `transport`, `beyond_distance`, `sensitivity` and `map_1T49` (signal and quantum minus classical) |
| `numbers.tex` | every number above as LaTeX macros (`\res{1T49}{beyond}` ...), so a write-up can quote them without retyping |
| `runs/<PDB>/` | each run's data: curve CSVs, the full ranking, the per-residue quantum-vs-classical comparison, every setting (`_parameters.json`, with a code version), every number (`_result.json`) and both walks at every γ (`_map.json`) |

## How to read the numbers

- **Beyond distance** is the ROC AUC for the known allosteric residues when each
  residue is only compared with residues the same number of contacts from the
  active site. Closeness alone scores exactly 0.5, so anything above 0.5 is
  information the walk adds.
- **p** comes from shuffling the allosteric labels within each distance shell
  10,000 times and taking the best AUC over all 18 noise levels each time, so it
  accounts for having picked the best γ.
- **Classical** is an ordinary continuous-time random walk on the same contacts
  at its best of 17 hopping rates. The quantum walk has to beat it for the
  quantum part to matter.
- **Noise gain** is how far the best γ sits above the better of the two ends
  (no noise, strongest noise); 0 means noise does not help.

## Summary

| PDB | Protein | beyond distance (γ) | p | classical | best random seed | reading |
|---|---|---|---|---|---|---|
| 1T49 | PTP1B | **0.727** (1.78) | **0.014** | 0.655 | 0.712 | significant, beats classical and every seed, best at moderate noise; holds in all 9 settings |
| 3LSW | GluA3 LBD | 0.746 (0.03) | 0.038 | 0.667 | 0.754 | significant, but no noise benefit and one seed matches it |
| 1IWH | haemoglobin α | 0.687 (0.56) | 0.176 | 0.627 | 0.628 | not significant (6 allosteric residues); depends on the cutoff |
| 3CSM | chorismate mutase | 0.496 (0.32) | 0.803 | 0.486 | 0.652 | no signal beyond distance |

Quantum minus classical (the walk compared with diffusion matched in spread, at the γ where the AUC
beyond distance peaks): its own AUC beyond distance is 0.686 (p = 0.013) for 1T49 and 0.736 (p = 0.022)
for 1IWH; 3LSW and 3CSM are not significant (p = 0.11).

The exact numbers are in `proteins.csv` and `sensitivity_*.csv`.
