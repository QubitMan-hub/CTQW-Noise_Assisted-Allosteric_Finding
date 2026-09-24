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
| `figures/beyond_distance.png` / `.svg` | distance-adjusted AUC versus noise for all four proteins (300 dpi) |
| `runs/<PDB>/` | each run's data: curve CSVs, the full ranking, every setting (`_parameters.json`, with a code version) and every number (`_result.json`) |

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

(Being generated: the first full run of make_results.py is in progress.)
