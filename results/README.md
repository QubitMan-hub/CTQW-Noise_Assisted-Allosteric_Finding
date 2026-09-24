# Results

Every number and figure here is produced by one script from the pipeline in this
repository, so they can be rebuilt and checked:

    python results/make_results.py            # rerun all walks (about 2 h on 4 cores), then tables + figures
    python results/make_results.py --reuse    # rebuild tables + figures from results/runs/ only

The pipeline defaults are used throughout: 8 Å Cβ contacts, hydropathy site
energies at scale 3, 18 dephasing rates (0 and 0.01 to 100), 5 random-energy
seeds and 10,000 label shuffles for p. The active and allosteric residues come
from the bundled ALLO table (Wu, Strömich & Yaliraki, Patterns 2022, Table S2).

## Which proteins

- **Development set** (1T49, 3LSW, 1IWH, 3CSM): picked while the method was
  built, so results on them may be optimistic.
- **Validation set** (ten proteins): drawn by `select_proteins.py`, whose rule
  was committed before any of them was run: one entry per PDB id, all labels on
  one chain, at least 5 allosteric residues, 100-350 residues, not a development
  protein; the 20 eligible entries are shuffled with seed 2026. The first six
  were run first, then the set was extended to the first ten of the same order.
  `selection.csv` lists all 20 in draw order. Some pairs are the same protein
  (3ZCW/4BBG, 3HFR/4B1F, and 3M3F with 3LSW), so they are not independent tests.

Every protein run is reported, including the ones where nothing is found.

## What is here

| File | What it holds |
|---|---|
| `proteins.csv` / `.md` | one row per protein: best distance-adjusted AUC, the γ where it peaks, the noise gain, p, both classical walks, the best random-energy seed, the quantum-minus-classical test |
| `sensitivity_1T49.*`, `sensitivity_1IWH.*` | the same numbers over contact cutoff 7/8/9 Å × energy scale 1/3/5 (no control) |
| `figures/*.png` / `.svg` | 300 dpi figures: `workflow`, `transport`, `beyond_distance`, `sensitivity` and `map_1T49` (signal and quantum minus classical) |
| `numbers.tex` | every number above as LaTeX macros (`\res{1T49}{beyond}` ...), so a write-up can quote them without retyping |
| `select_proteins.py`, `selection.csv` | the validation-set rule and its draw |
| `runs/<PDB>/` | each run's data: curve CSVs, the full ranking, the per-residue quantum-vs-classical comparison, every setting (`_parameters.json`, with a code version), every number (`_result.json`) and both walks at every γ (`_map.json`) |

## How to read the numbers

- **Beyond distance** is the ROC AUC for the known allosteric residues when each
  residue is only compared with residues the same number of contacts from the
  active site. Closeness alone scores exactly 0.5, so anything above 0.5 is
  information the walk adds.
- **p** comes from shuffling the allosteric labels within each distance shell
  10,000 times and taking the best AUC over all 18 noise levels each time, so it
  accounts for having picked the best γ.
- **Plain classical** is an ordinary continuous-time random walk on the same
  contacts at its best of 17 hopping rates.
- **Energy-weighted** is the classical walk the noisy quantum walk reduces to:
  same contacts, site energies and γ, hopping rates 2γ/(γ² + ΔE²), no
  interference. The quantum walk has to beat it for interference to matter.
- **Noise gain** is how far the best γ sits above the better of the two ends
  (no noise, strongest noise); 0 means noise does not help.

## Summary

| PDB | set | beyond distance (γ) | p | plain classical | energy-weighted | best random seed | reading |
|---|---|---|---|---|---|---|---|
| 1T49 | dev | **0.727** (1.78) | **0.014** | 0.655 | 0.728 | 0.712 | significant, best at moderate noise, beats plain classical and every seed; the energy-weighted walk matches it |
| 3LSW | dev | 0.746 (0.03) | 0.038 | 0.667 | 0.673 | 0.754 | significant and beats both classical walks, but no noise benefit and one seed matches it |
| 1IWH | dev | 0.687 (0.56) | 0.176 | 0.627 | 0.738 | 0.628 | not significant |
| 3CSM | dev | 0.496 (0.32) | 0.803 | 0.486 | 0.463 | 0.652 | nothing beyond distance |
| 2RD5 | val | 0.343 (0) | 0.997 | 0.364 | 0.586 | 0.455 | nothing beyond distance |
| 3HO6 | val | 0.585 (100) | 0.326 | 0.637 | 0.587 | 0.587 | not significant |
| 4B1F | val | 0.624 (100) | 0.227 | 0.620 | 0.624 | 0.684 | not significant |
| 4BBG | val | 0.560 (100) | 0.473 | 0.570 | 0.662 | 0.608 | not significant |
| 4PFK | val | 0.457 (56) | 0.887 | 0.457 | 0.478 | 0.565 | nothing beyond distance |
| 3PYY | val | **0.791** (0.02) | **< 0.001** | 0.590 | 0.776 | 0.625 | significant, beats plain classical and every seed; the energy-weighted walk nearly matches it |
| 3ZCW | val | 0.763 (100) | **< 0.001** | 0.777 | 0.764 | 0.763 | significant, but at the classical end; plain classical does as well |
| 3HFR | val | 0.576 (17.8) | 0.562 | 0.591 | 0.578 | 0.676 | not significant |
| 3M3F | val | 0.590 (100) | 0.321 | 0.595 | 0.590 | 0.631 | not significant |
| 3H30 | val | **0.840** (3.16) | **0.001** | 0.855 | 0.843 | 0.877 | significant, but plain classical and random energies do as well |

In the validation set 3 of 10 proteins are significant (chance alone gives that
about 1% of the time), so the walk from the active site does carry information
beyond distance in some proteins. But the quantum walk beats the plain classical
walk clearly in one of them (3PYY) and the energy-weighted walk in none: what
finds the site is energy-dependent classical hopping, not interference.

For 1T49 the result holds across all 9 cutoff × scale settings (p ≤ 0.024); the
energy-weighted walk matches it in all but one (7 Å, scale 1, where the quantum
walk's best, 0.860, is near the coherent end).

The exact numbers are in `proteins.csv` and `sensitivity_*.csv`.
