# Results

Every number and figure here is produced by scripts in this repository from the
pipeline's own runs, so they can be rebuilt and checked:

    python results/make_results.py            # rerun all walks (about 2.5 h on 4 cores), then tables + figures
    python results/make_results.py --reuse    # rebuild tables + figures from results/runs/ only (about a minute)
    python results/external_baselines.py      # optional: established predictors (needs `pip install prody`)
    python results/sensitivity_more.py        # optional: one setting at a time for the significant proteins

The pipeline defaults are used throughout: 8 Å Cβ contacts, hydropathy site
energies at scale 3, 18 dephasing rates (0 and 0.01 to 100), T = 30, 5
random-energy seeds and 10,000 label shuffles for p. The active and allosteric
residues come from the bundled ALLO table (Wu, Strömich & Yaliraki, Patterns
2022, Table S2).

## Which proteins, and when they were chosen

`PREREGISTRATION.md` has the full history with commit hashes.

- **Development set** (1T49, 3LSW, 1IWH, 3CSM): picked by hand while the method was
  built, so results on them may be optimistic.
- **Validation set** (ten proteins): drawn by `select_proteins.py`, whose rule was
  committed (607c3b6) before any of them was run; `selection.csv` has all 20
  eligible entries in their fixed shuffled order. The first six were run first.
  The decision to extend to the first ten of the same order was taken *after*
  seeing those six (37492fb): the order was fixed in advance, that stopping point
  was not.
- **Replication set** (2YHD, 4HO6, 2VD3, 3PXF, 1ZDS): the first structure of every
  protein not yet studied, from the rest of the draw order, declared final in
  `PREREGISTRATION.md` and committed (8877a94) before any was run.

Some held-out entries are the same protein (3ZCW/4BBG kinesin Eg5, 3HFR/4B1F
glutamate racemase, and 3M3F the same domain as development protein 3LSW), so the
15 held-out entries are 12 distinct new proteins. Counts below say which they use.

## What is here

| File | What it holds |
|---|---|
| `proteins.csv` / `.md` | one row per protein: best AUC beyond distance, its γ, a 95% bootstrap interval, p, both classical walks, the best random-energy seed, the quantum-minus-classical test |
| `external.csv` / `.md` | closeness, betweenness and perturbation-response scanning under the same test |
| `sensitivity_1T49.*`, `sensitivity_1IWH.*` | contact cutoff 7/8/9 Å × energy scale 1/3/5 (no control) |
| `sensitivity_more.*` | one setting at a time (T, cutoff, scale) for every significant protein |
| `figures/*.png` / `.svg` | 300 dpi figures: `workflow`, `transport`, `beyond_distance` (four key proteins), `beyond_distance_all`, `methods`, `sensitivity`, `map_1T49` |
| `numbers.tex` | every number above as LaTeX macros, so a write-up can quote them without retyping |
| `PREREGISTRATION.md`, `select_proteins.py`, `selection.csv` | how the held-out proteins were chosen |
| `runs/<PDB>/` | each run's data: curve CSVs, the full ranking, the per-residue quantum-vs-classical comparison, every setting (`_parameters.json`, with a code version), every number (`_result.json`) and both walks at every γ (`_map.json`) |

## How to read the numbers

- **Beyond distance** is the ROC AUC for the known allosteric residues when each
  residue is only compared with residues the same number of contacts from the
  active site. Anything that depends only on distance scores exactly 0.5.
- **p** comes from shuffling the allosteric labels within each distance shell
  10,000 times and taking the best AUC over all 18 noise levels each time, so it
  accounts for having picked the best γ.
- **95% interval**: residues resampled within each shell 2,000 times, at the chosen γ.
- **Plain classical** is an ordinary random walk on the same contacts at its best
  of 17 hopping rates. **Energy-weighted** is the classical walk the noisy quantum
  walk reduces to: same contacts, energies and γ, no interference.

## Summary

| PDB | set | beyond distance (γ) | 95% interval | p | plain classical | energy-weighted | best random seed |
|---|---|---|---|---|---|---|---|
| 1T49 | dev | **0.727** (1.78) | 0.576–0.851 | **0.014** | 0.655 | 0.728 | 0.712 |
| 3LSW | dev | 0.746 (0.03) | 0.480–0.949 | 0.038 | 0.667 | 0.673 | 0.754 |
| 1IWH | dev | 0.687 (0.56) | 0.354–0.951 | 0.176 | 0.627 | 0.738 | 0.628 |
| 3CSM | dev | 0.496 (0.32) | 0.296–0.696 | 0.803 | 0.486 | 0.463 | 0.652 |
| 2RD5 | val | 0.343 (0) | 0.221–0.473 | 0.997 | 0.364 | 0.586 | 0.455 |
| 3HO6 | val | 0.585 (100) | 0.455–0.715 | 0.326 | 0.637 | 0.587 | 0.587 |
| 4B1F | val | 0.624 (100) | 0.485–0.768 | 0.227 | 0.620 | 0.624 | 0.684 |
| 4BBG | val | 0.560 (100) | 0.327–0.777 | 0.473 | 0.570 | 0.662 | 0.608 |
| 4PFK | val | 0.457 (56) | 0.247–0.654 | 0.887 | 0.457 | 0.478 | 0.565 |
| 3PYY | val | **0.791** (0.02) | 0.674–0.890 | **< 0.001** | 0.590 | 0.776 | 0.625 |
| 3ZCW | val | 0.763 (100) | 0.669–0.843 | **< 0.001** | 0.777 | 0.764 | 0.763 |
| 3HFR | val | 0.576 (17.8) | 0.367–0.766 | 0.562 | 0.591 | 0.578 | 0.676 |
| 3M3F | val | 0.590 (100) | 0.393–0.778 | 0.321 | 0.595 | 0.590 | 0.631 |
| 3H30 | val | **0.840** (3.16) | 0.747–0.910 | **0.001** | 0.855 | 0.843 | 0.877 |
| 2YHD | rep | 0.567 (0.01) | 0.346–0.772 | 0.565 | 0.427 | 0.457 | 0.454 |
| 4HO6 | rep | 0.349 (1.78) | 0.212–0.496 | 0.995 | 0.388 | 0.459 | 0.415 |
| 2VD3 | rep | 0.332 (17.8) | 0.138–0.545 | 0.991 | 0.329 | 0.329 | 0.343 |
| 3PXF | rep | **0.734** (1.78) | 0.635–0.824 | **0.009** | 0.710 | 0.749 | 0.762 |
| 1ZDS | rep | 0.677 (10) | 0.538–0.803 | 0.086 | 0.704 | 0.779 | 0.746 |

What this adds up to:

- Counting each distinct protein once, 4 of the 12 held-out proteins are
  significant (chance alone gives that about 0.2% of the time), but only 1 of the
  5 pre-registered replication proteins is (3PXF; chance 23%).
- In no protein does the quantum walk beat both classical walks (by more than
  0.02) and every random-energy seed: interference adds nothing detectable.
- Described after the fact (not a pre-registered test), only 1T49 and 3PYY show an
  effect that needs hydropathy-dependent hopping and survives every control; in
  3ZCW, 3H30 and 3PXF a plain classical walk or random energies do as well or better.
- Against established predictors under the same test (`external.csv`), closeness
  centrality is significant in 6 of the 15 held-out entries, the quantum walk in 4,
  the energy-weighted walk in 5, the plain classical walk and betweenness in 4 and
  perturbation-response scanning in 3. The methods succeed on different proteins,
  and tested together they do not differ in how often they succeed (Cochran's Q,
  p = 0.83); no pairwise difference in AUC survives a Holm correction.
- Everything here is conditional on ALLO's definition of an allosteric site.
- The 95% intervals are wide (0.16 to 0.60), because each protein has only 5 to 19
  labelled allosteric residues.

- Changing one setting at a time for the six significant proteins
  (`sensitivity_more.csv`), 29 of 36 runs stay significant: the time window T
  hardly matters (12 of 12), the energy scale a little (10 of 12), the contact
  cutoff most (7 of 12). 1T49, 3PYY and 3H30 hold in every setting; 3LSW holds in
  only 2 of 6.

The exact numbers are in the CSV files.
