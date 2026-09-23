# Environment-assisted CTQW for allosteric site prediction

Predicting allosteric sites in proteins with a continuous-time quantum walk (CTQW) on
the protein's residue interaction network, with a tunable amount of dephasing noise.

**Hypothesis:** a pure quantum walk mostly reproduces classical centrality, and a pure
classical random walk misses interference effects. Adding controlled dephasing
(a quantum stochastic walk) interpolates between the two, and an intermediate noise
level should predict true allosteric residues best. This mirrors environment-assisted
quantum transport in photosynthesis. Predictions are validated against the Allosteric
Database (ASD).

## Pipeline status

1. Protein structure to residue interaction network: **done** (`src/rin_builder.py`)
2. Tunable walk (quantum stochastic walk) on the graph: code added (`src/walk_core.py`, `src/classical_qsw.py`)
3. Noise sweep and residue ranking: code added (`src/score_allostery.py`)
4. Benchmark against ASD labels and accuracy-vs-noise figure: code added
   (`src/run_benchmark.py`, `src/make_labels.py`); needs real labels in `data/benchmark/`

Optional quantum route (Classiq Qmod, small proteins only): `src/qmod_walk.py`,
`src/protein_to_qmod.py`. See `src/README.md` for the workflow and commands.

The scripts in `src/` import each other, so keep them in one folder and run them
from `src/`. Install with `pip install -r requirements.txt`; `classiq` is only
needed for the two Qmod scripts.

## Layout

```
src/                  pipeline scripts (all eight .py files) + README.md
requirements.txt      dependencies
data/structures/      input PDB or mmCIF files
data/benchmark/       ASD allosteric labels (added later)
outputs/              generated networks (graphml, npz, csv, png)
  example_1A8O/       worked example
notebooks/            analysis notebooks (added later)
docs/                 notes and paper drafts (PROJECT_BRIEF.md)
```

See `docs/PROJECT_BRIEF.md` for full context.
