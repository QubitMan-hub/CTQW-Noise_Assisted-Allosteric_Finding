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
2. Tunable walk (quantum stochastic walk) on the graph: not started
3. Noise sweep and residue ranking: not started
4. Benchmark against ASD labels and accuracy-vs-noise figure: not started

## Layout

```
src/                  source code (rin_builder.py)
data/structures/      input PDB or mmCIF files
data/benchmark/       ASD allosteric labels (added later)
outputs/              generated networks (graphml, npz, csv, png)
  example_1A8O/       worked example
notebooks/            analysis notebooks (added later)
docs/                 notes and paper drafts (PROJECT_BRIEF.md)
```

Dependencies: biopython, networkx, numpy, scipy, matplotlib.

See `docs/PROJECT_BRIEF.md` for full context.
