# Environment-assisted quantum walk for allosteric site prediction

This repo tests one idea: a quantum walk on a protein's residue network predicts
allosteric residues better when you add a moderate, tuned amount of noise
(dephasing) than when the walk is fully quantum or fully classical. The claim
mirrors environment-assisted transport in photosynthesis, where a middling
amount of noise moves energy through a network most efficiently.

## The heart of it

Everything about the physics lives in one file, `walk_core.py`: loading a
network, assigning residue site energies, building the walk Hamiltonian,
evolving the dephased walk, scoring residues, and the Pauli decomposition the
quantum route needs. Every other script is a thin command-line tool on top of
it, so they all agree by construction.

One modelling choice matters more than the rest: the residue site energies. A
plain contact graph with identical energies shows no noise benefit at all. The
default here assigns site energies by residue hydrophobicity (Kyte-Doolittle),
which is deterministic and biologically grounded, and there is a `random` mode
kept as a null-model control to compare against.

## The workflow

1. Build the network from a structure.
   `python rin_builder.py 1A8O.pdb`  (or use it inside the pipeline below)

2. Sanity-check the physics on one protein.
   `python classical_qsw.py NET.graphml --source A:151`
   Looks for the transport hump. Fast, exact, no quantum computer.

3. Make a labels file for a protein (curated active-site and allosteric residues).
   `python make_labels.py NET.graphml --source A:57,A:102 --allosteric A:196,A:203 --out labels.json`

4. Score one protein: does intermediate noise recover its known allosteric residues?
   `python score_allostery.py NET.graphml --labels labels.json`

5. The headline experiment: run across many labelled proteins and aggregate.
   `python run_benchmark.py --data-dir data/benchmark`
   Reports how often intermediate noise wins and the mean AUC-vs-noise curve.
   This is the result the paper reports.

6. Quantum route (optional, for the Classiq/hardware story on small proteins).
   `python protein_to_qmod.py 1A8O.pdb --source A:151`  ->  a .qmod file
   or `python qmod_walk.py NET.graphml --source A:151` from an existing network.

## Getting real labels

The scoring and benchmark steps need real allosteric labels. The cleanest source
is the ALLO benchmark supplementary (Wu et al. 2022) on Europe PMC, which curates
active-site and allosteric-site residues for 118 proteins from the Allosteric
Database. Turn each protein's curated residues into a labels JSON with
`make_labels.py`, put NAME.graphml and NAME.json together in one folder, and
point `run_benchmark.py` at it.

## Two routes, two size limits

The classical engine (`classical_qsw.py`, `score_allostery.py`,
`run_benchmark.py`) has no size limit and is what the paper's result rests on.
The Qmod route (`qmod_walk.py`, `protein_to_qmod.py`) is for the quantum-platform
and hardware story and only fits small proteins or single domains, because an
arbitrary Hamiltonian decomposes into up to N-squared Pauli terms. Use one chain
or one domain there; use the classical route for everything else.

## Honesty notes

- The example runs in this repo used synthetic labels to test the machinery.
  They are not biological results. Real labels give the real answer.
- The Qmod encoding's qubit bit-order should be checked once against the
  classical walk (noise off, small network, compare populations) before trusting
  quantum-route results.
- A null or negative result is still worth publishing: if intermediate noise
  does not help across real proteins, that is a real finding about these
  networks, and the honest framing is the strength of the work.

## Files

- `walk_core.py`      shared physics, scoring, and Pauli decomposition
- `rin_builder.py`    protein structure -> residue network
- `classical_qsw.py`  transport-vs-noise check for one protein
- `score_allostery.py` allosteric recovery vs noise for one protein
- `run_benchmark.py`  aggregate the result across many proteins
- `make_labels.py`    curated residue numbers -> labels JSON
- `qmod_walk.py`      network -> Classiq Qmod walk file
- `protein_to_qmod.py` structure -> Qmod, end to end
