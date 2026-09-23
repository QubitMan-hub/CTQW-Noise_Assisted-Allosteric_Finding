# Project brief: environment-assisted CTQW for allosteric site prediction

This document gives you the context for a research project and one small task to do
right now. Read the whole thing so future instructions make sense, but only do the
task in the section titled "Your task right now." Do not do anything beyond that yet.

## What this project is

We are building toward a research paper. The idea is to predict allosteric sites in
proteins using a continuous-time quantum walk (CTQW) on the protein's residue
interaction network, with a tunable amount of noise added to the walk.

Some background so the code choices make sense:

Allosteric sites are pockets on a protein away from the main active site. Binding
there changes the protein's behaviour from a distance. They are valuable drug targets
because they are more unique to each protein, but they are hard to predict.

A protein can be turned into a graph called a residue interaction network. Each amino
acid residue is a node, and two residues are connected by an edge if they are in
contact in the 3D structure. Allosteric communication is a signal travelling through
this network, so graph methods can be used to find the residues involved.

A CTQW is a quantum version of a random walk on that graph. The known result is that a
pure quantum walk mostly reproduces plain classical centrality, so on its own it does
not help much, and its interference can even trap the walker in one region.

The core hypothesis of the paper: if you add a controlled amount of dephasing noise,
you interpolate between the fully quantum walk and the fully classical random walk (a
quantum stochastic walk). There should be an intermediate noise level where the walk
predicts true allosteric residues better than either extreme. This mirrors
environment-assisted quantum transport seen in photosynthesis, where moderate noise
improves how energy moves through a network. The paper's key figure is prediction
accuracy plotted against noise level, expected to peak in the middle.

We will validate against the Allosteric Database (ASD), which lists experimentally
known allosteric residues.

## Pipeline stages and current status

1. Protein structure to residue interaction network. DONE. This is `rin_builder.py`.
2. Build the tunable walk (quantum stochastic walk) on the graph. NOT STARTED.
3. Sweep the noise level and rank residues at each level. NOT STARTED.
4. Benchmark rankings against ASD labels and produce the accuracy-vs-noise figure.
   NOT STARTED.

## What is already built

`rin_builder.py` converts a protein structure into a residue interaction network.

- Input: a PDB or mmCIF file, or a 4-character PDB ID it will download.
- Output: the graph as GraphML, a sparse adjacency matrix as .npz, a node table CSV
  (which fixes the row and column order of the adjacency matrix), an edge table CSV, a
  text summary, and a PNG picture.
- Contact definitions: `cb` (beta-carbon within 8 A, default), `ca` (alpha-carbon),
  and `heavy` (any heavy atoms within 5 A). Every edge also stores the real distance.
- It handles the messy real cases: multi-model NMR files, multiple chains, waters and
  other hetero groups, alternate conformations, missing atoms, insertion codes, and
  modified residues such as selenomethionine (MSE mapped to MET).
- Known limits, to address later: it uses amino acids only (no ligand nodes yet), and
  it works on a single static structure.

It was tested and works. Dependencies: biopython, networkx, numpy, scipy, matplotlib.
