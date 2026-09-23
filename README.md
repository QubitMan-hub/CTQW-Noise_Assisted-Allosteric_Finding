# Environment-assisted quantum walk, one protein at a time

Drop in a protein structure, get the quantum-walk circuit (Qmod) and a graph
showing whether signal transport through the protein peaks at an intermediate
amount of noise (the environment-assisted hump). No databases, no labels.

## Use it

    python run_protein.py 1A8O.pdb          # a file
    python run_protein.py 6LU7              # or a 4-character PDB ID

Or in the browser: `python web/app.py`, then open http://127.0.0.1:8000.

Outputs land in `output/`:

    NAME.graphml     the residue network
    NAME.qmod        the Classiq quantum-walk circuit
    NAME_hump.png    the graph: signal transport vs noise (look for the peak)
    NAME_hump.csv    the numbers behind it
    NAME_execute.txt how to run the .qmod on Classiq

## What the hump graph means

A quantum walk starts at a source residue and spreads through the contact
network. The graph sweeps the dephasing rate gamma from fully quantum (0) to
fully classical (large) and plots how much signal reaches a distant residue. A
peak in the middle means moderate noise moves the signal better than either
extreme, the same effect seen in photosynthesis.

The site energies matter: the default assigns them by residue hydrophobicity
(Kyte-Doolittle), which is what makes noise able to help. Without energy
differences between residues there is nothing for noise to fix and no hump.

## Options you might use

    --source A:151        start residue (default: the most-connected one)
    --target A:200        transport target (default: the farthest residue)
    --chains A            use only certain chains
    --method heavy        contact by heavy-atom distance (default cb: CB-CB 8A)
    --site-energy random  null-model control instead of hydropathy
    --scale 3             site-energy disorder strength
    --no-qmod             only make the hump graph
    --max-residues 200    skip the Qmod above this size (the graph still runs)

## Notes

- The Qmod is amplitude-encoded, so residue count barely affects qubit count
  (about 7 qubits for 70 residues). But an arbitrary Hamiltonian gives up to
  N-squared Pauli terms, so the circuit deepens fast; large proteins skip the
  Qmod by default (use one chain or a domain). The hump graph has no size limit.
- The Qmod encodes the noise-free walk. The noise (the hump) is added at run
  time on Classiq's density-matrix simulator with a dephasing channel; see
  NAME_execute.txt. The classical hump graph is the ground truth it should match.
- The source-encoding qubit order was checked against the classical walk (they
  match to ~2e-7); still worth one real Classiq run to confirm on hardware.

## Files

    run_protein.py   the one command (structure -> Qmod + hump graph)
    walk_core.py     the physics (walk, site energies, Pauli decomposition)
    rin_builder.py   structure -> residue network
