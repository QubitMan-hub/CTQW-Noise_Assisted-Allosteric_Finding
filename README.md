# Environment-assisted quantum walk, one protein at a time

Drop in a protein structure. The pipeline builds the residue contact network,
runs a dephased quantum walk across 25 noise levels, and answers two questions:

1. **Transport:** does an intermediate amount of noise move signal through the
   protein better than a fully quantum or a strongly dephased walk (the
   environment-assisted hump)?
2. **Allostery:** starting the walk at the active site, does intermediate noise
   rank the experimentally known allosteric residues higher (ROC AUC)? Labels
   come from the bundled ALLO benchmark table (118 proteins) or from you.

It also writes the Classiq Qmod circuit for the noise-free walk.

## Use it

    pip install -r requirements.txt
    python run_protein.py 1A8O.pdb                 # a file
    python run_protein.py 1IWH --control           # a PDB id, plus the null model

Or in the browser: `python web/app.py`, then open http://127.0.0.1:8000.

Outputs land in `output/`:

    NAME_hump.png/.svg/.csv        transport vs noise (300 dpi figure + data)
    NAME_allosteric.png/.svg/.csv  allosteric AUC vs noise, with baselines (if labels)
    NAME_ranking.csv               every residue ranked by the signal it receives
    NAME.qmod                      the Classiq circuit (12-decimal coefficients)
    NAME.graphml                   the residue network
    NAME_parameters.json           every setting, the labels used, a code version
    NAME_result.json               every number behind the figures
    NAME_sensitivity.csv           with --sensitivity

## What is measured

- **Transport** is the time-integrated signal reaching the *distal* residues,
  those at least half the network's radius (in contacts) from the source. That
  is averaged over dozens of residues, so it is much less arbitrary than a single
  "farthest" residue (which is still in the CSV).
- **Hump effect size:** peak minus the better of the two ends, in absolute and
  relative terms, and the range of noise over which the walk beats both ends.
- **Allosteric test:** the walk starts from all active-site residues together
  (one mixed initial state). Each residue is scored by the signal it receives,
  and the ROC AUC measures how well that ranks the known allosteric residues
  among all non-active-site residues. Two classical baselines are drawn for
  comparison: contact degree and graph proximity to the active site.
- **Null model (`--control`):** the same analysis with random site energies over
  5 seeds, drawn as a mean ± sd band, plus how the real hump's gain ranks among
  the seeds. A hump that random energies reproduce is not specific to the
  hydropathy model.
- **Sensitivity (`--sensitivity`):** re-runs with half and double the walk time
  and disorder 1.5 and 6, so you can check that a finding does not flip.
- **Qmod quality:** the Trotter fidelity of the written circuit against the exact
  evolution, reported in the output.

## Reading the results honestly

- Noise-assisted transport peaks are a generic property of disordered networks
  (Plenio & Huelga 2008; Rebentrost et al. 2009). The "most dephased" end is the
  quantum Zeno regime, where transport is suppressed, not a classical random
  walk. So a transport hump alone is expected; the control band shows whether it
  is anything more.
- The allosteric AUC curve is the test of the paper's claim. Compare its best
  value against the proximity baseline and the control band, not only against
  its own ends.

## Options you might use

    --source A:151          transport start residue (default: most connected)
    --chains A              chains to include (default: the labels' chains, else all)
    --labels none           skip the allosteric test (default: ALLO lookup by PDB id)
    --active A:57,A:102 --allosteric A:196,A:203   your own labels
    --site 2                which ALLO entry when a PDB has several (e.g. 1CE8_2)
    --control               null model over random site energies (--control-seeds 5)
    --sensitivity           robustness re-runs
    --site-energy random    random energies for the main run
    --scale 3               site-energy disorder strength
    --gammas 0,0.1,1,10     your own noise grid (must start at 0)
    --workers 4             parallel processes (default: CPU cores, max 8)
    --no-qmod               skip the circuit
    --max-residues 200      skip the Qmod above this size (the graphs still run)

## Speed and exactness

- With no noise the walk is solved exactly by diagonalising H.
- With noise it uses a sparse H, one mixed-state run for many source residues
  (exact, because the equation is linear), and the noise levels in parallel
  across CPU cores.
- `tests/test_walk_core.py` checks this against the original dense solver
  (agreement to about 1e-15 with noise; the old solver's own error without it).

A 250-residue protein takes a little over a minute on 4 cores for 25 noise
levels, versus several minutes for 12 levels before.

## Notes

- The Qmod is amplitude-encoded: about 7 qubits for 70 residues. An arbitrary
  Hamiltonian gives up to N-squared Pauli terms, so large proteins skip the Qmod
  by default (use one chain or a domain). The graphs have no size limit.
- The Qmod encodes the noise-free walk. Dephasing is added at run time on
  Classiq's density-matrix simulator; see NAME_execute.txt.
- The encoding (qubit order and coefficients) was checked offline against the
  classical walk; populations agree to about 1e-13. One real Classiq run is still
  worth doing.

## Labels and citation

`data/allo_labels.csv` is Supplementary Table S2 of Wu N, Strömich L, Yaliraki SN,
"Prediction of allosteric sites and signaling: insights from benchmarking
datasets", *Patterns* 3(1), 100408 (2022), doi:10.1016/j.patter.2021.100408
(CC BY 4.0), curated from the Allosteric Database. Cite it if you use it.

## Files

    run_protein.py   the one command; analyze() is shared with the web app
    walk_core.py     the physics: walk, metrics, Pauli decomposition, Trotter check
    labels.py        active-site / allosteric labels (ALLO table or your own)
    rin_builder.py   structure -> residue network
    web/             the local web app
    tests/           python tests/test_walk_core.py, python tests/test_labels.py
