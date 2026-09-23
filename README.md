# Environment-assisted quantum walk, one protein at a time

Give it any protein: a PDB id or a .pdb/.cif file. The pipeline builds the residue
contact network, runs a dephased quantum walk across 18 noise levels, and answers
two questions:

1. **Transport:** does an intermediate amount of noise move signal through the
   protein better than a fully quantum or a strongly dephased walk (the
   environment-assisted hump)?
2. **Allostery:** starting the walk at the active site, which residues receive
   the most signal (candidate allosteric residues), and, where the real
   allosteric residues are known, how well does it find them (ROC AUC)?

## What you need to give it

| You give | You get |
|---|---|
| just a PDB id or file | transport curve + residue ranking, walk from the most-connected residue |
| + the active site (`--active`) | walk from the active site: ranked candidate allosteric residues |
| + known allosteric residues (`--allosteric`) | the ROC AUC test of how well the walk finds them |

For the 118 proteins in the bundled ALLO benchmark table, the active site and
allosteric residues are filled in automatically from the PDB id. A PDB file does
not say where a protein's allosteric site is; that comes from experiments, which
is why a score is only possible where those residues are known.

## Use it

    pip install -r requirements.txt
    python run_protein.py 1A8O.pdb                 # a file
    python run_protein.py 4OBE --chains A          # any PDB id
    python run_protein.py 4OBE --chains A --active A:12,A:13,A:61   # with an active site
    python run_protein.py 1IWH                     # an ALLO protein (known sites looked up)
    python run_protein.py 1IWH --no-control        # quick look, skips the null model

Or in the browser: `python web/app.py`, then open http://127.0.0.1:8000.
The site (QubitMan, with a light/dark toggle) can also roll a random protein, or one with a known allosteric
site, and shows a short summary of any PDB entry before you run it. It needs an
internet connection to reach files.rcsb.org.

Outputs land in `output/`:

    NAME_hump.png/.svg/.csv        transport vs noise (300 dpi figure + data)
    NAME_allosteric.png/.svg/.csv  allosteric AUC vs noise, with baselines (if labels)
    NAME_adjusted.png/.svg/.csv    the same AUC among residues equally far from the active site
    NAME_ranking.csv               every residue ranked by the signal it receives
    NAME.graphml                   the residue network
    NAME_parameters.json           every setting, the labels used, a code version
    NAME_result.json               every number behind the figures

## What is measured

- **One walk per noise level.** It starts at the active site when one is known
  and feeds every result; otherwise it starts at `--source` (default: the
  most-connected residue).
- **Transport** is the time-integrated signal reaching the *distal* residues,
  those at least half the network's radius (in contacts) from the start. That is
  averaged over dozens of residues, so it is much less arbitrary than a single
  "farthest" residue.
- **Hump effect size:** peak minus the better of the two ends, in absolute and
  relative terms, and the range of noise over which the walk beats both ends.
- **Allosteric test:** the walk starts from all active-site residues together
  (one mixed initial state). Each residue is scored by the signal it receives,
  and the ROC AUC measures how well that ranks the known allosteric residues
  among all non-active-site residues. Two classical baselines are drawn for
  comparison: contact degree and graph proximity to the active site.
- **Beyond distance:** the raw signal mostly says how close a residue is to the
  active site, which proximity alone already captures. So each residue is also
  scored as a percentile among residues at the same hop distance from the active
  site (distances merged until every shell holds at least 5 residues), and the
  AUC is recomputed. Plain proximity scores exactly 0.5 here; anything above it
  is what the walk adds. Its baseline is contact degree, adjusted the same way.
- **Classical baseline:** an ordinary continuous-time random walk on the same
  contacts from the same start (dp/dt = -kLp, L the graph Laplacian), solved
  exactly for 17 hopping rates k from 0.01 to 100; its best rate is drawn as a
  baseline on both allosteric charts. Beating it is the evidence that the
  quantum part matters. It costs about a second.
- **Significance:** the allosteric labels are shuffled within each distance shell
  10,000 times (`--permutations`), and each shuffle takes its best score over all
  noise levels, so p accounts for picking the best gamma. The same test is run
  for the classical walk. No extra walks are needed.
- **Null model (on by default, `--no-control` skips it):** the same analysis with random site energies over
  5 seeds, drawn as a mean ± sd band, plus how the real hump's gain ranks among
  the seeds. A hump that random energies reproduce is not specific to the
  hydropathy model.

## Reading the results honestly

- Noise-assisted transport peaks are a generic property of disordered networks
  (Plenio & Huelga 2008; Rebentrost et al. 2009). The "most dephased" end is the
  quantum Zeno regime, where transport is suppressed, not a classical random
  walk. So a transport hump alone is expected; the control band shows whether it
  is anything more.
- The allosteric AUC curve is the test of the paper's claim. Compare its best
  value against the proximity baseline and the control band, not only against
  its own ends. The distance-adjusted curve is the stricter test: above 0.5 and
  above the control band means the walk sees something proximity does not.

## Options you might use

    --source A:151          start residue when no active site is known (default: most connected)
    --chains A              chains to include (default: the labels' chains, else all)
    --labels none           skip the allosteric test (default: ALLO lookup by PDB id)
    --active A:57,A:102     your active site (the walk starts here)
    --allosteric A:196      known allosteric residues (adds the AUC test)
    --site 2                which ALLO entry when a PDB has several (e.g. 1CE8_2)
    --control-seeds 20      more random-energy seeds for final figures (default 5)
    --no-control            skip the null model (about 4 times faster; --control-seeds 5)
    --site-energy random    random energies for the main run
    --scale 3               site-energy disorder strength
    --gammas 0,0.1,1,10     your own noise grid (must start at 0)
    --workers 4             parallel processes (default: CPU cores, max 8)

## Robustness check

    python sensitivity.py 1T49                  # contact cutoff 7,8,9 x energy scale 1,3,5
    python sensitivity.py 1T49 --control        # the same with the control (4x slower)

Reruns the protein over the grid and writes a table (CSV and Markdown) of the
distance-adjusted AUC, best gamma, whether noise helps, the p-value and the
classical baseline for every setting.

## Speed and exactness

- With no noise the walk is solved exactly by diagonalising H.
- With noise it uses a sparse H, one mixed-state run for many source residues
  (exact, because the equation is linear), the noise levels in parallel across
  CPU cores, and only the populations are stored (memory grows as N^2).
- `tests/test_walk_core.py` checks this against the original dense solver
  (agreement to about 1e-15 with noise; the old solver's own error without it).

Typical times on 4 cores with the default 5-seed control: 4OBE chain A (169
residues) about 2 min, 1IWH (141) about 70 s; with `--no-control` about 4 times
less. Time grows roughly with the square of the residue count.

## Notes

- Crystal structures often hold several copies of the protein. Use `--chains` to
  keep one, or the walk leaks across crystal contacts (4OBE has two KRAS copies).

## Labels and citation

`data/allo_labels.csv` is Supplementary Table S2 of Wu N, Strömich L, Yaliraki SN,
"Prediction of allosteric sites and signaling: insights from benchmarking
datasets", *Patterns* 3(1), 100408 (2022), doi:10.1016/j.patter.2021.100408
(CC BY 4.0), curated from the Allosteric Database. Cite it if you use it.

## Files

    run_protein.py   the one command; analyze() is shared with the web app
    walk_core.py     the physics: walk, metrics, parallel sweeps
    labels.py        active-site / allosteric residues (ALLO table or your own)
    rin_builder.py   structure -> residue network
    web/             the local web app
    tests/           python tests/test_walk_core.py, python tests/test_labels.py
