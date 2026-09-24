# Pre-registration of the replication set

Written and committed before any protein below was run. The commit that adds
this file is the time stamp.

## How the protein sets came about (the full history, including what was not pre-specified)

1. **Development set** (1T49, 3LSW, 1IWH, 3CSM): picked by hand while the method
   was built. Not pre-specified.
2. **Validation set, first stage** (2RD5, 3HO6, 4B1F, 4BBG, 4PFK, 3PYY): the rule in
   `select_proteins.py` and the draw in `selection.csv` were committed (607c3b6)
   before any of them was run.
3. **Validation set, second stage** (3ZCW, 3HFR, 3M3F, 3H30): the next four in the
   same fixed order. The decision to extend from six to ten was taken *after* the
   first six results were seen (commit 37492fb, after the six runs). The order was
   fixed in advance; the stopping point was not.

## Replication set (this file)

The remaining ten entries of the fixed draw order (`selection.csv`, rows 11-20) are
mostly other structures of proteins already studied. To add only new proteins, the
replication set is, in draw order, the first structure of every protein (by the ALLO
`protein` name) that does not already appear in the development or validation sets:

| draw order | PDB | protein |
|---|---|---|
| 11 | 2YHD | Androgen receptor |
| 12 | 4HO6 | Glucose-1-phosphate thymidylyltransferase |
| 15 | 2VD3 | ATP phosphoribosyltransferase |
| 16 | 3PXF | Cyclin-dependent kinase 2 |
| 17 | 1ZDS | Copper-containing nitrite reductase |

Skipped as repeats: 2YLO and 2QPY (androgen receptor, after 2YHD), 3K5V (ABL1, as
3PYY), 1NH8 (ATP phosphoribosyltransferase, after 2VD3), 2VVT (glutamate racemase,
as 4B1F and 3HFR). This exhausts the eligible pool: **this set is final and will not
be extended**, whatever it shows.

## Analysis plan

- Settings: the pipeline defaults, unchanged (8 Å Cβ contacts, hydropathy energies,
  s = 3, 18 noise levels, T = 30, 5 random-energy seeds, 10,000 shuffles, seed 0).
- Primary outcome per protein: the best AUC beyond distance over γ and its
  within-shell permutation p.
- Counting: by distinct protein. A protein with several structures counts once,
  and is significant if any of its structures is; the protein that repeats a
  development protein (3M3F, as 3LSW) is not counted among the held-out proteins.
- Reported for the replication set alone, and for all held-out proteins together
  (validation + replication): the number significant at p < 0.05, the one-sided
  binomial probability of at least that many by chance, and the number passing a
  Benjamini-Hochberg correction at 5%.
- Every protein is reported, including failures.
