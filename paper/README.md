# Figures and tables of the main text

This folder reproduces the figures and tables of the main text of the paper from the source tables in `data/`.

```bash
pip install matplotlib pandas numpy
python paper/scripts/make_figure1.py          # Fig. 1
python paper/scripts/make_figures_tables.py   # Figs. 2-4, Tables 2-3
python paper/scripts/make_table1.py           # Table 1
```

Outputs are written to `paper/figures/` (PNG and PDF) and `paper/tables/` (LaTeX). The copies in this repository
were produced by these commands and are identical to those in the paper.

| Output | Script | Source data |
|---|---|---|
| Fig. 1, method overview | `make_figure1.py` | `figure1_counts.json` (reference sizes) |
| Fig. 2 and Table 2, single-cell benchmark | `make_figures_tables.py` | `single_cell_benchmark.csv` |
| Fig. 3 and Table 3, bulk benchmark | `make_figures_tables.py` | `bulk_accuracy.csv`, `bulk_tests.csv`, `bulk_held_out_tests.csv` |
| Fig. 4, what drives a prediction | `make_figures_tables.py` | `interpretability_genes.csv` |
| Table 1, datasets | `make_table1.py` | `datasets.csv` |

## Source tables

- **`single_cell_benchmark.csv`**: one row per cohort (`gse181897`, `gse327907`), candidate panel (`P10`, `P37`,
  `P90`) and method. Top-1 accuracy per receiver and averaged over the five receivers (`top1_receiver_mean`) with its
  95% donor-cluster bootstrap interval (`ci_lo`, `ci_hi`). `n_receivers_r2s_hier_better` / `_worse`: the number of
  receivers in which R2S-hierarchical was significantly more / less accurate than the method (donor-paired sign-flip
  test, Holm-adjusted p < 0.05).
- **`bulk_accuracy.csv`**: study-mean top-1 accuracy with 95% study-bootstrap intervals, per bulk set and method.
  Sets: `reference_loso` (the CytoSig source studies, scored leave-one-study-out), `held_out_1`, `held_out_2`,
  `held_out_3` (studies outside the reference). `panel`: `full` (41 signals) or `CytoSig sub-panel` (30 signals).
- **`bulk_tests.csv`**: study-paired differences in accuracy (`diff`, 95% interval, sign-flip `p`) between
  R2S-hierarchical or R2S-pool and a comparator.
- **`bulk_held_out_tests.csv`**: the pre-specified comparisons on the held-out sets (R2S-hierarchical vs R2S-pool on
  41 signals; R2S-hierarchical vs published CytoSig on the 30 signals CytoSig covers). `p_holm` is the
  Holm-adjusted p value of the pre-registered test family on held-out set 1.
- **`interpretability_genes.csv`**: for each case (GSE181897 receiver and applied cytokine; TNF in held-out set 1),
  the 20 genes with the largest mean contribution to the coefficient of the applied signal, and whether each is among
  the 20 largest weights of the signal's consensus signature.
- **`datasets.csv`**, **`figure1_counts.json`**: dataset descriptions and reference sizes.

Cell types: `cd14_mono` (CD14 monocytes), `nk` (NK cells), `cd4_naive`, `cd8_naive` (naive CD4 and CD8 T cells),
`b_naive` (naive B cells).
