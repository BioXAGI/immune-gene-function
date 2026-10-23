# BayesPrism benchmark with single-cell ground truth (colorectal cancer)

`bayesprism_benchmark.R` is a self-contained R script (written for Linux) that tests how well
[BayesPrism](https://github.com/Danko-Lab/BayesPrism) recovers cell-type fractions and cell-type-specific
expression from bulk RNA-seq, using pseudo-bulk samples built from held-out patients' single cells.

- Data: Lee et al. 2020 colorectal 10x scRNA-seq (GEO GSE132465), downloaded automatically.
- Patients (not cells) are split into a reference set and a held-out test set.
- Pseudo-bulks: one per held-out sample (real composition) plus random Dirichlet mixtures.
- Scores: fractions (Pearson, Spearman, RMSE, bias vs true UMI and cell fractions, per cell type and state),
  expression (per-gene and per-sample correlation of BayesPrism Z vs true counts, normalized absolute error),
  and a leakage check (share of marker-gene reads assigned to tumor cells).

## Setup and run
```
sudo apt-get install -y r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev   # Debian/Ubuntu
Rscript bayesprism_benchmark.R --install=TRUE --help          # install missing R packages, show options
nohup Rscript bayesprism_benchmark.R --out=bp_benchmark > bp_benchmark.log 2>&1 &
```
Outputs in `--out`: `report.md`, `metrics_fractions.tsv`, `metrics_expression.tsv`, `leakage.tsv`,
truth and estimated fraction tables, `fractions_scatter.png`.

## Example result (quick settings: 8 held-out patients, 19 pseudo-bulks, 2 cores, 37 min)
| cell type | fraction Pearson | expression median per-gene r |
|---|---|---|
| Tumor | 0.99 | 0.98 |
| T cells | 0.99 | 0.84 |
| Myeloids | 0.98 | 0.95 |
| Stromal cells | 0.96 (overestimated, bias +0.07) | 0.93 |
| B cells | 0.99 (underestimated, bias -0.13) | 0.92 |

Immune-gene leakage into the tumor component matched the truth (PTPRC 0.7% estimated vs 0.7% true). This
is a best case: reference and pseudo-bulk share the 10x platform. On real bulk RNA-seq (TCGA) we observed
about 20% of immune-gene reads assigned to tumor cells, which reflects the platform difference.
