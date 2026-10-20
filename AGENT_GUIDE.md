# immune_gene_function.py: guide for AI agents and other models

## What it answers
For a gene X in tumor bulk RNA-seq: is X linked to **T-cell dysfunction** (TIDE-style), to **T-cell exclusion**
(fibroblast / myeloid compartments), or to a specific immune cell state? Works per cohort and pan-cancer.

Quick self-description: `python immune_gene_function.py --describe` (JSON schema of methods, columns, files,
interpretation rules). Every run writes `<out>/run_summary.json` (command, hits, top meta-analysis results).

## Minimal runs
```
python immune_gene_function.py --cohorts COADREAD LUAD SKCM --genes IDO1 TGFB1 --msi
python immune_gene_function.py --cohorts COADREAD --genes IDO1 --no-survival --pcorr-adjust full
python immune_gene_function.py --source custom --expr expr.tsv.gz --clinical clin.tsv --genes IDO1
```
First TCGA run downloads Xena files into `xena_cache/` (needs internet). Requires numpy, pandas, scipy,
lifelines, matplotlib. If Python reports an SSL certificate error on Windows, set `SSL_CERT_FILE` to certifi's path.

## Inputs
- TCGA cohorts by code (COADREAD, LUAD, SKCM, ... or ALL). One primary tumor per patient.
- Custom: `--expr` genes x samples, log-scale, first column gene symbol; optional `--clinical` with time/event
  columns (`--time-col`, `--event-col`) and optional `age`, `stage`, `msi`, `cohort`.
- Gene aliases (e.g. VSIR = C10orf54, LRRC33 = NRROS) are resolved automatically.

## Pipeline (per cohort, per gene X)
1. **Signature scores** = mean z-score of marker genes (see `--describe`, key `scores`). CTL = CD8A CD8B GZMA GZMB PRF1.
   A gene that belongs to a signature is removed from it before correlating (no self-correlation).
2. **Correlations (no outcome needed):** `rho_*` raw; `pcorr_*` partial (leukocyte PTPRC + CAF [+ MSI-H]);
   `pcorr_full_*` partial adjusted for leukocyte, CAF, exhaustion, Treg, macrophage, M2, monocyte, neutrophil,
   MDSC-like, VEGF [+ MSI-H]. `--decompose` shows which covariate removes how much.
3. **Dysfunction test (Cox):** `hazard ~ S + X + S*X + stage + age + CAF [+ MSI-H]`, S = CTL by default.
   Column `<EP>_HR_CTLxX` > 1 with small p: high X weakens the survival benefit of S. EP = OS or PFI.
4. FDR (Benjamini-Hochberg) within cohort; random-effects meta-analysis across cohorts (`meta_analysis.csv`).

## How to read results (rules for conclusions)
- Trust **meta-analysis q_random and cross-cohort consistency (I2)** over any single-cohort p near 0.05.
  Housekeeping genes (ACTB, GAPDH) reach |z| ~ 2 in single cohorts by chance; include them as negative controls.
- Positive controls that behave as expected: IDO1, TGFB1, PDCD1LG2 (CTLxX HR > 1 in meta-analysis).
- `pcorr_full_* >= 0.3` = link survives heavy adjustment (specific). Heavy adjustment can over-correct, so it is a lower bound.
- Myeloid / exhaustion / Treg scores are **not** in the default Cox model because they may be mediators
  (TIDE also left them out). Use `--adjust-scores` for sensitivity; `--interaction-score M2_LIKE` to test a
  different score x gene interaction.
- Prognosis in untreated TCGA is not immunotherapy response. TIDE itself restricts to cohorts where CTL
  predicts better survival; this script does not, so check `<EP>_HR_CTL` in the output before trusting a cohort.
- Correlation is not causation; bulk RNA cannot say which cell expresses a gene.

## Output files (`--out`, default `immune_gene_results/`)
`per_cohort_results.csv`, `meta_analysis.csv`, `decomposition.csv` (--decompose), `score_prognosis.csv`
(--score-prognosis), `heatmap_*.png`, `km_*.png` (--km), `run_summary.json`.

## Column naming
`rho_<SCORE>`, `pcorr_<SCORE>`, `pcorr_full_<SCORE>`, `<EP>_HR_CTLxX|_p_CTLxX|_q_CTLxX`, `<EP>_HR_X|_p_X`,
`<EP>_HR_CTL_Xlow|_Xhigh`, `<EP>_n|_events|_covariates`, `interpretation`, `in_signatures`.
When `--interaction-score` is not CTL, the columns still carry the `CTL` names but refer to that score
(recorded in column `interaction_score`).

## Reference
Jiang P et al. Signatures of T cell dysfunction and exclusion predict cancer immunotherapy response.
Nat Med 24:1550-1558 (2018).
