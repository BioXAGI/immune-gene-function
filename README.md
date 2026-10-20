# immune-gene-function

Infer whether a gene is linked to T-cell dysfunction, T-cell exclusion, or a specific immune/stromal
cell state from bulk tumor RNA-seq (TCGA via UCSC Xena, or your own data). Per cohort and pan-cancer.

- `immune_gene_function.py`: main pipeline. Signature scores, raw/partial correlations, a TIDE-style
  CTL x gene Cox interaction test (Jiang et al., Nat Med 2018), meta-analysis, plots.
  Run `python immune_gene_function.py --describe` for a JSON schema of methods and outputs.
- `AGENT_GUIDE.md`: usage and interpretation rules, written for other AI agents and models.
- `analyses/`: follow-up scripts for specific questions (OS gene models, gene-gene synergy, TNF/BIRC/KRAS,
  NOTCH, STK11/PTGES, etc.). They import `immune_gene_function`, so run them from the repo root with
  `PYTHONPATH=.`. Scripts using cBioPortal need `SSL_CERT_FILE` set to certifi's path on some Windows setups.
  `notch_tp53.py`, `stk11_ptges_luad.py`, `klk11_tp53.py`, `lrrc33_tp53.py`, `gene_panel_*.py` are earlier
  one-off scripts with hard-coded `output/...` paths; treat them as legacy.
- `literature_search/`: PubMed and ClinicalTrials.gov query scripts.

## Quick start
```
pip install -r requirements.txt
python immune_gene_function.py --cohorts COADREAD LUAD SKCM --genes IDO1 TGFB1 --msi
```
First TCGA run downloads Xena files into `xena_cache/` (about 1 GB for all cohorts).

## Caveats
TCGA cohorts have few events, scores are marker-gene averages (not cell counts), prognosis in untreated
tumors is not immunotherapy response, and a single-cohort p near 0.05 is weak evidence. Include negative
controls (ACTB, GAPDH) and rely on meta-analysis and replication.
