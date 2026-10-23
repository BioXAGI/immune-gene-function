# immune-gene-function

Infer whether a gene is linked to T-cell dysfunction, T-cell exclusion, or a specific immune/stromal
cell state from bulk tumor RNA-seq (TCGA via UCSC Xena, or your own data). Per cohort and pan-cancer.

- `immune_gene_function.py`: main pipeline. Signature scores, raw/partial correlations, a TIDE-style
  CTL x gene Cox interaction test (Jiang et al., Nat Med 2018), meta-analysis, plots.
  Run `python immune_gene_function.py --describe` for a JSON schema of methods and outputs.
- `state_context.py`: links EcoTyper cell states (or any per-sample score table) to immune/stromal context and
  to other cell states, pooled across TCGA cohorts with 95% confidence intervals and an ecotype-consistency
  check. Survival-free and standalone (numpy, pandas, scipy only). Also tests user hypothesis signatures
  (`--signature`, `--genes`) and precomputed scores such as ssGSEA/GSVA (`--signature-scores`, compared with the
  mean z-score by default). `python state_context.py describe` prints its JSON schema.
- `bayesprism/bayesprism_benchmark.R`: standalone BayesPrism benchmark (Linux) using colorectal single-cell
  pseudo-bulks with known fractions and cell-type expression. See `bayesprism/README.md`.
- `AGENT_GUIDE.md`: usage and interpretation rules, written for other AI agents and models.

## Quick start
```
pip install -r requirements.txt
python immune_gene_function.py --cohorts COADREAD LUAD SKCM --genes IDO1 TGFB1 --msi
python immune_gene_function.py --scan --cohorts BRCA HNSC SKCM --msi            # genome-wide CTL x gene scan
python state_context.py build --ecotyper-dir <EcoTyper>/Carcinoma/Carcinoma_Fractions/Cell_States/discovery
python state_context.py run --query-prefix Epi_ --msi                          # cell state -> immune context
```
First TCGA run downloads Xena files into `xena_cache/` (about 1 GB for all cohorts).

## Caveats
TCGA cohorts have few events, scores are marker-gene averages (not cell counts), prognosis in untreated
tumors is not immunotherapy response, and a single-cohort p near 0.05 is weak evidence. Include negative
controls (ACTB, GAPDH) and rely on meta-analysis and replication.
