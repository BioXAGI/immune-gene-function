# CRC MSI subgroup analysis

`--msi-cohort pooled` (default) retains the pooled expression cohort, including samples
without an MSI annotation. `--msi` separately requests the historical MSI adjustment;
it does **not** select an MSS cohort. Inspect `msi_adjusted`,
`msi_adjustment_status`, and the actual endpoint covariates before interpreting a run.
As before, the generic pooled pipeline only includes MSI adjustment when >=70% have
recognized annotations and both groups are present. The comparison runner below uses
recognized annotations only and requires both groups; it cannot silently become unadjusted.

## Selection and annotation

- `--msi-cohort mss`: MSS plus MSI-L by default. `--msi-low exclude` selects strictly MSS.
- `--msi-cohort msi-h`: only MSI-H. Unknown or unrecognized labels never enter either subgroup.
- `--msi-missing keep` is the pooled default; `exclude` removes unknown annotations;
  `error` refuses a dataset containing any unknown annotation, including in subgroup runs.
- Recognized labels are case-insensitive MSS/stable/microsatellite stable, MSI-L/MSIL/L/low,
  MSI-H/MSIH/H/high, and microsatellite instability-high/low; underscores become hyphens.
  Free text, ambiguous mixed labels, and numeric codes are unknown. Custom data must use
  a string `msi` column with these labels; numeric encodings require explicit recoding.
- The TCGA loader chooses the candidate MSI column with the most recognized annotations
  (ties follow `MSI_COLUMNS` order). It no longer discards annotation columns with <=20
  observations. It records the selected column; it does not reconcile disagreements among
  different clinical annotation sources. Inspect the source before using external MSI data.
- MSI status is not direct proof of mismatch repair protein proficiency or deficiency.
- TCGA subgroup selection is limited to COADREAD/COAD/READ; custom data are assumed to
  be the CRC cohort supplied by the investigator. Targeted selection is not supported by
  `--scan`; combinations that could silently ignore selection are rejected.

Selection happens before signature scoring in the main CLI. The comparison runner holds
pooled-reference signature construction fixed across groups, preserving a common score
reference. This is deliberately distinct from separate main-CLI runs that re-score each subgroup.
The existing loader chooses one tumor per patient, preferring primary tumors but retaining
its existing blood-derived/metastatic fallbacks. Custom input must contain one independent
sample per patient. The comparison runner rejects duplicate sample identifiers.

## Reproduce

Use Python 3.12. Run from the repository root:

```bash
python -m venv .venv
# Linux/macOS:
. .venv/bin/activate
# Windows PowerShell instead: .venv\Scripts\Activate.ps1
python -m pip install -r analysis/requirements.txt
python -m unittest discover -s tests -v

# Original pooled behavior, with requested MSI adjustment:
python immune_gene_function.py --cohorts COADREAD --genes CD276 C1GALT1 --msi --out results/pooled
# Explicit MSS (including MSI-L); no varying MSI covariate in this subgroup:
python immune_gene_function.py --cohorts COADREAD --genes CD276 C1GALT1 --msi-cohort mss --out results/mss
# Strict MSS sensitivity:
python immune_gene_function.py --cohorts COADREAD --genes CD276 C1GALT1 --msi-cohort mss --msi-low exclude --out results/strict_mss
# MSI-H-only exploratory cohort; sparse survival fits may be skipped:
python immune_gene_function.py --cohorts COADREAD --genes CD276 C1GALT1 --msi-cohort msi-h --out results/msih

# Complete predeclared panel, basic/full correlations and OS/PFI CTL interactions:
python compare_crc_msi.py --config analysis/crc_glycosylation.json --out results/crc_msi
# Limited correlation-only run, retaining model-free uncertainty and BH correction:
python compare_crc_msi.py --config analysis/crc_glycosylation.json --no-survival --bootstrap 100 --out results/crc_pilot
# Supply a local CRC log-expression matrix and clinical annotations instead:
python compare_crc_msi.py --config analysis/crc_glycosylation.json --expr expr.tsv.gz --clinical clinical.tsv --out results/crc_custom
```

Expression is genes x samples, log scale; clinical is indexed by the same sample IDs,
with `msi`, optional `age` and `stage`, and `OS.time`/`OS`, `PFI.time`/`PFI` when available.
Times must be positive finite values and events 0/1. No outcome, MSI, or covariate imputation
is performed. Custom sample IDs must already identify distinct patients.

The configuration fixes the 19-gene hypothesis/control panel, eight signatures, both
adjustment sets, endpoints, seed, 1,000 bootstrap replicates, and model thresholds. For
strict-MSS sensitivity, copy the configuration and change `msi_low` to `exclude`; the
pooled reference still includes MSI-L. All panel genes, including controls, enter the
same declared FDR families. ACTB/GAPDH are calibration genes, not guaranteed nulls.

## Outputs and statistical contract

The main CLI writes `cohort_selection.json` before fitting, even for an empty subgroup.
It reports input/selected/excluded counts, label counts, source, and available endpoint
sample/event counts. Per-gene `OS_n`/`OS_events` etc. are model-specific complete-case
counts and may be smaller. Each attempted Cox model has `status`, number of parameters,
actual covariates, and an interaction HR/95% CI only if fitting succeeds.

Main-CLI correlation outputs include coefficient, p, approximate Fisher-z 95% CI, actual
n, status and adjustment variables. The residual-correlation t test uses n-k-2 degrees
of freedom, correcting the previous residual Pearson p-value that ignored k covariates.
This changes historical partial-correlation p-values; effect estimates for nonsingular
models are unchanged. Fisher-z intervals and t p-values are approximations for ranked data.
BH q-values are computed per cohort and adjustment set across the reported genes x all
reported signatures; the comparison runner instead uses its explicitly predeclared panel.

The comparison runner produces:

| File | Meaning |
| --- | --- |
| `estimates.csv` | MSS and pooled MSI-adjusted partial correlations and CTL interaction log HR/HR, CIs, p/q, n/events, model status |
| `mss_minus_pooled.csv` | MSS minus pooled estimates; paired bootstrap CI, valid/requested resamples, status |
| `cohort_selection.json` | Annotation and selection audit plus endpoint event counts, including MSI-H counts |
| `run_manifest.json` | Resolved configuration, effective flags, seed, versions, source/input SHA256, git commit and completion/failure status |

BH families in the comparison are (1) all 19 genes x 8 signatures per adjustment set and
cohort; (2) all 19 CTL interactions per endpoint and cohort. Failed or missing hypotheses
count as p=1 in the denominator and receive no reported q. OS and PFI are separate
exploratory endpoint families, not a single confirmatory family. There is no additional
cross-family or cross-cohort multiplicity claim. Full adjustment may include mediators
or induce collider bias; it is a sensitivity model, not proof of a causal direct effect.

The pooled comparator includes recognized MSI annotations only, eliminating a silent
unadjusted fallback. Both Cox fits use the same available covariate set (>=70% coverage in
the annotated pool) and its common complete cases. Missing or constant mandatory MSI
adjustment prevents the pooled fit. Gene and CTL standardization is within each fit as
in the original pipeline, then interaction coefficients/SEs are converted to the same
pooled gene-SD x CTL-SD scale for the comparison. The existing ridge penalty 0.01 remains;
penalized Wald intervals/p-values are exploratory. Different predictor variances can
still change the effective regularization; test penalty sensitivity before relying on
small differences. Age/stage/CAF are handled as in the original model (stage as ordinal).

The MSS subset and pooled reference overlap. Do not use an independent Fisher-z test,
compare p-value significance labels, or meta-analyze them as independent cohorts.
Paired MSI-stratified patient resampling fits both models to the same draw, preserving
that dependence and observed MSI composition. Signature construction and pooled SD units
are held fixed; rank residuals and Cox fits are refit in each draw. Percentile delta CIs
are **descriptive**, with no delta p/q or formal between-group claim. At least 100 valid
replicates and >=80% of requested replicates are required to report a delta CI. Use the
configured 1,000 or more for interpretation; 100 is only a limited execution check.
No interval is emitted when those conditions fail; `--bootstrap 0` is point estimates only.

Cox guards default to >=40 complete cases, >=20 events and >=10 events per fitted
parameter (including gene, CTL, interaction and all adjustment columns). Constant
predictors, singular/ill-conditioned designs, convergence warnings, exceptions, and
nonfinite/invalid SEs suppress estimates and produce a status. Split-half descriptive
models have separate counts and guards. These are pragmatic screens, not evidence that
proportional hazards or causal assumptions hold. Review PH diagnostics, influential
observations, sparse-event imbalance, and external replication before interpretation.
These new guards cover the targeted CTL interaction path; older optional score-prognosis
and genome-wide scan models retain their existing separate methods.

Bulk RNA cannot identify the expressing cell, glycosylation activity, membrane protein
abundance, tumor recognition, or drug target dependence. TCGA survival is prognosis,
not immunotherapy response. The attached DeLLphi-305 trial slide is not an input to this
analysis and does not validate these CRC associations.
