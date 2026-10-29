# Limited TCGA COADREAD correlation run, 2026-10-29

Executed with Python 3.12, the committed 19-gene/8-signature configuration, basic and full
adjustments, seed 20261029, and 100 paired bootstrap replicates. The downloaded Xena
HiSeqV2 matrix supplied 380 selected tumors: MSS 254, MSI-L 63, MSI-H 53, unknown 10.
MSS includes MSI-L (n=317); the pooled MSI-adjusted comparator has 370 annotated tumors.
These counts describe this Xena matrix, not every TCGA colorectal sample in other releases.

Exact command (output location here differs from the scratch execution path in the manifest):

```bash
python compare_crc_msi.py --config analysis/crc_glycosylation.json --no-survival --bootstrap 100 --out analysis/results/crc_msi_pilot
```

All 608 correlation estimates passed the numerical checks. All 304 paired differences had
100/100 valid bootstrap replicates. Individual CIs use the approximate Fisher-z method;
delta CIs use paired percentile bootstrap and are descriptive. The bootstrap count is
small; reproduce with the default 1,000 before interpreting interval boundaries.

For illustration, basic-adjusted CTL correlations were:

| Gene | MSS r (95% CI) | Pooled MSI-adjusted r (95% CI) | MSS minus pooled (bootstrap 95% CI) |
| --- | --- | --- | --- |
| CD276 | 0.013 (-0.098, 0.123) | 0.017 (-0.086, 0.119) | -0.004 (-0.039, 0.031) |
| C1GALT1 | 0.005 (-0.105, 0.116) | -0.026 (-0.128, 0.077) | 0.031 (0.001, 0.075) |

The C1GALT1 basic-CTL delta interval is close to zero and is one of 304 descriptive
comparisons without delta multiple-testing inference. It is not evidence for an MSS-specific
dependency. Full-adjusted results differ: read both model families and their covariates in
the CSVs rather than selecting a favorable analysis. This is exploratory bulk-RNA evidence.

OS/PFI availability before model covariate filtering:

| Cohort | OS n / events | PFI n / events |
| --- | --- | --- |
| Annotated pooled | 362 / 84 | 360 / 99 |
| MSS + MSI-L | 311 / 72 | 309 / 87 |
| MSI-H | 51 / 12 | 51 / 12 |

No survival model was fitted: `lifelines` was unavailable and package installation was
blocked in the execution environment. No CTL interaction estimates or immunotherapy-response
conclusions are claimed. MSI-H has only 12 events per endpoint and would fail the default
>=20-event guard even before the events-per-parameter guard. Full reproduction commands,
the direct dependency pins, model safeguards, and interpretation rules are in
`docs/CRC_MSI_ANALYSIS.md`. Input and source hashes are in `run_manifest.json`.

Synthetic validation: 26 tests collected, 24 passed; two actual lifelines-fit tests skipped
because the dependency was unavailable. Successful-summary and convergence-warning output
paths were additionally exercised with explicit mocked fitters. The GitHub workflow installs
the dependency and runs the full synthetic suite; its result must be checked separately.
