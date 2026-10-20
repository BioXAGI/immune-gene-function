"""Effect of BIRC-family genes on OS/PFI inside TNF-high vs TNF-low tumors (median split per cohort).
Cox per stratum: hazard ~ gene (z-scored within stratum) + stage + age [+ MSI-H].
Formal test: hazard ~ gene + TNF + gene*TNF + covariates on the full cohort (continuous TNF).
Meta-analysis across cohorts per stratum and for the interaction."""
import os, warnings
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
import immune_gene_function as m
warnings.filterwarnings("ignore")
GENES = ["BIRC2", "BIRC3", "XIAP", "BIRC5"]
COHORTS = ["COADREAD", "LUAD", "SKCM", "PAAD"]
OUT = "tnf_strata_birc"; os.makedirs(OUT, exist_ok=True)

rows = []
for c in COHORTS:
    df = m.load_tcga(c, GENES + ["TNF"], "xena_cache")
    miss = [g for g in GENES + ["TNF"] if g not in df]
    if miss: print(c, "missing", miss)
    msi = ["MSI_H"] if df["MSI_H"].notna().mean() >= 0.7 and df["MSI_H"].nunique() > 1 else []
    clin = m.usable_covariates(df, ["stage", "age"] + msi)
    for ep in ["OS", "PFI"]:
        T = f"{ep}.time"
        for g in [x for x in GENES if x in df]:
            d = df[[T, ep, g, "TNF"] + clin].dropna(); d = d[d[T] > 0].copy()
            hi = d["TNF"] >= d["TNF"].median()
            for lab, sub in [("TNF-low", d[~hi]), ("TNF-high", d[hi])]:
                sub = sub.copy(); sub["X"] = m.zscore(sub[g])
                s = CoxPHFitter(penalizer=0.01).fit(sub[[T, ep, "X"] + clin], T, ep).summary.loc["X"]
                rows.append(dict(cohort=c, endpoint=ep, gene=g, group=lab, n=len(sub), events=int(sub[ep].sum()),
                                 HR=s["exp(coef)"], CI_low=s["exp(coef) lower 95%"], CI_high=s["exp(coef) upper 95%"],
                                 logHR=s["coef"], se=s["se(coef)"], p=s["p"]))
            z = d.copy(); z["X"], z["Tn"] = m.zscore(z[g]), m.zscore(z["TNF"]); z["XxT"] = z["X"] * z["Tn"]
            s = CoxPHFitter(penalizer=0.01).fit(z[[T, ep, "X", "Tn", "XxT"] + clin], T, ep).summary.loc["XxT"]
            rows.append(dict(cohort=c, endpoint=ep, gene=g, group="interaction", n=len(z), events=int(z[ep].sum()),
                             HR=s["exp(coef)"], CI_low=s["exp(coef) lower 95%"], CI_high=s["exp(coef) upper 95%"],
                             logHR=s["coef"], se=s["se(coef)"], p=s["p"]))
res = pd.DataFrame(rows); res.to_csv(f"{OUT}/cohort_results.csv", index=False, float_format="%.4g")
meta = pd.DataFrame([dict(endpoint=ep, gene=g, group=gr, **m.meta(x["logHR"], x["se"]))
                     for (ep, g, gr), x in res.groupby(["endpoint", "gene", "group"])])
meta.to_csv(f"{OUT}/meta_analysis.csv", index=False, float_format="%.4g")

for ep in ["OS", "PFI"]:
    print(f"\n===== {ep}: HR per SD of BIRC gene (p) =====")
    for g in GENES:
        for gr in ["TNF-low", "TNF-high", "interaction"]:
            x = res[(res.endpoint == ep) & (res.gene == g) & (res.group == gr)].set_index("cohort")
            if x.empty: continue
            mm = meta[(meta.endpoint == ep) & (meta.gene == g) & (meta.group == gr)].iloc[0]
            mt = f"META {mm.HR_random:.2f} ({mm.CI_low:.2f}-{mm.CI_high:.2f}) p={mm.p_random:.3f} I2={mm.I2:.0%}" if pd.notna(mm.get("HR_random")) else ""
            print(f"  {g:6s} {gr:11s} " + "  ".join(f"{k}: {x.loc[k,'HR']:.2f} ({x.loc[k,'p']:.2f})" for k in x.index) + f" | {mt}")
print("\nn / events per stratum (TNF-low | TNF-high):")
for ep in ["OS"]:
    for c in COHORTS:
        x = res[(res.endpoint == ep) & (res.gene == "BIRC2") & (res.cohort == c)].set_index("group")
        print(f"  {c}: low n={x.loc['TNF-low','n']} ev={x.loc['TNF-low','events']} | high n={x.loc['TNF-high','n']} ev={x.loc['TNF-high','events']}")
