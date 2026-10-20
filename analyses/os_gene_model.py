"""Overall-survival models for a small gene set (default ERAP1 ERAP2 CALR PDIA3) across TCGA cohorts.

Models per cohort (all continuous z-scored genes, hazard ratio per 1 SD):
  M0 univariable: one gene at a time
  M1 clinical:    all genes + stage + age [+ MSI-H]
  M2 immune:      M1 + CTL + CAF
  S  score:       mean z of the genes (sign-free) + clinical covariates, and S + CTL + CAF
Reports HR, 95% CI, p, events, concordance index, and a random-effects meta-analysis across cohorts.
Usage: python os_gene_model.py [--genes G1 G2 ..] [--cohorts C1 C2 ..] [--endpoint OS] [--out DIR]
"""
import argparse
import os
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

import immune_gene_function as m

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--genes", nargs="+", default=["ERAP1", "ERAP2", "CALR", "PDIA3"])
ap.add_argument("--cohorts", nargs="+", default=["COADREAD", "LUAD", "SKCM", "PAAD"])
ap.add_argument("--endpoint", default="OS")
ap.add_argument("--out", default="os_gene_model")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
EP, T = a.endpoint, f"{a.endpoint}.time"


def fit(d, terms):
    cph = CoxPHFitter(penalizer=0.01).fit(d[[T, EP] + terms], T, EP)
    return cph, cph.summary


rows = []
for c in a.cohorts:
    df = m.add_scores(m.load_tcga(c, a.genes, "xena_cache"))
    msi = ["MSI_H"] if df["MSI_H"].notna().mean() >= 0.7 and df["MSI_H"].nunique() > 1 else []
    clin = m.usable_covariates(df, ["stage", "age"] + msi)
    d = df[[T, EP, "CTL", "CAF"] + a.genes + clin].dropna()
    d = d[d[T] > 0].copy()
    for col in a.genes + ["CTL", "CAF"]:
        d[col] = m.zscore(d[col])
    d["SCORE"] = m.zscore(d[a.genes].mean(axis=1))
    ev = int(d[EP].sum())

    def add(model, summ, cph, names):
        for n in names:
            r = summ.loc[n]
            rows.append(dict(cohort=c, model=model, term=n, n=len(d), events=ev, HR=r["exp(coef)"],
                             CI_low=r["exp(coef) lower 95%"], CI_high=r["exp(coef) upper 95%"],
                             logHR=r["coef"], se=r["se(coef)"], p=r["p"], C_index=cph.concordance_index_))

    for g in a.genes:
        cph, s = fit(d, [g] + clin)
        add("M0 gene + clinical", s, cph, [g])
    cph, s = fit(d, a.genes + clin)
    add("M1 all genes + clinical", s, cph, a.genes)
    cph, s = fit(d, a.genes + clin + ["CTL", "CAF"])
    add("M2 all genes + clinical + CTL + CAF", s, cph, a.genes + ["CTL", "CAF"])
    cph, s = fit(d, ["SCORE"] + clin)
    add("S score + clinical", s, cph, ["SCORE"])
    cph, s = fit(d, ["SCORE"] + clin + ["CTL", "CAF"])
    add("S score + clinical + CTL + CAF", s, cph, ["SCORE"])
    cph0, _ = fit(d, clin)
    rows.append(dict(cohort=c, model="clinical only", term="(baseline)", n=len(d), events=ev,
                     C_index=cph0.concordance_index_))
    cc = d[a.genes].corr(method="spearman")
    print(f"[{c}] n={len(d)} events={ev} clinical={clin}; gene-gene Spearman range "
          f"{cc.values[np.triu_indices(len(a.genes), 1)].min():+.2f}..{cc.values[np.triu_indices(len(a.genes), 1)].max():+.2f}")

res = pd.DataFrame(rows)
res.to_csv(os.path.join(a.out, "cohort_models.csv"), index=False, float_format="%.4g")

meta_rows = []
for (model, term), g in res.dropna(subset=["se"]).groupby(["model", "term"]):
    meta_rows.append(dict(model=model, term=term, **m.meta(g["logHR"], g["se"])))
meta = pd.DataFrame(meta_rows)
meta.to_csv(os.path.join(a.out, "meta_analysis.csv"), index=False, float_format="%.4g")

pd.set_option("display.width", 200)
for model in ["M0 gene + clinical", "M1 all genes + clinical", "M2 all genes + clinical + CTL + CAF",
              "S score + clinical", "S score + clinical + CTL + CAF"]:
    print(f"\n== {model}: HR per SD (p) by cohort | meta-analysis ==")
    sub = res[res.model == model]
    for t in sub.term.unique():
        x = sub[sub.term == t].set_index("cohort")
        mm = meta[(meta.model == model) & (meta.term == t)]
        mt = (f"META {mm.HR_random.iloc[0]:.2f} ({mm.CI_low.iloc[0]:.2f}-{mm.CI_high.iloc[0]:.2f}) "
              f"p={mm.p_random.iloc[0]:.3f} I2={mm.I2.iloc[0]:.0%}") if len(mm) and "HR_random" in mm and pd.notna(mm.HR_random.iloc[0]) else ""
        print(f"  {t:6s} " + "  ".join(f"{k}: {x.loc[k, 'HR']:.2f} ({x.loc[k, 'p']:.2f})" for k in x.index) + f" | {mt}")
print("\n== Concordance index ==")
print(res.pivot_table(index="model", columns="cohort", values="C_index", aggfunc="first").round(3).to_string())
print("\nevents:", res.groupby("cohort")["events"].first().to_dict())
