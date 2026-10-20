"""CTL x gene Cox interaction with and without IFN-gamma adjustment (ADAR vs controls)."""
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

import immune_gene_function as m

warnings.filterwarnings("ignore")
GENES = ["ADAR", "IDO1", "ACTB", "GAPDH"]
COHORTS = ["COADREAD", "LUAD", "SKCM"]


def fit(d, time, event, x, extra):
    d = d.copy()
    d["CTL"] = m.zscore(d["CTL"])
    d["X"] = m.zscore(d[x])
    d["CTLxX"] = d["CTL"] * d["X"]
    terms = ["CTL", "X", "CTLxX"] + extra
    s = CoxPHFitter(penalizer=0.01).fit(d[[time, event] + terms], time, event).summary
    return s.loc["CTLxX", "exp(coef)"], s.loc["CTLxX", "p"]


rows = []
for c in COHORTS:
    df = m.add_scores(m.load_tcga(c, GENES, "xena_cache"))
    msi = ["MSI_H"] if df["MSI_H"].notna().mean() >= 0.7 else []
    cov = m.usable_covariates(df, ["stage", "age", "CAF"] + msi)
    df["IFNG_z"] = m.zscore(df["IFNG"])
    df["CTLxIFN"] = m.zscore(df["CTL"]) * df["IFNG_z"]
    for g in GENES:
        # IFN-independent part of the gene (IFNG score removed from IDO1 would be circular,
        # so the IFNG signature used here excludes the gene itself)
        ifn = m.score_without(df, "IFNG", g)
        df["IFN_g"] = m.zscore(ifn)
        df["CTLxIFN_g"] = m.zscore(df["CTL"]) * df["IFN_g"]
        ok = df[[g, "IFN_g"]].dropna()
        b = np.polyfit(ok["IFN_g"], ok[g], 1)
        df["resid"] = df[g] - np.polyval(b, df["IFN_g"])
        rho = df[[g, "IFN_g"]].corr(method="spearman").iloc[0, 1]
        for ep in ["OS", "PFI"]:
            t = f"{ep}.time"
            d = df[[t, ep, "CTL", g, "resid", "IFN_g", "CTLxIFN_g"] + cov].dropna()
            d = d[d[t] > 0]
            specs = {"base": (g, []), "+IFN": (g, ["IFN_g"]),
                     "+IFN+CTLxIFN": (g, ["IFN_g", "CTLxIFN_g"]), "IFN-independent": ("resid", cov)}
            for name, (x, extra) in specs.items():
                extra = [e for e in extra if e not in cov] + (cov if name != "IFN-independent" else [])
                hr, p = fit(d, t, ep, x, extra if name != "IFN-independent" else cov + ["IFN_g"])
                rows.append(dict(cohort=c, gene=g, endpoint=ep, model=name, HR=hr, p=p, rho_IFN=rho))
r = pd.DataFrame(rows)
r.to_csv("adar_ifn_check.csv", index=False, float_format="%.4g")
for ep in ["OS", "PFI"]:
    print(f"\n== {ep}: CTLxX interaction HR (p) ==")
    for g in GENES:
        for c in COHORTS:
            x = r[(r.gene == g) & (r.cohort == c) & (r.endpoint == ep)].set_index("model")
            print(f"  {g:5s} {c:9s} rho(IFN)={x.rho_IFN.iloc[0]:+.2f} | " +
                  " | ".join(f"{k} {x.loc[k, 'HR']:.2f} ({x.loc[k, 'p']:.2f})" for k in x.index))
