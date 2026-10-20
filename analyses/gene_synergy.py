"""Pairwise (and three-way) gene x gene interaction on survival, per TCGA cohort + meta-analysis.
Cox: hazard ~ A + B + A*B + clinical [+ CTL + CAF]; genes z-scored. Interaction HR > 1 means the joint
effect is larger than the product of the separate effects (hazard-scale synergy); < 1 means antagonism /
buffering. Usage: python gene_synergy.py [--genes BIRC2 TNF KRAS] [--cohorts ...] [--endpoint OS]"""
import argparse, itertools, os, warnings
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
import immune_gene_function as m
warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--genes", nargs="+", default=["BIRC2", "TNF", "KRAS"])
ap.add_argument("--cohorts", nargs="+", default=["COADREAD", "LUAD", "SKCM", "PAAD"])
ap.add_argument("--endpoints", nargs="+", default=["OS", "PFI"])
ap.add_argument("--out", default="gene_synergy")
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)

rows = []
for c in a.cohorts:
    df = m.add_scores(m.load_tcga(c, a.genes, "xena_cache"))
    msi = ["MSI_H"] if df["MSI_H"].notna().mean() >= 0.7 and df["MSI_H"].nunique() > 1 else []
    clin = m.usable_covariates(df, ["stage", "age"] + msi)
    for ep in a.endpoints:
        T = f"{ep}.time"
        d = df[[T, ep, "CTL", "CAF"] + a.genes + clin].dropna(); d = d[d[T] > 0].copy()
        for g in a.genes + ["CTL", "CAF"]:
            d[g] = m.zscore(d[g])
        for adj_name, adj in [("clinical", clin), ("clinical+CTL+CAF", clin + ["CTL", "CAF"])]:
            for A, B in itertools.combinations(a.genes, 2):
                d[f"{A}x{B}"] = d[A] * d[B]
                terms = [A, B, f"{A}x{B}"] + adj
                s = CoxPHFitter(penalizer=0.01).fit(d[[T, ep] + terms], T, ep).summary
                for t in [A, B, f"{A}x{B}"]:
                    r = s.loc[t]
                    rows.append(dict(cohort=c, endpoint=ep, adjust=adj_name, pair=f"{A}+{B}", term=t,
                                     kind="interaction" if "x" in t and t == f"{A}x{B}" else "main",
                                     events=int(d[ep].sum()), HR=r["exp(coef)"], logHR=r["coef"],
                                     se=r["se(coef)"], p=r["p"]))
            if len(a.genes) == 3 and adj_name == "clinical":
                A, B, C = a.genes
                d["three"] = d[A] * d[B] * d[C]
                terms = a.genes + [f"{A}x{B}", f"{A}x{C}", f"{B}x{C}", "three"] + adj
                d[f"{A}x{C}"] = d[A] * d[C]; d[f"{B}x{C}"] = d[B] * d[C]
                s = CoxPHFitter(penalizer=0.01).fit(d[[T, ep] + terms], T, ep).summary
                r = s.loc["three"]
                rows.append(dict(cohort=c, endpoint=ep, adjust=adj_name, pair="three-way", term="three",
                                 kind="interaction", events=int(d[ep].sum()), HR=r["exp(coef)"],
                                 logHR=r["coef"], se=r["se(coef)"], p=r["p"]))
res = pd.DataFrame(rows); res.to_csv(os.path.join(a.out, "cohort_interactions.csv"), index=False, float_format="%.4g")
meta = []
for (ep, adj, pair, term, kind), g in res.groupby(["endpoint", "adjust", "pair", "term", "kind"]):
    meta.append(dict(endpoint=ep, adjust=adj, pair=pair, term=term, kind=kind, **m.meta(g["logHR"], g["se"])))
meta = pd.DataFrame(meta); meta.to_csv(os.path.join(a.out, "meta_analysis.csv"), index=False, float_format="%.4g")

for ep in a.endpoints:
    print(f"\n===== {ep}: interaction HR (p) by cohort | random-effects meta =====")
    for adj in ["clinical", "clinical+CTL+CAF"]:
        print(f"-- adjusted for {adj}")
        x = res[(res.endpoint == ep) & (res.adjust == adj) & (res.kind == "interaction")]
        for pair in x.pair.unique():
            y = x[x.pair == pair].set_index("cohort")
            mm = meta[(meta.endpoint == ep) & (meta.adjust == adj) & (meta.pair == pair) & (meta.kind == "interaction")]
            mt = f"META {mm.HR_random.iloc[0]:.2f} ({mm.CI_low.iloc[0]:.2f}-{mm.CI_high.iloc[0]:.2f}) p={mm.p_random.iloc[0]:.3f} I2={mm.I2.iloc[0]:.0%}" if len(mm) and "HR_random" in mm and pd.notna(mm.HR_random.iloc[0]) else ""
            print(f"  {pair:11s} " + "  ".join(f"{k}: {y.loc[k,'HR']:.2f} ({y.loc[k,'p']:.3f})" for k in y.index) + f" | {mt}")
