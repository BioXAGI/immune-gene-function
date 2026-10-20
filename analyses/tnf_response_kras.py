"""TNF response proxies by KRAS mutation status (TCGA RNA-seq + cBioPortal mutation calls)."""
import json, os, time, urllib.error, urllib.request, warnings
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
from scipy import stats
import immune_gene_function as m
warnings.filterwarnings("ignore")

COHORTS = {"COADREAD": "coadread", "LUAD": "luad", "PAAD": "paad"}
TNF_TARGETS = ["NFKBIA", "TNFAIP3", "ICAM1", "CXCL8", "CCL2", "TRAF1", "NFKB1", "RELB", "SOD2", "CXCL1", "CXCL2", "VCAM1"]
DEATH = ["CFLAR", "CASP8", "RIPK1", "RIPK3", "MLKL", "FADD", "TNFRSF1A", "XIAP", "BIRC2", "BIRC3"]
OUT = "tnf_response_kras"; os.makedirs(OUT, exist_ok=True)
B = "https://www.cbioportal.org/api"


def api(path, body=None):
    for i in range(5):
        try:
            req = urllib.request.Request(B + path, data=json.dumps(body).encode() if body else None,
                                         headers={"Content-Type": "application/json", "Accept": "application/json"})
            return json.load(urllib.request.urlopen(req, timeout=180))
        except urllib.error.HTTPError as e:
            if e.code < 500 or i == 4:
                raise
            time.sleep(3 * (i + 1))


def kras_mut(study):
    S = f"{study}_tcga_pan_can_atlas_2018"
    seq = set(api(f"/sample-lists/{S}_sequenced")["sampleIds"])
    mu = api(f"/molecular-profiles/{S}_mutations/mutations/fetch", {"sampleListId": f"{S}_all", "entrezGeneIds": [3845]})
    mut = {x["sampleId"] for x in mu if not x.get("mutationType", "").startswith(("Silent", "Intron", "3'", "5'", "RNA"))}
    return {s: float(s in mut) for s in seq}


def ols_interaction(y, x, k, cov):
    """y ~ x + k + x*k + cov on ranks; returns beta, p for x*k and slopes of x in k=0 and k=1."""
    d = pd.concat([y, x, k] + cov, axis=1).dropna()
    r = d.rank()
    names = ["y", "x", "k"] + [f"c{i}" for i in range(len(cov))]
    r.columns = names
    r["xk"] = r["x"] * r["k"]
    X = np.column_stack([np.ones(len(r)), r["x"], r["k"], r["xk"]] + [r[n] for n in names[3:]])
    b, *_ = np.linalg.lstsq(X, r["y"].values, rcond=None)
    res = r["y"].values - X @ b
    dof = len(r) - X.shape[1]
    cov_b = np.linalg.inv(X.T @ X) * (res @ res / dof)
    t = b[3] / np.sqrt(cov_b[3, 3])
    return b[3], 2 * stats.t.sf(abs(t), dof)


rows, drows, srows = [], [], []
for c, study in COHORTS.items():
    genes = ["TNF"] + TNF_TARGETS + DEATH
    df = m.add_scores(m.load_tcga(c, genes, "xena_cache"))
    km = kras_mut(study)
    df["KRAS_mut"] = pd.Series({k: v for k, v in km.items() if k in df.index})
    df = df[df["KRAS_mut"].notna()].copy()
    tg = [g for g in TNF_TARGETS if g in df]
    df["TNF_SIG"] = pd.concat([m.zscore(df[g]) for g in tg], axis=1).mean(axis=1)
    msi = ["MSI_H"] if df["MSI_H"].notna().mean() >= 0.7 and df["MSI_H"].nunique() > 1 else []
    adj = ["LEUKOCYTE", "CAF"]
    print(f"[{c}] n={len(df)}, KRAS-mut {int(df.KRAS_mut.sum())}, TNF-target genes used: {len(tg)}")
    for lab, sub in [("KRAS-WT", df[df.KRAS_mut == 0]), ("KRAS-mut", df[df.KRAS_mut == 1])]:
        raw = stats.spearmanr(sub["TNF"], sub["TNF_SIG"], nan_policy="omit")
        pc = m.partial_spearman(sub, "TNF", "TNF_SIG", adj)
        rows.append(dict(cohort=c, group=lab, n=len(sub), rho_raw=raw[0], rho_partial=pc[0], p_partial=pc[1]))
    b, p = ols_interaction(df["TNF_SIG"], df["TNF"], df["KRAS_mut"], [df["LEUKOCYTE"], df["CAF"]])
    rows.append(dict(cohort=c, group="difference (KRAS-mut vs WT slope)", n=len(df), rho_raw=np.nan, rho_partial=b, p_partial=p))
    for g in [x for x in DEATH if x in df]:
        a, w = df.loc[df.KRAS_mut == 1, g], df.loc[df.KRAS_mut == 0, g]
        drows.append(dict(cohort=c, gene=g, diff_mut_minus_wt=a.median() - w.median(), p=stats.mannwhitneyu(a, w).pvalue))
    cov = m.usable_covariates(df, ["stage", "age"] + msi)
    for ep in ["OS", "PFI"]:
        T = f"{ep}.time"
        for lab, sub in [("KRAS-WT", df[df.KRAS_mut == 0]), ("KRAS-mut", df[df.KRAS_mut == 1])]:
            d = sub[[T, ep, "TNF"] + cov].dropna(); d = d[d[T] > 0].copy()
            if d[ep].sum() < 10: continue
            d["X"] = m.zscore(d["TNF"])
            s = CoxPHFitter(penalizer=0.01).fit(d[[T, ep, "X"] + cov], T, ep).summary.loc["X"]
            srows.append(dict(cohort=c, endpoint=ep, group=lab, n=len(d), events=int(d[ep].sum()),
                              HR=s["exp(coef)"], logHR=s["coef"], se=s["se(coef)"], p=s["p"]))
rd, dd, sd = pd.DataFrame(rows), pd.DataFrame(drows), pd.DataFrame(srows)
rd.to_csv(f"{OUT}/tnf_signaling.csv", index=False, float_format="%.4g")
dd.to_csv(f"{OUT}/death_machinery.csv", index=False, float_format="%.4g")
sd.to_csv(f"{OUT}/tnf_survival_by_kras.csv", index=False, float_format="%.4g")

print("\nTNF expression vs TNF-target signature (Spearman; partial = adjusted for leukocyte + CAF)")
print(rd.round(3).to_string(index=False))
print("\nDeath machinery, KRAS-mutant minus wild-type (median log2, Mann-Whitney p)")
for g in DEATH:
    x = dd[dd.gene == g]
    if len(x): print(f"  {g:9s} " + "  ".join(f"{r.cohort}: {r.diff_mut_minus_wt:+.2f} (p={r.p:.2f})" for _, r in x.iterrows()))
print("\nTNF expression and survival (HR per SD, p) by KRAS group")
for ep in ["OS", "PFI"]:
    for lab in ["KRAS-WT", "KRAS-mut"]:
        x = sd[(sd.endpoint == ep) & (sd.group == lab)].set_index("cohort")
        mm = m.meta(x["logHR"], x["se"])
        mt = f"META {mm['HR_random']:.2f} ({mm['CI_low']:.2f}-{mm['CI_high']:.2f}) p={mm['p_random']:.3f}" if "HR_random" in mm else ""
        print(f"  {ep} {lab:9s} " + "  ".join(f"{k}: {x.loc[k,'HR']:.2f} ({x.loc[k,'p']:.2f}; ev {x.loc[k,'events']})" for k in x.index) + f" | {mt}")
