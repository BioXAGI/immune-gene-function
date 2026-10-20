"""TNF-response and cell-death scores in TCGA with composition adjustment.
Scores (mean z of genes): TNF_RESP, DEATH_RESIST, DEATH_EXEC, NET_SENS (= EXEC - RESIST).
Composition covariates: LEUKOCYTE (PTPRC), CAF, ENDOTHELIAL, EPITHELIAL (EPCAM KRT8 KRT18 KRT19 CDH1; purity proxy), MSI-H.
Analyses
  1. score ~ KRAS_mut (+ composition), standardized beta, unadjusted and adjusted; meta across cohorts
  2. TNF expression -> TNF_RESP slope, adjusted, by KRAS group and difference
  3. survival: hazard ~ TNF + SCORE + TNF*SCORE + stage + age [+ MSI-H] + composition; meta across cohorts
"""
import json, os, time, urllib.error, urllib.request, warnings
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
from scipy import stats
import immune_gene_function as m
warnings.filterwarnings("ignore")

SETS = {
    "TNF_RESP": ["NFKBIA", "TNFAIP3", "ICAM1", "CXCL8", "CCL2", "TRAF1", "NFKB1", "RELB", "SOD2", "CXCL1", "CXCL2", "VCAM1"],
    "DEATH_RESIST": ["CFLAR", "BIRC2", "BIRC3", "XIAP", "BCL2L1", "MCL1", "BCL2"],
    "DEATH_EXEC": ["RIPK1", "RIPK3", "MLKL", "CASP8", "FADD", "TNFRSF1A"],
    "EPITHELIAL": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"],
}
SCORES = ["TNF_RESP", "DEATH_RESIST", "DEATH_EXEC", "NET_SENS"]
COHORTS = {"COADREAD": "coadread", "LUAD": "luad", "PAAD": "paad", "SKCM": "skcm"}
COMP = ["LEUKOCYTE", "CAF", "ENDOTHELIAL", "EPITHELIAL"]
OUT = "tnf_death_scores"; os.makedirs(OUT, exist_ok=True)
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


def ols(d, y, x, covs):
    d = d[[y, x] + covs].dropna()
    Y = m.zscore(d[y]).values
    X = np.column_stack([np.ones(len(d)), d[x].astype(float)] + [m.zscore(d[c]) for c in covs])
    b = np.linalg.lstsq(X, Y, rcond=None)[0]
    res = Y - X @ b; dof = len(d) - X.shape[1]
    se = np.sqrt(np.linalg.inv(X.T @ X)[1, 1] * (res @ res / dof))
    return b[1], se, 2 * stats.t.sf(abs(b[1] / se), dof), len(d)


kr, tr, sr = [], [], []
for c, study in COHORTS.items():
    allg = ["TNF"] + sum(SETS.values(), [])
    df = m.add_scores(m.load_tcga(c, allg, "xena_cache"))
    for name, gl in SETS.items():
        have = [g for g in gl if g in df and df[g].std() > 0]
        df[name] = pd.concat([m.zscore(df[g]) for g in have], axis=1).mean(axis=1)
        if len(have) < len(gl): print(f"[{c}] {name}: missing {[g for g in gl if g not in have]}")
    df["NET_SENS"] = m.zscore(df["DEATH_EXEC"]) - m.zscore(df["DEATH_RESIST"])
    df["KRAS_mut"] = pd.Series(kras_mut(study))
    msi = ["MSI_H"] if df["MSI_H"].notna().mean() >= 0.7 and df["MSI_H"].nunique() > 1 else []
    comp = COMP + msi
    d = df[df["KRAS_mut"].notna()].copy()
    print(f"[{c}] n={len(df)}; with KRAS calls {len(d)}, mutant {int(d.KRAS_mut.sum())}")
    if d.KRAS_mut.sum() >= 20:
        for s in SCORES:
            for lab, cv in [("unadjusted", []), ("composition-adjusted", comp)]:
                b, se, p, n = ols(d, s, "KRAS_mut", cv)
                kr.append(dict(cohort=c, score=s, model=lab, beta_SD=b, se=se, p=p, n=n))
        for lab, sub in [("KRAS-WT", d[d.KRAS_mut == 0]), ("KRAS-mut", d[d.KRAS_mut == 1])]:
            b, se, p, n = ols(sub, "TNF_RESP", "TNF", [x for x in comp])
            tr.append(dict(cohort=c, group=lab, beta_SD=b, se=se, p=p, n=n))
        d["TNFxK"] = m.zscore(d["TNF"]) * d["KRAS_mut"]
        dd = d[["TNF_RESP", "TNF", "KRAS_mut", "TNFxK"] + comp].dropna(); dd = dd.assign(TNF=m.zscore(dd["TNF"]))
        dd["TNFxK"] = dd["TNF"] * dd["KRAS_mut"]
        Y = m.zscore(dd["TNF_RESP"]).values; X = np.column_stack([np.ones(len(dd)), dd.TNF, dd.KRAS_mut, dd.TNFxK] + [m.zscore(dd[x]) for x in comp])
        bb = np.linalg.lstsq(X, Y, rcond=None)[0]; rs = Y - X @ bb; dof = len(dd) - X.shape[1]
        se = np.sqrt(np.linalg.inv(X.T @ X)[3, 3] * (rs @ rs / dof))
        tr.append(dict(cohort=c, group="difference (mut - WT)", beta_SD=bb[3], se=se, p=2 * stats.t.sf(abs(bb[3] / se), dof), n=len(dd)))
    clin = m.usable_covariates(df, ["stage", "age"] + msi)
    for ep in ["OS", "PFI"]:
        T = f"{ep}.time"
        for s in SCORES:
            z = df[[T, ep, "TNF", s] + COMP + clin].dropna(); z = z[z[T] > 0].copy()
            z["TNF_z"], z["S"] = m.zscore(z["TNF"]), m.zscore(z[s]); z["TxS"] = z["TNF_z"] * z["S"]
            terms = ["TNF_z", "S", "TxS"] + clin + ["LEUKOCYTE", "CAF", "EPITHELIAL"]
            for cc in ["LEUKOCYTE", "CAF", "EPITHELIAL"]: z[cc] = m.zscore(z[cc])
            su = CoxPHFitter(penalizer=0.01).fit(z[[T, ep] + terms], T, ep).summary
            for t in ["TNF_z", "S", "TxS"]:
                sr.append(dict(cohort=c, endpoint=ep, score=s, term=t, events=int(z[ep].sum()), HR=su.loc[t, "exp(coef)"],
                               logHR=su.loc[t, "coef"], se=su.loc[t, "se(coef)"], p=su.loc[t, "p"]))
kd, td, sd = pd.DataFrame(kr), pd.DataFrame(tr), pd.DataFrame(sr)
kd.to_csv(f"{OUT}/kras_scores.csv", index=False, float_format="%.4g"); td.to_csv(f"{OUT}/tnf_to_response.csv", index=False, float_format="%.4g")
sd.to_csv(f"{OUT}/survival_models.csv", index=False, float_format="%.4g")


def pool(x, col="beta_SD"):
    mm = m.meta(x[col].values if col == "logHR" else np.exp(x[col].values), x["se"].values) if False else None
    return mm


print("\n1) KRAS-mutant minus wild-type, in SD of the score (beta; p)")
for s in SCORES:
    for lab in ["unadjusted", "composition-adjusted"]:
        x = kd[(kd.score == s) & (kd.model == lab)].set_index("cohort")
        mm = m.meta(x["beta_SD"].values, x["se"].values)  # fixed/random on beta treated as log-scale
        lg = np.log(mm["HR_random"]) if "HR_random" in mm else np.nan
        print(f"  {s:13s} {lab:21s} " + "  ".join(f"{k}: {x.loc[k,'beta_SD']:+.2f} ({x.loc[k,'p']:.3f})" for k in x.index) +
              (f" | META {lg:+.2f} p={mm['p_random']:.3f}" if "p_random" in mm else ""))
print("\n2) TNF expression -> TNF-response score slope (adjusted for composition)")
print(td.round(3).to_string(index=False))
print("\n3) Survival: HR (p) per SD; TxS = TNF x score interaction; adjusted for stage, age, leukocyte, CAF, epithelial")
for ep in ["OS", "PFI"]:
    for s in SCORES:
        for t in ["TNF_z", "S", "TxS"]:
            x = sd[(sd.endpoint == ep) & (sd.score == s) & (sd.term == t)].set_index("cohort")
            mm = m.meta(x["logHR"], x["se"])
            mt = f"META {mm['HR_random']:.2f} ({mm['CI_low']:.2f}-{mm['CI_high']:.2f}) p={mm['p_random']:.3f}" if "HR_random" in mm else ""
            print(f"  {ep} {s:13s} {t:5s} " + "  ".join(f"{k}: {x.loc[k,'HR']:.2f} ({x.loc[k,'p']:.2f})" for k in x.index) + f" | {mt}")
