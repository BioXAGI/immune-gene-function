"""KRAS-mutation survival models in TCGA (mutation calls from cBioPortal PanCancer Atlas).
Per cohort and endpoint (Cox, stage + age [+ MSI-H] adjusted):
  A  KRAS mutant vs wild type
  B  allele groups (G12D, G12V, G12C, G13x, other) vs wild type (groups with >= MIN_ALLELE tumors)
  C  gene x KRAS-mutation interaction for BIRC2, BIRC3, TNF (gene z-scored)
  D  BIRC2 x TNF interaction within KRAS-mutant and within KRAS-wild-type tumors
Also compares BIRC2/BIRC3/TNF expression by KRAS status. Random-effects meta-analysis across cohorts."""
import json, os, time, urllib.error, urllib.request, warnings
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
from scipy import stats
import immune_gene_function as m
warnings.filterwarnings("ignore")

COHORTS = {"COADREAD": "coadread", "LUAD": "luad", "PAAD": "paad", "SKCM": "skcm"}
GENES = ["BIRC2", "BIRC3", "TNF"]
MIN_ALLELE, OUT = 10, "kras_model"
os.makedirs(OUT, exist_ok=True)
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


def kras_status(study):
    S = f"{study}_tcga_pan_can_atlas_2018"
    seq = set(api(f"/sample-lists/{S}_sequenced")["sampleIds"])
    mu = api(f"/molecular-profiles/{S}_mutations/mutations/fetch",
             {"sampleListId": f"{S}_all", "entrezGeneIds": [3845]})
    alle = {}
    for x in mu:
        if x.get("mutationType", "").startswith(("Silent", "Intron", "3'", "5'", "RNA")):
            continue
        alle.setdefault(x["sampleId"], []).append(x.get("proteinChange", "?"))
    def grp(ch):
        s = ",".join(ch)
        for k in ["G12D", "G12V", "G12C", "G12A", "G12R", "G12S"]:
            if k in s: return k
        if "G13" in s: return "G13x"
        if "Q61" in s: return "Q61x"
        return "other"
    return {sid: ("WT" if sid not in alle else grp(alle[sid])) for sid in seq}


def cox(d, T, ep, terms):
    return CoxPHFitter(penalizer=0.01).fit(d[[T, ep] + terms], T, ep).summary


rows, expr_rows = [], []
for c, study in COHORTS.items():
    df = m.load_tcga(c, GENES, "xena_cache")
    st = kras_status(study)
    df["KRAS_group"] = pd.Series({k: v for k, v in st.items() if k in df.index})
    df = df[df["KRAS_group"].notna()].copy()
    df["KRAS_mut"] = (df["KRAS_group"] != "WT").astype(float)
    cnt = df["KRAS_group"].value_counts().to_dict()
    msi = ["MSI_H"] if df["MSI_H"].notna().mean() >= 0.7 and df["MSI_H"].nunique() > 1 else []
    clin = m.usable_covariates(df, ["stage", "age"] + msi)
    print(f"[{c}] tumors with mutation calls: {len(df)}; KRAS mutant {int(df.KRAS_mut.sum())} "
          f"({100*df.KRAS_mut.mean():.0f}%); groups {cnt}")
    for g in GENES:
        a, b = df.loc[df.KRAS_mut == 1, g].dropna(), df.loc[df.KRAS_mut == 0, g].dropna()
        if len(a) > 3 and len(b) > 3:
            expr_rows.append(dict(cohort=c, gene=g, n_mut=len(a), n_wt=len(b), diff=a.median() - b.median(),
                                  p=stats.mannwhitneyu(a, b).pvalue))
    for ep in ["OS", "PFI"]:
        T = f"{ep}.time"
        d0 = df[df[T].notna() & df[ep].notna() & (df[T] > 0)]
        d0 = d0.dropna(subset=clin + GENES).copy()
        if d0["KRAS_mut"].nunique() < 2 or d0.loc[d0.KRAS_mut == 1, ep].sum() < 5:
            continue
        for g in GENES:
            d0[g] = m.zscore(d0[g])
        ev = int(d0[ep].sum())
        def put(model, term, s):
            rows.append(dict(cohort=c, endpoint=ep, model=model, term=term, events=ev, n=len(d0),
                             HR=s["exp(coef)"], logHR=s["coef"], se=s["se(coef)"], p=s["p"]))
        s = cox(d0, T, ep, ["KRAS_mut"] + clin); put("A KRAS mut vs WT", "KRAS_mut", s.loc["KRAS_mut"])
        d1 = d0.copy(); terms = []
        for k, n in d1["KRAS_group"].value_counts().items():
            if k != "WT" and n >= MIN_ALLELE:
                d1[f"k_{k}"] = (d1["KRAS_group"] == k).astype(float); terms.append(f"k_{k}")
        if terms:
            s = cox(d1, T, ep, terms + clin)
            for t in terms: put("B allele vs WT", t, s.loc[t])
        for g in GENES:
            d0["XxK"] = d0[g] * d0["KRAS_mut"]
            s = cox(d0, T, ep, [g, "KRAS_mut", "XxK"] + clin)
            put("C gene x KRASmut", f"{g}: effect in KRAS-WT", s.loc[g]); put("C gene x KRASmut", f"{g}: interaction", s.loc["XxK"])
        for lab, sub in [("KRAS-mut", d0[d0.KRAS_mut == 1]), ("KRAS-WT", d0[d0.KRAS_mut == 0])]:
            if sub[ep].sum() < 10: continue
            sub = sub.copy(); sub["BxT"] = sub["BIRC2"] * sub["TNF"]
            s = cox(sub, T, ep, ["BIRC2", "TNF", "BxT"] + clin)
            put("D BIRC2 x TNF", f"{lab}: BIRC2", s.loc["BIRC2"]); put("D BIRC2 x TNF", f"{lab}: TNF", s.loc["TNF"])
            put("D BIRC2 x TNF", f"{lab}: BIRC2xTNF", s.loc["BxT"])
res = pd.DataFrame(rows); res.to_csv(f"{OUT}/cohort_models.csv", index=False, float_format="%.4g")
pd.DataFrame(expr_rows).to_csv(f"{OUT}/expression_by_kras.csv", index=False, float_format="%.4g")
meta = pd.DataFrame([dict(endpoint=ep, model=mo, term=t, **m.meta(x["logHR"], x["se"]))
                     for (ep, mo, t), x in res.groupby(["endpoint", "model", "term"])])
meta.to_csv(f"{OUT}/meta_analysis.csv", index=False, float_format="%.4g")

print("\nExpression, KRAS-mutant minus wild-type (median log2; Mann-Whitney p):")
e = pd.DataFrame(expr_rows)
for g in GENES:
    print(f"  {g:6s} " + "  ".join(f"{r.cohort}: {r['diff']:+.2f} (p={r.p:.2f})" for _, r in e[e.gene == g].iterrows()))
for ep in ["OS", "PFI"]:
    print(f"\n===== {ep}: HR (p) by cohort | meta =====")
    for mo in ["A KRAS mut vs WT", "B allele vs WT", "C gene x KRASmut", "D BIRC2 x TNF"]:
        print(f"-- {mo}")
        x = res[(res.endpoint == ep) & (res.model == mo)]
        for t in x.term.unique():
            y = x[x.term == t].set_index("cohort")
            mm = meta[(meta.endpoint == ep) & (meta.model == mo) & (meta.term == t)].iloc[0]
            mt = f"META {mm.HR_random:.2f} ({mm.CI_low:.2f}-{mm.CI_high:.2f}) p={mm.p_random:.3f}" if pd.notna(mm.get("HR_random")) else f"(k={mm.k})"
            print(f"  {t:24s} " + "  ".join(f"{k}: {y.loc[k,'HR']:.2f} ({y.loc[k,'p']:.2f})" for k in y.index) + f" | {mt}")
