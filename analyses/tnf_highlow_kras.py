"""HR of TNF-high vs TNF-low (median split within cohort and KRAS group) for OS and PFI.
Group difference: Cox on the whole cohort with TNF_high + KRAS_mut + TNF_high*KRAS_mut.
Also reports the continuous per-SD HR. Melanoma has ~2% KRAS mutants, so it is shown as KRAS-WT-like only."""
import json, os, time, urllib.error, urllib.request, warnings
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
import immune_gene_function as m
warnings.filterwarnings("ignore")
COHORTS = {"COADREAD": "coadread", "LUAD": "luad", "PAAD": "paad", "SKCM": "skcm"}
OUT = "tnf_highlow_kras"; os.makedirs(OUT, exist_ok=True)
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


def cox(d, T, ep, terms):
    return CoxPHFitter(penalizer=0.01).fit(d[[T, ep] + terms], T, ep).summary


rows = []
for c, study in COHORTS.items():
    df = m.load_tcga(c, ["TNF"], "xena_cache")
    df["KRAS_mut"] = pd.Series(kras_mut(study))
    df = df[df["KRAS_mut"].notna()].copy()
    msi = ["MSI_H"] if df["MSI_H"].notna().mean() >= 0.7 and df["MSI_H"].nunique() > 1 else []
    cov = m.usable_covariates(df, ["stage", "age"] + msi)
    for ep in ["OS", "PFI"]:
        T = f"{ep}.time"
        d = df[[T, ep, "TNF", "KRAS_mut"] + cov].dropna(); d = d[d[T] > 0].copy()
        for lab, k in [("KRAS-WT", 0.0), ("KRAS-mut", 1.0)]:
            sub = d[d.KRAS_mut == k].copy()
            if sub[ep].sum() < 10:
                continue
            sub["TNF_high"] = (sub["TNF"] >= sub["TNF"].median()).astype(float)
            s = cox(sub, T, ep, ["TNF_high"] + cov).loc["TNF_high"]
            rows.append(dict(cohort=c, endpoint=ep, group=lab, n=len(sub), events=int(sub[ep].sum()), term="TNF high vs low",
                             HR=s["exp(coef)"], logHR=s["coef"], se=s["se(coef)"], p=s["p"]))
        if d.loc[d.KRAS_mut == 1, ep].sum() >= 10:
            d["TNF_high"] = d.groupby("KRAS_mut")["TNF"].transform(lambda x: (x >= x.median()).astype(float))
            d["HxK"] = d["TNF_high"] * d["KRAS_mut"]
            s = cox(d, T, ep, ["TNF_high", "KRAS_mut", "HxK"] + cov).loc["HxK"]
            rows.append(dict(cohort=c, endpoint=ep, group="difference", n=len(d), events=int(d[ep].sum()),
                             term="interaction (KRAS-mut vs WT)", HR=s["exp(coef)"], logHR=s["coef"], se=s["se(coef)"], p=s["p"]))
res = pd.DataFrame(rows); res.to_csv(f"{OUT}/tnf_highlow.csv", index=False, float_format="%.4g")
meta = {}
for ep in ["OS", "PFI"]:
    print(f"\n===== {ep}: TNF-high vs TNF-low, HR (p; events) =====")
    for grp in ["KRAS-WT", "KRAS-mut", "difference"]:
        x = res[(res.endpoint == ep) & (res.group == grp)].set_index("cohort")
        if x.empty: continue
        mm = m.meta(x["logHR"], x["se"])
        mt = f"META {mm['HR_random']:.2f} ({mm['CI_low']:.2f}-{mm['CI_high']:.2f}) p={mm['p_random']:.3f} I2={mm['I2']:.0%}" if "HR_random" in mm else ""
        print(f"  {grp:10s} " + "  ".join(f"{k}: {x.loc[k,'HR']:.2f} ({x.loc[k,'p']:.2f}; {x.loc[k,'events']})" for k in x.index) + f" | {mt}")
    sub = res[(res.endpoint == ep) & (res.group.isin(["KRAS-WT", "KRAS-mut"])) & (res.cohort != "SKCM")]
    for grp in ["KRAS-WT", "KRAS-mut"]:
        x = sub[sub.group == grp]; mm = m.meta(x["logHR"], x["se"])
        print(f"  [excl. melanoma] {grp:8s} META {mm['HR_random']:.2f} ({mm['CI_low']:.2f}-{mm['CI_high']:.2f}) p={mm['p_random']:.3f}")
