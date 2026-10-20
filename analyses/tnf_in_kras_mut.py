import json, time, urllib.error, urllib.request, warnings
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
import immune_gene_function as m
warnings.filterwarnings("ignore")
B="https://www.cbioportal.org/api"
def api(p,b=None):
    for i in range(5):
        try:
            r=urllib.request.Request(B+p,data=json.dumps(b).encode() if b else None,headers={"Content-Type":"application/json","Accept":"application/json"})
            return json.load(urllib.request.urlopen(r,timeout=180))
        except urllib.error.HTTPError as e:
            if e.code<500 or i==4: raise
            time.sleep(3*(i+1))
def kras(study):
    S=f"{study}_tcga_pan_can_atlas_2018"; seq=set(api(f"/sample-lists/{S}_sequenced")["sampleIds"])
    mu={x["sampleId"] for x in api(f"/molecular-profiles/{S}_mutations/mutations/fetch",{"sampleListId":f"{S}_all","entrezGeneIds":[3845]}) if not x.get("mutationType","").startswith(("Silent","Intron","3'","5'","RNA"))}
    return {s:float(s in mu) for s in seq}
EPI=["EPCAM","KRT8","KRT18","KRT19","CDH1"]
MODELS={"stage+age(+MSI)":[], "+composition":["LEUKOCYTE","CAF","EPITHELIAL"], "+composition+CTL":["LEUKOCYTE","CAF","EPITHELIAL","CTL"]}
rows=[]
for c,s in {"COADREAD":"coadread","LUAD":"luad","PAAD":"paad"}.items():
    df=m.add_scores(m.load_tcga(c,["TNF"]+EPI,"xena_cache"))
    df["EPITHELIAL"]=pd.concat([m.zscore(df[g]) for g in EPI if g in df],axis=1).mean(axis=1)
    df["K"]=pd.Series(kras(s)); df=df[df.K==1].copy()
    msi=["MSI_H"] if df["MSI_H"].notna().mean()>=0.7 and df["MSI_H"].nunique()>1 else []
    clin=m.usable_covariates(df,["stage","age"]+msi)
    for ep in ["OS","PFI"]:
        T=f"{ep}.time"
        for name,extra in MODELS.items():
            d=df[[T,ep,"TNF"]+clin+extra].dropna(); d=d[d[T]>0].copy()
            d["X"]=m.zscore(d["TNF"])
            for e in extra: d[e]=m.zscore(d[e])
            su=CoxPHFitter(penalizer=0.01).fit(d[[T,ep,"X"]+clin+extra],T,ep).summary.loc["X"]
            hi=(d["TNF"]>=d["TNF"].median()).astype(float); d2=d.assign(H=hi)
            sh=CoxPHFitter(penalizer=0.01).fit(d2[[T,ep,"H"]+clin+extra],T,ep).summary.loc["H"]
            rows.append(dict(cohort=c,endpoint=ep,model=name,n=len(d),events=int(d[ep].sum()),HR_perSD=su["exp(coef)"],logHR=su["coef"],se=su["se(coef)"],p=su["p"],HR_high_vs_low=sh["exp(coef)"],p_hl=sh["p"]))
r=pd.DataFrame(rows); r.to_csv("tnf_in_kras_mut.csv",index=False,float_format="%.4g")
for ep in ["OS","PFI"]:
    print(f"\n== {ep}, KRAS-mutant only: TNF HR per SD (p) | high vs low (p) ==")
    for name in MODELS:
        x=r[(r.endpoint==ep)&(r.model==name)].set_index("cohort"); mm=m.meta(x.logHR,x.se)
        print(f"  {name:20s} "+"  ".join(f"{k}: {x.loc[k,'HR_perSD']:.2f} ({x.loc[k,'p']:.2f}) [{x.loc[k,'HR_high_vs_low']:.2f}; ev {x.loc[k,'events']}]" for k in x.index)+f" | META {mm['HR_random']:.2f} ({mm['CI_low']:.2f}-{mm['CI_high']:.2f}) p={mm['p_random']:.3f} I2={mm['I2']:.0%}")
