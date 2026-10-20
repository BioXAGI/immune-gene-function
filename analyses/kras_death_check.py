import json, urllib.request, warnings
import numpy as np, pandas as pd
from scipy import stats
import immune_gene_function as m
warnings.filterwarnings("ignore")
G=["CFLAR","BIRC2","BIRC3","XIAP","CASP8","RIPK1","RIPK3","MLKL","FADD","TNFRSF1A"]
S="coadread_tcga_pan_can_atlas_2018"; B="https://www.cbioportal.org/api"
def api(p,b=None):
    r=urllib.request.Request(B+p,data=json.dumps(b).encode() if b else None,headers={"Content-Type":"application/json","Accept":"application/json"})
    return json.load(urllib.request.urlopen(r,timeout=180))
seq=set(api(f"/sample-lists/{S}_sequenced")["sampleIds"])
mu={x["sampleId"] for x in api(f"/molecular-profiles/{S}_mutations/mutations/fetch",{"sampleListId":f"{S}_all","entrezGeneIds":[3845]}) if not x.get("mutationType","").startswith(("Silent","Intron","3'","5'","RNA"))}
df=m.add_scores(m.load_tcga("COADREAD",G,"xena_cache"))
df["K"]=pd.Series({s:float(s in mu) for s in seq}); df=df[df.K.notna()].copy()
def ols(y,covs):
    d=df[[y,"K"]+covs].dropna(); X=np.column_stack([np.ones(len(d)),d.K]+[d[c] for c in covs]); yv=d[y].values
    b=np.linalg.lstsq(X,yv,rcond=None)[0]; res=yv-X@b; dof=len(d)-X.shape[1]
    se=np.sqrt(np.linalg.inv(X.T@X)[1,1]*(res@res/dof)); return b[1],2*stats.t.sf(abs(b[1]/se),dof),len(d)
sd=lambda y: df[y].std()
rows=[]
for g in G:
    r={"gene":g}
    a=df[df.K==1][g];w=df[df.K==0][g]; r["raw"]=a.median()-w.median(); r["p_raw"]=stats.mannwhitneyu(a,w).pvalue
    for name,cv in [("MSI",["MSI_H"]),("MSI+leuk+CAF",["MSI_H","LEUKOCYTE","CAF"])]:
        b,p,n=ols(g,cv); r[f"b_{name}"]=b; r[f"p_{name}"]=p
    mss=df[df.MSI_H==0]; a=mss[mss.K==1][g];w=mss[mss.K==0][g]; r["MSS_diff"]=a.median()-w.median(); r["p_MSS"]=stats.mannwhitneyu(a,w).pvalue
    r["sd"]=sd(g); rows.append(r)
t=pd.DataFrame(rows)
p=t.p_raw.values; o=np.argsort(p); q=np.empty(len(p)); q[o]=np.minimum.accumulate((p[o]*len(p)/(np.arange(len(p))+1))[::-1])[::-1]; t["q_raw"]=np.minimum(q,1)
print(f"n={len(df)}, KRAS-mut {int(df.K.sum())}, MSI-H {int(df.MSI_H.sum())}")
print(t.round(3).to_string(index=False))
