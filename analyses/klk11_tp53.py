import json, urllib.request, statistics as st, math
B="https://www.cbioportal.org/api"
S="coadread_tcga_pan_can_atlas_2018"
def post(p,b):
    r=urllib.request.Request(B+p,data=json.dumps(b).encode(),headers={"Content-Type":"application/json","Accept":"application/json"})
    return json.load(urllib.request.urlopen(r,timeout=120))
def get(p): return json.load(urllib.request.urlopen(urllib.request.Request(B+p,headers={"Accept":"application/json"}),timeout=120))
genes={"KLK11":11012,"TP53":7157,"KRAS":3845,"NRAS":4893,"BRAF":673}
expr=post(f"/molecular-profiles/{S}_rna_seq_v2_mrna/molecular-data/fetch",{"sampleListId":f"{S}_all","entrezGeneIds":[11012]})
ex={d["sampleId"]:math.log2(d["value"]+1) for d in expr if d.get("value") is not None}
muts=post(f"/molecular-profiles/{S}_mutations/mutations/fetch",{"sampleListId":f"{S}_all","entrezGeneIds":[7157,3845,4893,673]})
mut={}
for m in muts: mut.setdefault(m["sampleId"],set()).add({v:k for k,v in genes.items()}[m["entrezGeneId"]])
seq=set(get(f"/sample-lists/{S}_sequenced")["sampleIds"])
samples=[s for s in ex if s in seq]
def grp(f): return [ex[s] for s in samples if f(mut.get(s,set()))]
def ranksum(a,b):
    allv=sorted([(v,0) for v in a]+[(v,1) for v in b]); r=[0]*len(allv); i=0
    while i<len(allv):
        j=i
        while j+1<len(allv) and allv[j+1][0]==allv[i][0]: j+=1
        for k in range(i,j+1): r[k]=(i+j)/2+1
        i=j+1
    R=sum(r[k] for k in range(len(allv)) if allv[k][1]==0); n1,n2=len(a),len(b)
    U=R-n1*(n1+1)/2; mu=n1*n2/2; sd=math.sqrt(n1*n2*(n1+n2+1)/12); z=(U-mu)/sd
    return z, math.erfc(abs(z)/math.sqrt(2))
def rep(name,a,b,la,lb):
    z,p=ranksum(a,b)
    print(f"{name}: {la} n={len(a)} median={st.median(a):.2f} | {lb} n={len(b)} median={st.median(b):.2f} | diff={st.median(a)-st.median(b):+.2f} log2 | p={p:.2e}")
wt=lambda m:not({"KRAS","NRAS","BRAF"}&m)
rep("All CRC", grp(lambda m:"TP53" in m), grp(lambda m:"TP53" not in m),"TP53mut","TP53wt")
rep("RAS/BRAF-WT only", grp(lambda m:"TP53" in m and wt(m)), grp(lambda m:"TP53" not in m and wt(m)),"TP53mut","TP53wt")
rep("TP53mut: WT vs RAS/BRAF-mut", grp(lambda m:"TP53" in m and wt(m)), grp(lambda m:"TP53" in m and not wt(m)),"RAS/BRAF-WT","RAS/BRAF-mut")
# normal tissue not in pan-can; overall distribution
allv=list(ex.values()); print(f"KLK11 log2(RSEM+1) overall: median={st.median(allv):.2f}, IQR={sorted(allv)[len(allv)//4]:.2f}-{sorted(allv)[3*len(allv)//4]:.2f}, n={len(allv)}")

# Confounder check: restrict to chromosomally unstable (CIN, i.e. MSS) tumors via TCGA molecular subtype
sub={}
for d in get(f"/studies/{S}/clinical-data?clinicalDataType=PATIENT&attributeId=SUBTYPE&pageSize=10000"):
    sub[d["patientId"]]=d["value"]
pid=lambda s:s[:12]
from collections import Counter
print("Subtypes:",Counter(sub.values()))
cin=[s for s in samples if sub.get(pid(s),"").endswith("CIN")]
samples_bak=samples; samples=cin
rep("CIN/MSS only", grp(lambda m:"TP53" in m), grp(lambda m:"TP53" not in m),"TP53mut","TP53wt")
rep("CIN/MSS, RAS/BRAF-WT", grp(lambda m:"TP53" in m and wt(m)), grp(lambda m:"TP53" not in m and wt(m)),"TP53mut","TP53wt")
samples=samples_bak
for t in sorted(set(sub.values())):
    v=[ex[s] for s in samples if sub.get(pid(s))==t]
    if v: print(f"  {t}: n={len(v)} median KLK11={st.median(v):.2f}")
