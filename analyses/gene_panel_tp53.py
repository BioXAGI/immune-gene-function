import json, urllib.request, urllib.error, statistics as st, math, csv
from collections import Counter
B="https://www.cbioportal.org/api"; S="coadread_tcga_pan_can_atlas_2018"
import time
def _open(req):
    for i in range(5):
        try: return json.load(urllib.request.urlopen(req,timeout=180))
        except urllib.error.HTTPError as e:
            if e.code<500 or i==4: raise
            time.sleep(3*(i+1))
def post(p,b):
    return _open(urllib.request.Request(B+p,data=json.dumps(b).encode(),headers={"Content-Type":"application/json","Accept":"application/json"}))
def get(p): return _open(urllib.request.Request(B+p,headers={"Accept":"application/json"}))

PANEL=["PLAGL2","REEP1","ASXL1","POFUT1","SMTNL2","CAMKV","TRPM6","WASF3","SPATA25"]
IMMUNE=["CD8A","GZMA","PRF1"]
DRIVERS=["TP53","KRAS","NRAS","BRAF"]
ent={}
for x in PANEL+IMMUNE+DRIVERS:
    try: ent[x]=get(f"/genes/{x}")["entrezGeneId"]
    except urllib.error.HTTPError: pass
sym={v:k for k,v in ent.items()}
missing=[x for x in PANEL if x not in ent]
if missing: print("Not found:",missing)

expr=post(f"/molecular-profiles/{S}_rna_seq_v2_mrna/molecular-data/fetch",
          {"sampleListId":f"{S}_all","entrezGeneIds":[ent[x] for x in PANEL+IMMUNE if x in ent]})
E={}
for d in expr:
    if d.get("value") is not None: E.setdefault(sym[d["entrezGeneId"]],{})[d["sampleId"]]=math.log2(d["value"]+1)
muts=post(f"/molecular-profiles/{S}_mutations/mutations/fetch",
          {"sampleListId":f"{S}_all","entrezGeneIds":[ent[x] for x in PANEL+DRIVERS if x in ent]})
mut={}
for m in muts: mut.setdefault(m["sampleId"],set()).add(sym[m["entrezGeneId"]])
seq=set(get(f"/sample-lists/{S}_sequenced")["sampleIds"])
sub={d["patientId"]:d["value"] for d in get(f"/studies/{S}/clinical-data?clinicalDataType=PATIENT&attributeId=SUBTYPE&pageSize=10000")}
samples=[s for s in E["CD8A"] if s in seq]
cyt={s:(E["CD8A"][s]+E["GZMA"][s]+E["PRF1"][s])/3 for s in samples}
mss=[s for s in samples if sub.get(s[:12],"").endswith("CIN")]

def ranksum(a,b):
    v=sorted([(x,0) for x in a]+[(x,1) for x in b]); r=[0]*len(v); i=0
    while i<len(v):
        j=i
        while j+1<len(v) and v[j+1][0]==v[i][0]: j+=1
        for k in range(i,j+1): r[k]=(i+j)/2+1
        i=j+1
    n1,n2=len(a),len(b); U=sum(r[k] for k in range(len(v)) if v[k][1]==0)-n1*(n1+1)/2
    z=(U-n1*n2/2)/math.sqrt(n1*n2*(n1+n2+1)/12); return math.erfc(abs(z)/math.sqrt(2))
def rank(x):
    o=sorted(range(len(x)),key=lambda i:x[i]); r=[0]*len(x)
    for p,i in enumerate(o): r[i]=p
    return r
def spearman(a,b):
    ra,rb=rank(a),rank(b); ma,mb=st.mean(ra),st.mean(rb)
    num=sum((x-ma)*(y-mb) for x,y in zip(ra,rb)); den=math.sqrt(sum((x-ma)**2 for x in ra)*sum((y-mb)**2 for y in rb))
    return num/den if den else 0
def cmp(gene,ss):
    a=[E[gene][s] for s in ss if s in E[gene] and "TP53" in mut.get(s,set())]
    b=[E[gene][s] for s in ss if s in E[gene] and "TP53" not in mut.get(s,set())]
    return st.median(a)-st.median(b), ranksum(a,b)

rows=[]
n=len(samples)
print(f"TCGA CRC sequenced samples with RNA: {n} (MSS/CIN: {len(mss)})\n")
print(f"{'gene':8s} {'median':>6s} {'mut%':>5s} {'TP53 d(all)':>11s} {'p':>8s} {'TP53 d(MSS)':>11s} {'p':>8s} {'rho_immune':>10s} {'MSI med':>7s} {'CIN med':>7s}")
for gname in PANEL:
    if gname not in E: continue
    vals=[E[gname][s] for s in samples if s in E[gname]]
    d1,p1=cmp(gname,samples); d2,p2=cmp(gname,mss)
    rho=spearman([E[gname][s] for s in samples],[cyt[s] for s in samples])
    msi=st.median([E[gname][s] for s in samples if sub.get(s[:12],"").endswith("MSI")])
    cin=st.median([E[gname][s] for s in mss])
    mp=100*sum(1 for s in samples if gname in mut.get(s,set()))/n
    print(f"{gname:8s} {st.median(vals):6.2f} {mp:5.1f} {d1:+11.2f} {p1:8.1e} {d2:+11.2f} {p2:8.1e} {rho:+10.2f} {msi:7.2f} {cin:7.2f}")
    rows.append([gname,round(st.median(vals),2),round(mp,1),round(d1,2),f"{p1:.2e}",round(d2,2),f"{p2:.2e}",round(rho,2),round(msi,2),round(cin,2)])
with open("output/gene_panel_tp53.csv","w",newline="",encoding="utf-8-sig") as f:
    w=csv.writer(f); w.writerow(["gene","median_log2_expr","mutation_pct","TP53mut_vs_wt_diff_all","p_all","TP53mut_vs_wt_diff_MSS","p_MSS","spearman_vs_cytolytic_score","median_MSI","median_CIN"]); w.writerows(rows)
print("\nSaved output/gene_panel_tp53.csv")
