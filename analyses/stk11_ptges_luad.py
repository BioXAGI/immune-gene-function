import json, math, statistics as st, urllib.request, urllib.error, time, csv
B="https://www.cbioportal.org/api"; S="luad_tcga_pan_can_atlas_2018"
def _open(req):
    for i in range(5):
        try: return json.load(urllib.request.urlopen(req,timeout=180))
        except urllib.error.HTTPError as e:
            if e.code<500 or i==4: raise
            time.sleep(3*(i+1))
post=lambda p,b:_open(urllib.request.Request(B+p,data=json.dumps(b).encode(),headers={"Content-Type":"application/json","Accept":"application/json"}))
get=lambda p:_open(urllib.request.Request(B+p,headers={"Accept":"application/json"}))
EXPR=["PTGES","PTGS2","PTGS1","PTGES2","PTGES3","HPGD","PTGER2","PTGER4","STK11","CD8A","GZMA","PRF1","CD274","CXCL1","CXCL2","IL6","CSF3","CXCR2","CD163","FOXP3"]
MUT=["STK11","KRAS","KEAP1","TP53"]
ent={g:get(f"/genes/{g}")["entrezGeneId"] for g in set(EXPR+MUT)}; sym={v:k for k,v in ent.items()}
E={}
for d in post(f"/molecular-profiles/{S}_rna_seq_v2_mrna/molecular-data/fetch",{"sampleListId":f"{S}_all","entrezGeneIds":[ent[g] for g in EXPR]}):
    if d.get("value") is not None: E.setdefault(sym[d["entrezGeneId"]],{})[d["sampleId"]]=math.log2(d["value"]+1)
mut={}
for m in post(f"/molecular-profiles/{S}_mutations/mutations/fetch",{"sampleListId":f"{S}_all","entrezGeneIds":[ent[g] for g in MUT]}):
    mut.setdefault(m["sampleId"],set()).add(sym[m["entrezGeneId"]])
cna={}
for d in post(f"/molecular-profiles/{S}_gistic/discrete-copy-number/fetch?discreteCopyNumberEventType=HOMDEL",{"sampleListId":f"{S}_cnaseq","entrezGeneIds":[ent["STK11"]]}):
    cna.setdefault(d["sampleId"],set()).add("STK11")
seq=set(get(f"/sample-lists/{S}_sequenced")["sampleIds"])
S_=[s for s in E["PTGES"] if s in seq]
stk=lambda s:"STK11" in mut.get(s,set()) or "STK11" in cna.get(s,set())
def ranksum(a,b):
    v=sorted([(x,0) for x in a]+[(x,1) for x in b]); r=[0]*len(v); i=0
    while i<len(v):
        j=i
        while j+1<len(v) and v[j+1][0]==v[i][0]: j+=1
        for k in range(i,j+1): r[k]=(i+j)/2+1
        i=j+1
    n1,n2=len(a),len(b); U=sum(r[k] for k in range(len(v)) if v[k][1]==0)-n1*(n1+1)/2
    return math.erfc(abs((U-n1*n2/2)/math.sqrt(n1*n2*(n1+n2+1)/12))/math.sqrt(2))
def rk(x):
    o=sorted(range(len(x)),key=lambda i:x[i]); r=[0]*len(x)
    for p,i in enumerate(o): r[i]=p
    return r
def rho(a,b):
    ra,rb=rk(a),rk(b); ma,mb=st.mean(ra),st.mean(rb)
    return sum((x-ma)*(y-mb) for x,y in zip(ra,rb))/math.sqrt(sum((x-ma)**2 for x in ra)*sum((y-mb)**2 for y in rb))
def cmp(g,A,Bb,la,lb):
    a=[E[g][s] for s in A]; b=[E[g][s] for s in Bb]
    return f"{la} n={len(a)} vs {lb} n={len(b)}: diff {st.median(a)-st.median(b):+.2f}, p={ranksum(a,b):.1e}"
mt=[s for s in S_ if stk(s)]; wt=[s for s in S_ if not stk(s)]
kr=lambda s:"KRAS" in mut.get(s,set())
print(f"LUAD samples {len(S_)}; STK11 altered {len(mt)} ({100*len(mt)/len(S_):.0f}%)\n")
print("PTGES  STK11 alt vs wt:          ", cmp("PTGES",mt,wt,"STK11","WT"))
print("PTGES  KRAS-mut only:            ", cmp("PTGES",[s for s in mt if kr(s)],[s for s in wt if kr(s)],"KL","K"))
print("PTGES  KRAS-wt only:             ", cmp("PTGES",[s for s in mt if not kr(s)],[s for s in wt if not kr(s)],"STK11","WT"))
print("PTGES  KEAP1-wt only:            ", cmp("PTGES",[s for s in mt if "KEAP1" not in mut.get(s,set())],[s for s in wt if "KEAP1" not in mut.get(s,set())],"STK11","WT"))
print("PTGES  KEAP1 mut vs wt (all):    ", cmp("PTGES",[s for s in S_ if "KEAP1" in mut.get(s,set())],[s for s in S_ if "KEAP1" not in mut.get(s,set())],"KEAP1","WT"))
print("\nOther genes, STK11 altered vs WT:")
rows=[]
for g in EXPR:
    a=[E[g][s] for s in mt]; b=[E[g][s] for s in wt]; d=st.median(a)-st.median(b); p=ranksum(a,b)
    r=rho([E["PTGES"][s] for s in S_],[E[g][s] for s in S_])
    rows.append([g,round(d,2),f"{p:.1e}",round(r,2)])
    print(f"  {g:7s} diff {d:+.2f} p={p:.1e} | rho with PTGES {r:+.2f}")
with open("output/stk11_ptges_luad.csv","w",newline="",encoding="utf-8-sig") as f:
    w=csv.writer(f); w.writerow(["gene","STK11alt_vs_wt_diff_log2","p","spearman_with_PTGES"]); w.writerows(rows)
