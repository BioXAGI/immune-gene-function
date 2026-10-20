import gzip, json, math, statistics as st, urllib.request, urllib.error, time, csv
LIG=["DLL1","DLL3","DLL4","JAG1","JAG2"]; REC=["NOTCH1","NOTCH2","NOTCH3","NOTCH4"]
TGT=["HES1","HEY1","NRARP","ATOH1","POFUT1"]; GENES=LIG+REC+TGT
X={}
with gzip.open("output/coadread_hiseqv2.gz","rt") as f:
    hdr=f.readline().rstrip("\n").split("\t")[1:]
    for line in f:
        p=line.rstrip("\n").split("\t")
        if p[0] in GENES: X[p[0]]=dict(zip(hdr,map(float,p[1:])))
tum=[s for s in hdr if s[13:15]=="01"]; nor=[s for s in hdr if s[13:15]=="11"]
pairs=[(s[:12]+"-01",s) for s in nor if s[:12]+"-01" in X["DLL1"]]

B="https://www.cbioportal.org/api"; S="coadread_tcga_pan_can_atlas_2018"
def _open(req):
    for i in range(5):
        try: return json.load(urllib.request.urlopen(req,timeout=180))
        except urllib.error.HTTPError as e:
            if e.code<500 or i==4: raise
            time.sleep(3*(i+1))
post=lambda p,b:_open(urllib.request.Request(B+p,data=json.dumps(b).encode(),headers={"Content-Type":"application/json","Accept":"application/json"}))
get=lambda p:_open(urllib.request.Request(B+p,headers={"Accept":"application/json"}))
seq={s[:12] for s in get(f"/sample-lists/{S}_sequenced")["sampleIds"]}
tp={m["sampleId"][:12] for m in post(f"/molecular-profiles/{S}_mutations/mutations/fetch",{"sampleListId":f"{S}_all","entrezGeneIds":[7157]})}
sub={d["patientId"]:d["value"] for d in get(f"/studies/{S}/clinical-data?clinicalDataType=PATIENT&attributeId=SUBTYPE&pageSize=10000")}

def ranksum(a,b):
    v=sorted([(x,0) for x in a]+[(x,1) for x in b]); r=[0]*len(v); i=0
    while i<len(v):
        j=i
        while j+1<len(v) and v[j+1][0]==v[i][0]: j+=1
        for k in range(i,j+1): r[k]=(i+j)/2+1
        i=j+1
    n1,n2=len(a),len(b); U=sum(r[k] for k in range(len(v)) if v[k][1]==0)-n1*(n1+1)/2
    z=(U-n1*n2/2)/math.sqrt(n1*n2*(n1+n2+1)/12); return math.erfc(abs(z)/math.sqrt(2))
def signrank(d):
    d=[x for x in d if x!=0]; n=len(d); o=sorted(range(n),key=lambda i:abs(d[i])); r=[0]*n
    for p,i in enumerate(o): r[i]=p+1
    W=sum(r[i] for i in range(n) if d[i]>0); mu=n*(n+1)/4; sd=math.sqrt(n*(n+1)*(2*n+1)/24)
    return math.erfc(abs((W-mu)/sd)/math.sqrt(2))

T=[s for s in tum if s[:12] in seq]
mut=[s for s in T if s[:12] in tp]; wt=[s for s in T if s[:12] not in tp]
mss=lambda L:[s for s in L if sub.get(s[:12],"").endswith("CIN")]
print(f"Tumors {len(tum)} | normals {len(nor)} | matched pairs {len(pairs)} | sequenced tumors TP53mut {len(mut)} / wt {len(wt)}\n")
print(f"{'gene':7s} {'normal':>6s} {'tumor':>6s} {'pair d':>7s} {'%pairs up':>9s} {'p':>8s} | {'TP53 d':>6s} {'p':>8s} | {'MSS d':>6s} {'p':>8s}")
rows=[]
for g in GENES:
    x=X[g]; d=[x[t]-x[n] for t,n in pairs]
    pn=signrank(d); up=100*sum(1 for v in d if v>0)/len(d)
    a=[x[s] for s in mut]; b=[x[s] for s in wt]; p1=ranksum(a,b)
    am=[x[s] for s in mss(mut)]; bm=[x[s] for s in mss(wt)]; p2=ranksum(am,bm)
    r=[g,st.median(x[s] for s in nor),st.median(x[s] for s in tum),st.median(d),up,pn,
       st.median(a)-st.median(b),p1,st.median(am)-st.median(bm),p2]
    rows.append(r)
    print(f"{g:7s} {r[1]:6.2f} {r[2]:6.2f} {r[3]:+7.2f} {r[4]:8.0f}% {r[5]:8.1e} | {r[6]:+6.2f} {r[7]:8.1e} | {r[8]:+6.2f} {r[9]:8.1e}")
with open("output/notch_tp53.csv","w",newline="",encoding="utf-8-sig") as f:
    w=csv.writer(f); w.writerow(["gene","median_normal","median_tumor","paired_median_diff","pct_pairs_tumor_higher","p_paired","TP53mut_vs_wt_diff","p_TP53","TP53_diff_MSS_only","p_TP53_MSS"])
    for r in rows: w.writerow([r[0]]+[round(v,3) if isinstance(v,float) and abs(v)>=1e-3 else f"{v:.2e}" if isinstance(v,float) else v for v in r[1:]])
print("\nValues: log2(RSEM+1). Saved output/notch_tp53.csv")
