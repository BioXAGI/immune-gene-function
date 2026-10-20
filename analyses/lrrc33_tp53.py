import gzip, json, math, statistics as st, urllib.request, urllib.error, time
G=["LRRC33","TGFB1","LRRC32","ITGB8","ITGAV","CD68","CD163","CD14","ITGAM","FOXP3","CD8A","SMAD7","SERPINE1"]
X={}
with gzip.open("output/coadread_hiseqv2.gz","rt") as f:
    hdr=f.readline().rstrip("\n").split("\t")[1:]
    for line in f:
        p=line.rstrip("\n").split("\t")
        if p[0] in G: X[p[0]]=dict(zip(hdr,map(float,p[1:])))
print("found:",sorted(X))
tum=[s for s in hdr if s[13:15]=="01"]; nor=[s for s in hdr if s[13:15]=="11"]
pairs=[(s[:12]+"-01",s) for s in nor if s[:12]+"-01" in hdr]
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
    return math.erfc(abs((U-n1*n2/2)/math.sqrt(n1*n2*(n1+n2+1)/12))/math.sqrt(2))
def signrank(d):
    d=[x for x in d if x!=0]; n=len(d); o=sorted(range(n),key=lambda i:abs(d[i])); r=[0]*n
    for p,i in enumerate(o): r[i]=p+1
    W=sum(r[i] for i in range(n) if d[i]>0); return math.erfc(abs((W-n*(n+1)/4)/math.sqrt(n*(n+1)*(2*n+1)/24))/math.sqrt(2))
def rk(x):
    o=sorted(range(len(x)),key=lambda i:x[i]); r=[0]*len(x)
    for p,i in enumerate(o): r[i]=p
    return r
def rho(a,b):
    ra,rb=rk(a),rk(b); ma,mb=st.mean(ra),st.mean(rb)
    return sum((x-ma)*(y-mb) for x,y in zip(ra,rb))/math.sqrt(sum((x-ma)**2 for x in ra)*sum((y-mb)**2 for y in rb))
T=[s for s in tum if s[:12] in seq]; mut=[s for s in T if s[:12] in tp]; wt=[s for s in T if s[:12] not in tp]
mss=lambda L:[s for s in L if sub.get(s[:12],"").endswith("CIN")]
x=X["LRRC33"]; d=[x[t]-x[n] for t,n in pairs]
print(f"\nNRROS(LRRC33): normal {st.median(x[s] for s in nor):.2f} tumor {st.median(x[s] for s in tum):.2f}; paired d {st.median(d):+.2f} ({100*sum(v>0 for v in d)/len(d):.0f}% up) p={signrank(d):.1e}")
for g in ["LRRC33","TGFB1","LRRC32","ITGB8"]:
    a=[X[g][s] for s in mut]; b=[X[g][s] for s in wt]; am=[X[g][s] for s in mss(mut)]; bm=[X[g][s] for s in mss(wt)]
    print(f"{g:7s} TP53mut-wt all {st.median(a)-st.median(b):+.2f} p={ranksum(a,b):.1e} | MSS {st.median(am)-st.median(bm):+.2f} p={ranksum(am,bm):.1e}")
print("\nSpearman with NRROS across tumors:")
for g in ["TGFB1","CD68","CD163","CD14","ITGAM","ITGB8","LRRC32","FOXP3","CD8A","SERPINE1","SMAD7"]:
    print(f"  {g:8s} {rho([X["LRRC33"][s] for s in tum],[X[g][s] for s in tum]):+.2f}")
