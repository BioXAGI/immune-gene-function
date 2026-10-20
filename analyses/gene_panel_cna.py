import json, urllib.request, urllib.error, time, statistics as st
from collections import Counter
B="https://www.cbioportal.org/api"; S="coadread_tcga_pan_can_atlas_2018"
def _open(req):
    for i in range(5):
        try: return json.load(urllib.request.urlopen(req,timeout=180))
        except urllib.error.HTTPError as e:
            if e.code<500 or i==4: raise
            time.sleep(3*(i+1))
post=lambda p,b:_open(urllib.request.Request(B+p,data=json.dumps(b).encode(),headers={"Content-Type":"application/json","Accept":"application/json"}))
get=lambda p:_open(urllib.request.Request(B+p,headers={"Accept":"application/json"}))
PANEL=["PLAGL2","REEP1","ASXL1","POFUT1","SMTNL2","CAMKV","TRPM6","WASF3","SPATA25"]
info={x:get(f"/genes/{x}") for x in PANEL+["TP53"]}
ent={x:info[x]["entrezGeneId"] for x in info}; sym={v:k for k,v in ent.items()}
cna=post(f"/molecular-profiles/{S}_gistic/discrete-copy-number/fetch?discreteCopyNumberEventType=ALL",
         {"sampleListId":f"{S}_cnaseq","entrezGeneIds":list(ent.values())})
C={}
for d in cna: C.setdefault(sym[d["entrezGeneId"]],{})[d["sampleId"]]=d["alteration"]
muts=post(f"/molecular-profiles/{S}_mutations/mutations/fetch",{"sampleListId":f"{S}_all","entrezGeneIds":[ent["TP53"]]})
tp=set(m["sampleId"] for m in muts)
print(f"{'gene':8s} {'chr':6s} {'gain% TP53mut':>13s} {'gain% TP53wt':>12s} {'loss% TP53mut':>13s} {'loss% TP53wt':>12s}")
for g in PANEL:
    d=C.get(g,{}); ss=list(d)
    m=[d[s] for s in ss if s in tp]; w=[d[s] for s in ss if s not in tp]
    pct=lambda v,f:100*sum(1 for x in v if f(x))/len(v)
    ch=info[g].get("chromosome") or ""
    print(f"{g:8s} {str(ch):6s} {pct(m,lambda x:x>0):13.0f} {pct(w,lambda x:x>0):12.0f} {pct(m,lambda x:x<0):13.0f} {pct(w,lambda x:x<0):12.0f}")
