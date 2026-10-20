import json, urllib.request, urllib.parse, time, xml.etree.ElementTree as ET
E="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
def q(u):
    for i in range(4):
        try: return urllib.request.urlopen(u,timeout=60).read()
        except Exception: time.sleep(2*(i+1))
    return b""
Q={"STK11_PGE2":'(STK11[tiab] OR LKB1[tiab]) AND (PTGES[tiab] OR mPGES-1[tiab] OR PGE2[tiab] OR "prostaglandin E2"[tiab] OR COX-2[tiab] OR PTGS2[tiab] OR cyclooxygenase[tiab])',
   "PTGES_lung":'(PTGES[tiab] OR mPGES-1[tiab] OR "microsomal prostaglandin E synthase"[tiab]) AND (lung[tiab]) AND (cancer[tiab] OR carcinoma[tiab] OR tumor[tiab])',
   "STK11_immune":'(STK11[tiab] OR LKB1[tiab]) AND (lung[tiab]) AND (neutrophil*[tiab] OR "immune checkpoint"[tiab] OR "cold"[tiab] OR immunosuppress*[tiab])',
   "PGE2_ICB":'(PGE2[tiab] OR "prostaglandin E2"[tiab] OR EP2[tiab] OR EP4[tiab]) AND ("anti-PD-1"[tiab] OR "checkpoint"[tiab]) AND (cancer[tiab])',
   "KEAP1_PGE2":'(KEAP1[tiab] OR NRF2[tiab] OR NFE2L2[tiab]) AND (PTGES[tiab] OR mPGES-1[tiab] OR PGE2[tiab]) AND (cancer[tiab] OR tumor[tiab])'}
out={}
for k,t in Q.items():
    r=json.loads(q(f"{E}/esearch.fcgi?db=pubmed&retmode=json&retmax=7&sort=relevance&term="+urllib.parse.quote(t)))
    ids=r["esearchresult"]["idlist"]; arts=[]
    if ids:
        x=ET.fromstring(q(f"{E}/efetch.fcgi?db=pubmed&retmode=xml&id={','.join(ids)}"))
        for a in x.findall(".//PubmedArticle"):
            arts.append({"pmid":a.findtext(".//PMID"),"year":(a.findtext(".//PubDate/Year") or a.findtext(".//PubDate/MedlineDate") or "")[:4],
              "journal":a.findtext(".//ISOAbbreviation") or "","title":"".join(a.find(".//ArticleTitle").itertext()),
              "abstract":" ".join("".join(s.itertext()) for s in a.findall(".//AbstractText"))})
    out[k]={"count":r["esearchresult"]["count"],"top":arts}; time.sleep(0.4)
json.dump(out,open("output/pubmed_stk11_ptges.json","w",encoding="utf-8"),indent=1,ensure_ascii=False)
for k,v in out.items():
    print(f"\n##### {k} (total {v['count']})")
    for a in v["top"]: print(f"- {a['pmid']} {a['year']} {a['journal']} | {a['title']}\n   {a['abstract'][-550:]}")
