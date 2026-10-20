import json, urllib.request, urllib.parse, time, xml.etree.ElementTree as ET
E="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
def q(u):
    for i in range(4):
        try: return urllib.request.urlopen(u,timeout=60).read()
        except Exception: time.sleep(2*(i+1))
    return b""
Q={
 "thymus":'(DLL4[tiab] OR "Delta-like 4"[tiab]) AND (thymus[tiab] OR thymic[tiab] OR "T cell development"[tiab] OR "T-cell development"[tiab])',
 "effector_T_GVHD":'(DLL4[tiab] OR "Delta-like 4"[tiab]) AND ("graft-versus-host"[tiab] OR alloreactive[tiab] OR "fibroblastic reticular"[tiab])',
 "Tfh_germinal":'(DLL4[tiab] OR "Delta-like 4"[tiab]) AND ("follicular helper"[tiab] OR "germinal center"[tiab] OR Tfh[tiab])',
 "tumor_immunity":'(DLL4[tiab] OR "Delta-like 4"[tiab]) AND (tumor[tiab] OR cancer[tiab]) AND ("T cell"[tiab] OR "T cells"[tiab] OR macrophage*[tiab] OR "myeloid-derived"[tiab] OR "regulatory T"[tiab] OR immunosuppress*[tiab])',
 "PD1_combo":'(DLL4[tiab] OR "Delta-like 4"[tiab] OR demcizumab[tiab] OR navicixizumab[tiab]) AND (PD-1[tiab] OR PD-L1[tiab] OR pembrolizumab[tiab] OR "checkpoint"[tiab])',
 "clinical_safety":'(demcizumab[tiab] OR navicixizumab[tiab] OR enoticumab[tiab] OR "anti-DLL4"[tiab]) AND (phase[tiab] OR patients[tiab])',
}
out={}
for k,t in Q.items():
    r=json.loads(q(f"{E}/esearch.fcgi?db=pubmed&retmode=json&retmax=8&sort=relevance&term="+urllib.parse.quote(t)))
    ids=r["esearchresult"]["idlist"]; arts=[]
    if ids:
        x=ET.fromstring(q(f"{E}/efetch.fcgi?db=pubmed&retmode=xml&id={','.join(ids)}"))
        for a in x.findall(".//PubmedArticle"):
            arts.append({"pmid":a.findtext(".//PMID"),"year":(a.findtext(".//PubDate/Year") or a.findtext(".//PubDate/MedlineDate") or "")[:4],
              "journal":a.findtext(".//ISOAbbreviation") or "","title":"".join(a.find(".//ArticleTitle").itertext()),
              "abstract":" ".join("".join(s.itertext()) for s in a.findall(".//AbstractText"))})
    out[k]={"count":r["esearchresult"]["count"],"top":arts}; time.sleep(0.4)
json.dump(out,open("output/pubmed_dll4_immune.json","w",encoding="utf-8"),indent=1,ensure_ascii=False)
for k,v in out.items():
    print(f"\n##### {k} (total {v['count']})")
    for a in v["top"]:
        print(f"- {a['pmid']} {a['year']} {a['journal']} | {a['title']}\n   {a['abstract'][-650:]}")
