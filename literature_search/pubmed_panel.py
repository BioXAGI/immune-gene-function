import json, urllib.request, urllib.parse, time, xml.etree.ElementTree as ET
E="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
def q(url):
    for i in range(4):
        try: return urllib.request.urlopen(url,timeout=60).read()
        except Exception:
            time.sleep(2*(i+1))
    return b""
GENES={"PLAGL2":"PLAGL2","REEP1":"REEP1","ASXL1":"ASXL1","POFUT1":"POFUT1","SMTNL2":"SMTNL2",
       "CAMKV":"CAMKV","TRPM6":"TRPM6","WASF3":"WASF3 OR WAVE3","SPATA25":"SPATA25"}
CTX={"colorectal":'(colorectal OR colon OR rectal)',"immune":'(immune OR immunity OR "T cell" OR macrophage OR "tumor microenvironment" OR immunotherapy)',"TP53":'(TP53 OR p53)'}
out={}
for g,term in GENES.items():
    out[g]={}
    for c,ct in CTX.items():
        full=f'({term}[tiab]) AND {ct}[tiab]' if " OR " not in term else f'(({term.replace(" OR ","[tiab] OR ")}[tiab])) AND {ct}[tiab]'
        r=json.loads(q(f"{E}/esearch.fcgi?db=pubmed&retmode=json&retmax=6&sort=relevance&term="+urllib.parse.quote(full)))
        ids=r["esearchresult"]["idlist"]; n=r["esearchresult"]["count"]
        arts=[]
        if ids:
            x=ET.fromstring(q(f"{E}/efetch.fcgi?db=pubmed&retmode=xml&id={','.join(ids)}"))
            for a in x.findall(".//PubmedArticle"):
                t="".join(a.find(".//ArticleTitle").itertext())
                y=(a.findtext(".//PubDate/Year") or a.findtext(".//PubDate/MedlineDate") or "")[:4]
                j=a.findtext(".//ISOAbbreviation") or ""
                pm=a.findtext(".//PMID")
                ab=" ".join("".join(s.itertext()) for s in a.findall(".//AbstractText"))
                arts.append({"pmid":pm,"year":y,"journal":j,"title":t,"abstract":ab[:900]})
        out[g][c]={"count":n,"top":arts}
        time.sleep(0.4)
    print(g,{c:out[g][c]["count"] for c in CTX})
json.dump(out,open("output/pubmed_panel.json","w",encoding="utf-8"),indent=1,ensure_ascii=False)
