import json, urllib.request, urllib.parse, time, xml.etree.ElementTree as ET
E="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
def q(u):
    for i in range(4):
        try: return urllib.request.urlopen(u,timeout=60).read()
        except Exception: time.sleep(2*(i+1))
    return b""
Q={"KRAS_cFLIP":'(KRAS[tiab] OR "oncogenic Ras"[tiab]) AND (cFLIP[tiab] OR CFLAR[tiab] OR "FLICE-inhibitory"[tiab]) AND (TNF[tiab] OR "TNF-related"[tiab] OR TRAIL[tiab] OR apoptosis[tiab])',
"KRAS_IAP_SMAC":'(KRAS[tiab] OR "RAS mutation"[tiab]) AND ("SMAC mimetic"[tiab] OR "IAP antagonist"[tiab] OR birinapant[tiab] OR cIAP[tiab] OR BIRC2[tiab] OR BIRC3[tiab])',
"KRAS_TNF_death_CRC":'(KRAS[tiab]) AND (colorectal[tiab] OR colon[tiab]) AND ("TNF-induced"[tiab] OR TNFα[tiab] OR TNF[tiab]) AND (necroptosis[tiab] OR apoptosis[tiab] OR "cell death"[tiab])',
"KRAS_RIPK1_TNF":'(KRAS[tiab] OR "oncogenic RAS"[tiab]) AND (RIPK1[tiab] OR necroptosis[tiab]) AND (TNF[tiab])'}
out={}
for k,t in Q.items():
    r=json.loads(q(f"{E}/esearch.fcgi?db=pubmed&retmode=json&retmax=8&sort=relevance&term="+urllib.parse.quote(t)))
    ids=r["esearchresult"]["idlist"]; arts=[]
    if ids:
        x=ET.fromstring(q(f"{E}/efetch.fcgi?db=pubmed&retmode=xml&id={','.join(ids)}"))
        for a in x.findall(".//PubmedArticle"):
            arts.append({"pmid":a.findtext(".//PMID"),"year":(a.findtext(".//PubDate/Year") or a.findtext(".//PubDate/MedlineDate") or "")[:4],"journal":a.findtext(".//ISOAbbreviation") or "","title":"".join(a.find(".//ArticleTitle").itertext()),"abstract":" ".join("".join(s.itertext()) for s in a.findall(".//AbstractText"))})
    out[k]={"count":r["esearchresult"]["count"],"top":arts}; time.sleep(0.4)
json.dump(out,open("pubmed_kras_tnf.json","w",encoding="utf-8"),indent=1,ensure_ascii=False)
for k,v in out.items():
    print(f"\n#### {k} (total {v['count']})")
    for a in v["top"]: print(f"- {a['pmid']} {a['year']} {a['journal']} | {a['title']}\n   {a['abstract'][-420:]}")
