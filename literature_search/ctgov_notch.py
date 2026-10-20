import json, urllib.request, urllib.parse, csv
API="https://clinicaltrials.gov/api/v2/studies"
TERMS=["Notch","gamma-secretase","DLL4","DLL3","JAG1","JAG2","navicixizumab","nirogacestat","crenigacestat","AL102","CB-103",
       "demcizumab","tarextumab","brontictuzumab","tarlatamab","fucose","POFUT1","RO4929097","PF-03084014","MK-0752","LY3039478","BMS-906024"]
def fetch(params):
    out=[]; tok=None
    while True:
        p=dict(params,pageSize=100,format="json")
        if tok: p["pageToken"]=tok
        d=json.load(urllib.request.urlopen(API+"?"+urllib.parse.urlencode(p),timeout=60))
        out+=d.get("studies",[]); tok=d.get("nextPageToken")
        if not tok: return out
seen={}
for t in TERMS:
    for s in fetch({"query.cond":"colorectal OR colon OR rectal OR gastrointestinal OR \"solid tumor\" OR \"advanced cancer\"","query.term":t}):
        ps=s["protocolSection"]; nct=ps["identificationModule"]["nctId"]
        seen.setdefault(nct,{"s":ps,"hits":set()})["hits"].add(t)
CRC=("colorectal","colon","rectal","rectum","bowel")
rows=[]
for nct,v in seen.items():
    ps=v["s"]; idm=ps["identificationModule"]; st=ps["statusModule"]; cm=ps.get("conditionsModule",{})
    conds=cm.get("conditions",[]); title=idm.get("briefTitle","")
    arms=ps.get("armsInterventionsModule",{}).get("interventions",[])
    ivs="; ".join(f'{i.get("name")}' for i in arms)
    blob=(title+" "+" ".join(conds)+" "+" ".join(cm.get("keywords",[]))).lower()
    elig=ps.get("eligibilityModule",{}).get("eligibilityCriteria","").lower()
    crc_specific=any(k in blob for k in CRC); crc_mentioned=crc_specific or any(k in elig for k in CRC)
    txt=(title+" "+ivs+" "+ps.get("descriptionModule",{}).get("briefSummary","")).lower()
    notch_rel=any(k.lower() in txt for k in TERMS)
    if not notch_rel or not crc_mentioned: continue
    rows.append({"nct":nct,"status":st.get("overallStatus"),"phase":",".join(ps.get("designModule",{}).get("phases",[]) or []),
        "start":st.get("startDateStruct",{}).get("date",""),"title":title,"interventions":ivs,"conditions":"; ".join(conds[:6]),
        "crc_in_title_or_conditions":crc_specific,"sponsor":ps.get("sponsorCollaboratorsModule",{}).get("leadSponsor",{}).get("name",""),
        "matched_terms":", ".join(sorted(v["hits"]))})
ACTIVE={"RECRUITING","NOT_YET_RECRUITING","ACTIVE_NOT_RECRUITING","ENROLLING_BY_INVITATION"}
rows.sort(key=lambda r:(r["status"] not in ACTIVE, not r["crc_in_title_or_conditions"], r["start"]),reverse=False)
with open("output/ctgov_notch_crc.csv","w",newline="",encoding="utf-8-sig") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
for r in rows:
    flag="ACTIVE" if r["status"] in ACTIVE else r["status"]
    print(f'{r["nct"]} | {flag} | {r["phase"]} | {r["start"]} | CRC-specific={r["crc_in_title_or_conditions"]} | {r["sponsor"]}\n   {r["title"]}\n   -> {r["interventions"][:160]}')
print(f"\n{len(rows)} trials; active: {sum(r['status'] in ACTIVE for r in rows)}")
