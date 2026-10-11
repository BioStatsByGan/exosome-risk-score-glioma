"""第三轮文献扩充：为前言与讨论的每个论点准备足量真实研究型文献。"""
import json, os, time, urllib.parse, urllib.request

ROOT = "/mnt/f/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
TOPICS = {
 "exosome_sig_pancancer": '(TITLE:"exosome" OR TITLE:"exosomal" OR TITLE:"extracellular vesicle") AND (TITLE:"prognostic signature" OR TITLE:"risk signature" OR TITLE:"risk model" OR TITLE:"prognostic model")',
 "exosomal_mirna_glioma": '(TITLE:"exosomal microRNA" OR TITLE:"exosomal miRNA" OR TITLE:"serum exosomal") AND (TITLE:"glioma" OR TITLE:"glioblastoma")',
 "exosome_gbm_immune": '(TITLE:"extracellular vesicles" OR TITLE:"exosomes") AND (TITLE:"glioblastoma" OR TITLE:"glioma") AND (TITLE:"immune" OR TITLE:"macrophage" OR TITLE:"microglia" OR TITLE:"immunosuppression")',
 "exosome_treatment": '(TITLE:"exosome" OR TITLE:"extracellular vesicles") AND (TITLE:"blood-brain barrier" OR TITLE:"drug delivery" OR TITLE:"therapy") AND (TITLE:"glioma" OR TITLE:"glioblastoma" OR TITLE:"brain")',
 "mozin_ezrin_cancer": '(TITLE:"moesin" OR TITLE:"ezrin" OR TITLE:"ERM proteins") AND (TITLE:"cancer" OR TITLE:"glioma" OR TITLE:"invasion" OR TITLE:"metastasis")',
 "ppia_cancer": '(TITLE:"cyclophilin A" OR TITLE:"PPIA" OR TITLE:"peptidyl-prolyl isomerase") AND (TITLE:"cancer" OR TITLE:"tumour" OR TITLE:"tumor" OR TITLE:"glioma")',
 "eif2s1_cancer": '(TITLE:"eIF2alpha" OR TITLE:"eIF2α" OR TITLE:"EIF2S1" OR TITLE:"integrated stress response") AND (TITLE:"cancer" OR TITLE:"glioma" OR TITLE:"tumour" OR TITLE:"tumor")',
 "cct3_cancer": '(TITLE:"CCT3" OR TITLE:"chaperonin containing TCP1" OR TITLE:"TRiC") AND (TITLE:"cancer" OR TITLE:"tumour" OR TITLE:"tumor" OR TITLE:"prognosis")',
 "snx_cancer": '(TITLE:"sorting nexin" OR TITLE:"SNX18") AND (TITLE:"cancer" OR TITLE:"endosomal" OR TITLE:"tumour" OR TITLE:"tumor")',
 "aqp1_cancer": '(TITLE:"aquaporin 1" OR TITLE:"AQP1") AND (TITLE:"glioma" OR TITLE:"cancer" OR TITLE:"tumour" OR TITLE:"tumor")',
 "myc_glioma": '(TITLE:"MYC") AND (TITLE:"glioma" OR TITLE:"glioblastoma") AND (TITLE:"prognosis" OR TITLE:"progression")',
 "utp20_ribosome": '(TITLE:"UTP20" OR TITLE:"ribosome biogenesis") AND (TITLE:"cancer" OR TITLE:"glioma" OR TITLE:"prognosis" OR TITLE:"tumour")',
 "gng12_cancer": '(TITLE:"GNG12" OR TITLE:"G protein gamma") AND (TITLE:"cancer" OR TITLE:"tumour" OR TITLE:"tumor" OR TITLE:"glioma")',
 "adgre5_cancer": '(TITLE:"ADGRE5" OR TITLE:"CD97") AND (TITLE:"cancer" OR TITLE:"glioma" OR TITLE:"tumour" OR TITLE:"tumor")',
 "lrrk2_cancer": '(TITLE:"LRRK2") AND (TITLE:"cancer" OR TITLE:"glioma" OR TITLE:"tumour" OR TITLE:"tumor")',
 "acy3_cancer": '(TITLE:"aminoacylase" OR TITLE:"ACY1" OR TITLE:"ACY3") AND (TITLE:"cancer" OR TITLE:"tumour" OR TITLE:"tumor")',
 "ezh2_broad": '(TITLE:"EZH2") AND (TITLE:"prognosis" OR TITLE:"prognostic") AND (TITLE:"cancer" OR TITLE:"tumour" OR TITLE:"tumor")',
 "stupp_standard": '(TITLE:"radiotherapy plus concomitant and adjuvant temozolomide") OR (TITLE:"Stupp") AND (TITLE:"glioblastoma")',
 "gbm_recurrence": '(TITLE:"glioblastoma") AND (TITLE:"recurrence" OR TITLE:"relapse") AND (TITLE:"survival" OR TITLE:"outcome")',
 "idh_wildtype_gbm": '(TITLE:"IDH-wildtype" OR TITLE:"IDH wild-type") AND (TITLE:"glioblastoma") AND (TITLE:"outcome" OR TITLE:"prognosis" OR TITLE:"survival")',
 "model_validation_methods": '(TITLE:"decision curve analysis" OR TITLE:"net reclassification" OR TITLE:"calibration") AND (TITLE:"prediction model" OR TITLE:"prognostic model")',
 "ml_survival_glioma": '(TITLE:"machine learning") AND (TITLE:"glioma" OR TITLE:"glioblastoma") AND (TITLE:"prognosis" OR TITLE:"survival")',
 "exosome_ncrna_cancer": '(TITLE:"exosomal" OR TITLE:"exosome") AND (TITLE:"lncRNA" OR TITLE:"circRNA" OR TITLE:"long non-coding") AND (TITLE:"prognosis" OR TITLE:"prognostic")',
 "glioma_grade_prognosis": '(TITLE:"WHO grade" OR TITLE:"histological grade") AND (TITLE:"glioma") AND (TITLE:"survival" OR TITLE:"prognosis")',
 "exosome_isolation_biomarker": '(TITLE:"liquid biopsy") AND (TITLE:"glioma" OR TITLE:"glioblastoma" OR TITLE:"brain tumour" OR TITLE:"brain tumor")',
}

def search(query, size=9):
    q = f'({query}) AND (SRC:MED) NOT (PUB_TYPE:"Review") NOT (PUB_TYPE:"Editorial")'
    url = (f"{API}?query={urllib.parse.quote(q)}&format=json&resultType=core"
           f"&pageSize={size}&sort=CITED%20desc")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r).get("resultList", {}).get("result", [])

def fmt(a):
    auth = a.get("authorString", "").rstrip(".")
    parts = [x.strip() for x in auth.split(",") if x.strip()]
    if len(parts) > 3: auth = ", ".join(parts[:3]) + ", et al"
    ji = a.get("journalInfo", {}); j = ji.get("journal", {}).get("title", "")
    vol = ji.get("volume", ""); pages = a.get("pageInfo", "")
    return {"pmid": a.get("pmid", ""), "doi": a.get("doi", ""), "title": a.get("title", "").strip(),
            "journal": j, "year": a.get("pubYear", ""), "volume": vol, "pages": pages,
            "authors": auth, "cited": a.get("citedByCount", 0),
            "ref": f"{auth}. {a.get('title','').strip().rstrip('.')}. {j}, {a.get('pubYear','')}"
                   + (f", {vol}" if vol else "") + (f": {pages}" if pages else "")
                   + (f". doi: {a.get('doi')}" if a.get("doi") else "") + (f". PMID: {a.get('pmid')}" if a.get("pmid") else "")}

out = {}
for k, q in TOPICS.items():
    try:
        items = [fmt(x) for x in search(q) if x.get("pmid")]
        seen, uniq = set(), []
        for it in items:
            t = it["title"].lower()[:55]
            if t in seen: continue
            seen.add(t); uniq.append(it)
        out[k] = uniq
        print(f"\n=== {k} ({len(uniq)}) ===")
        for it in uniq[:6]:
            print(f"  {it['pmid']} | {it['year']} | {it['journal'][:26]:26s} | {it['title'][:86]}")
    except Exception as e:
        print(f"[{k}] ERR {e}")
        out[k] = []
    time.sleep(0.3)

json.dump(out, open(os.path.join(ROOT, "稿", "参考文献_扩充.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
print("\n总计", sum(len(v) for v in out.values()), "篇")
