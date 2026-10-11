"""按主题检索真实研究型文献（Europe PMC REST，剔除综述），落 JSON + 生成 GB/T 7714 风格清单。"""
import json, os, time, urllib.parse, urllib.request

ROOT = "/mnt/f/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
OUT = os.path.join(ROOT, "稿")
os.makedirs(OUT, exist_ok=True)
API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

TOPICS = {
 "exosome_glioma_model": '(glioma OR glioblastoma) AND (exosome OR exosomal) AND (prognostic model OR risk score) AND (SRC:MED)',
 "exosome_liquid_biopsy": '(glioma OR glioblastoma) AND extracellular vesicle AND (biomarker OR liquid biopsy) AND (SRC:MED)',
 "exosome_machinery": '(RAB27A OR exosome secretion) AND (tumour OR tumor OR cancer) AND (SRC:MED)',
 "shap_interpretability": '(SHAP OR "Shapley additive explanations") AND (cancer OR survival OR prognosis) AND (machine learning) AND (SRC:MED)',
 "idh_glioma": 'IDH mutation AND glioma AND (prognosis OR survival) AND (SRC:MED)',
 "mgmt_temozolomide": 'MGMT promoter methylation AND glioblastoma AND temozolomide AND (SRC:MED)',
 "who_cns5": '"WHO CNS5" OR "WHO 2021 classification" AND glioma AND molecular AND (SRC:MED)',
 "nomogram_glioma": 'nomogram AND glioma AND (overall survival OR prediction model) AND (SRC:MED)',
 "ezh2_glioma": 'EZH2 AND (glioma OR glioblastoma) AND (proliferation OR prognosis) AND (SRC:MED)',
 "glioma_tme": 'glioma AND (tumour microenvironment OR tumor microenvironment) AND (macrophage OR immune infiltration) AND (prognosis) AND (SRC:MED)',
 "exosome_tcga_signature": 'exosome-related genes AND (TCGA OR prognostic signature) AND cancer AND (SRC:MED)',
 "lasso_cox": 'lasso AND Cox AND (prognostic signature OR risk model) AND cancer AND (SRC:MED)',
 "neuronal_glioma": '(neuron-glioma OR synaptic) AND glioma AND (progression OR plasticity) AND (SRC:MED)',
 "glioma_biomarker_validation": 'glioma AND prognostic signature AND external validation AND (CGGA OR "Chinese Glioma Genome Atlas") AND (SRC:MED)',
}

def search(query, page_size=10):
    q = f'({query}) AND (PUB_TYPE:"Journal Article") NOT (PUB_TYPE:"Review") NOT (PUB_TYPE:"Editorial")'
    url = f"{API}?query={urllib.parse.quote(q)}&format=json&resultType=core&pageSize={page_size}&sort=CITED%20desc"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r).get("resultList", {}).get("result", [])

def fmt(a):
    auth = a.get("authorString", "").rstrip(".")
    parts = [x.strip() for x in auth.split(",") if x.strip()]
    if len(parts) > 3:
        auth = ", ".join(parts[:3]) + ", et al"
    j = a.get("journalInfo", {}).get("journal", {}).get("title") or a.get("bookOrReportDetails", {}).get("publisher", "")
    vol = a.get("journalInfo", {}).get("volume", "")
    pages = a.get("pageInfo", "")
    return {
        "pmid": a.get("pmid", ""), "doi": a.get("doi", ""), "title": a.get("title", "").strip(),
        "journal": j, "year": a.get("pubYear", ""), "volume": vol, "pages": pages,
        "authors": auth, "type": ",".join(a.get("pubTypeList", {}).get("pubType", []) if isinstance(a.get("pubTypeList", {}).get("pubType", []), list) else []) ,
        "cited": a.get("citedByCount", 0),
        "ref": f"{auth}. {a.get('title','').strip().rstrip('.')}. {j}, {a.get('pubYear','')}"
               + (f", {vol}" if vol else "") + (f": {pages}" if pages else "")
               + (f". doi: {a.get('doi')}" if a.get("doi") else "") + (f". PMID: {a.get('pmid')}" if a.get("pmid") else ""),
    }

store = {}
for key, q in TOPICS.items():
    try:
        res = search(q)
        items = [fmt(x) for x in res if x.get("pmid")]
        # 去重：同标题只留引用最高的
        seen, uniq = set(), []
        for it in items:
            t = it["title"].lower()[:60]
            if t in seen: continue
            seen.add(t); uniq.append(it)
        store[key] = uniq
        print(f"[{key}] {len(uniq)} 篇")
        time.sleep(0.4)
    except Exception as e:
        print(f"[{key}] ERROR {type(e).__name__} {e}")
        store[key] = []

json.dump(store, open(os.path.join(OUT, "参考文献_原始.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
n = sum(len(v) for v in store.values())
print("总计", n, "篇")
