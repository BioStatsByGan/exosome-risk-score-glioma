"""精确取文献：先按标题核实地标文献，再按窄主题补足（Europe PMC）。"""
import json, os, time, urllib.parse, urllib.request

OUT = "/mnt/f/0000_Bioinfo/222_tcgaSHAP_胶质瘤/稿"
API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

LANDMARKS = [
 "Comprehensive genomic characterization defines human glioblastoma genes and core pathways",
 "IDH1 and IDH2 mutations in gliomas",
 "MGMT gene silencing and benefit from temozolomide in glioblastoma",
 "The 2021 WHO Classification of Tumors of the Central Nervous System: a summary",
 "The somatic genomic landscape of glioblastoma",
 "Molecular Profiling Reveals Biologically Discrete Subsets and Pathways of Progression in Diffuse Glioma",
 "Comprehensive, Integrative Genomic Analysis of Diffuse Lower-Grade Gliomas",
 "Cancer exosomes perform cell-independent microRNA biogenesis and promote tumorigenesis",
 "Tumour exosome integrins determine organotropic metastasis",
 "Rab27a and Rab27b control different steps of the exosome secretion pathway",
 "Electrical and synaptic integration of glioma into neural circuits",
 "A unifying explanation for diverse metabolic regulation by 2-hydroxyglutarate",
 "MBTD1 is associated with prognosis and immune infiltration in glioma",
 "Tumor-associated macrophages in glioma: friend or foe",
 "TRIPOD statement",
 "A unified approach to interpreting model predictions",
 "From local explanations to global understanding with explainable AI for trees",
 "Explaining prediction models and individual predictions with feature contributions",
 "Chinese Glioma Genome Atlas (CGGA): A Comprehensive Resource with Functional Genomic Data from Chinese Glioma Patients",
 "An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analysis",
 "Extracellular vesicles in cancer: exosomes, microvesicles and the emerging role of large oncosomes",
 "Minimal information for studies of extracellular vesicles 2018 (MISEV2018)",
 "Glioblastoma multiforme: a review of where we have been and where we are going",
 "Immune and stromal cell associations and prognosis in glioma",
 "Somatic mutations in the chromatin remodeling pathway in glioblastoma",
]

TOPICS = {
 "exosome_signature_glioma": '(TITLE:"exosome" OR TITLE:"exosomal" OR TITLE:"extracellular vesicle") AND TITLE:"glioma" AND (TITLE:"prognostic" OR TITLE:"signature" OR TITLE:"risk")',
 "exosome_signature_tcga": '(TITLE:"exosome") AND (TITLE:"signature" OR TITLE:"risk model") AND (TITLE:"prognosis" OR TITLE:"prognostic") AND cancer',
 "shap_clinical": '(TITLE:"SHAP" OR ABSTRACT:"Shapley") AND (TITLE:"prognosis" OR TITLE:"survival" OR TITLE:"prediction") AND (TITLE:"cancer" OR TITLE:"tumour" OR TITLE:"tumor" OR TITLE:"glioma")',
 "lasso_cox_model": '(TITLE:"Lasso") AND (TITLE:"Cox") AND (TITLE:"prognostic" OR TITLE:"signature")',
 "cgga_validation": 'ABSTRACT:"CGGA" AND (TITLE:"prognostic model" OR TITLE:"signature") AND glioma',
 "rab27a_cancer": '(TITLE:"RAB27A" OR TITLE:"Rab27a") AND (TITLE:"cancer" OR TITLE:"tumour" OR TITLE:"tumor" OR TITLE:"glioma" OR TITLE:"metastasis")',
 "ezh2_cancer": 'TITLE:"EZH2" AND (TITLE:"glioma" OR TITLE:"glioblastoma")',
 "exosome_brain": '(TITLE:"exosome" OR TITLE:"extracellular vesicles") AND (TITLE:"brain" OR TITLE:"glioblastoma" OR TITLE:"glioma")',
 "gbm_statistics": 'TITLE:"CBTRUS" AND (TITLE:"brain" OR TITLE:"central nervous system")',
 "glioma_nomogram": 'TITLE:"nomogram" AND (TITLE:"glioma" OR TITLE:"glioblastoma")',
}

def search(query, size=12, core=True):
    url = (f"{API}?query={urllib.parse.quote(query + ' AND (SRC:MED)')}&format=json"
           f"&resultType={'core' if core else 'lite'}&pageSize={size}&sort=CITED%20desc")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r).get("resultList", {}).get("result", [])

def fmt(a):
    auth = a.get("authorString", "").rstrip(".")
    parts = [x.strip() for x in auth.split(",") if x.strip()]
    if len(parts) > 3: auth = ", ".join(parts[:3]) + ", et al"
    ji = a.get("journalInfo", {})
    j = ji.get("journal", {}).get("title", "") or ""
    vol = ji.get("volume", ""); pages = a.get("pageInfo", "")
    types = a.get("pubTypeList", {}).get("pubType", [])
    types = types if isinstance(types, list) else [types]
    return {"pmid": a.get("pmid", ""), "doi": a.get("doi", ""), "title": a.get("title", "").strip(),
            "journal": j, "year": a.get("pubYear", ""), "volume": vol, "pages": pages,
            "authors": auth, "types": ",".join(types), "cited": a.get("citedByCount", 0),
            "is_review": any("review" in t.lower() for t in types),
            "ref": f"{auth}. {a.get('title','').strip().rstrip('.')}. {j}"
                   + (f", {a.get('pubYear','')}" if a.get("pubYear") else "")
                   + (f", {vol}" if vol else "") + (f": {pages}" if pages else "")
                   + (f". doi: {a.get('doi')}" if a.get("doi") else "")
                   + (f". PMID: {a.get('pmid')}" if a.get("pmid") else "")}

result = {"landmarks": [], "topics": {}}
for t in LANDMARKS:
    q = f'TITLE:"{t}"'
    try:
        res = search(q, size=2)
        if res:
            it = fmt(res[0]); it["query"] = t
            result["landmarks"].append(it)
            print(f"[L] {it['pmid']} {it['year']} {it['journal'][:32]:32s} {it['title'][:70]}  {'REVIEW' if it['is_review'] else ''}")
        else:
            print(f"[L] 未命中: {t[:60]}")
    except Exception as e:
        print(f"[L] ERR {t[:40]} {e}")
    time.sleep(0.3)

for k, q in TOPICS.items():
    try:
        res = search(q, size=12)
        items = [fmt(x) for x in res if x.get("pmid")]
        seen, uniq = set(), []
        for it in items:
            key = it["title"].lower()[:55]
            if key in seen: continue
            seen.add(key); uniq.append(it)
        result["topics"][k] = uniq
        print(f"\n=== {k} ({len(uniq)}) ===")
        for it in uniq[:8]:
            print(f"  {it['pmid']} | {it['year']} | {'' if not it['is_review'] else 'REVIEW '}{it['journal'][:30]:30s} | {it['title'][:88]}")
    except Exception as e:
        print(f"[{k}] ERR {e}")
    time.sleep(0.3)

json.dump(result, open(os.path.join(OUT, "参考文献_精查.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print("\nlandmarks:", len(result["landmarks"]), "| topics:", sum(len(v) for v in result["topics"].values()))
