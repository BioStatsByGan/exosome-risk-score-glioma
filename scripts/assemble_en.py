"""组装英文手稿：拼接四部分 + 英文三线表 + 引用编号（Vancouver 风格文献清单）。"""
import csv, json, os, re, urllib.parse, urllib.request

ROOT = "/mnt/f/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
EN = os.path.join(ROOT, "稿", "en")
OUT = os.path.join(ROOT, "稿")

parts = ["part1_abstract_intro.md", "part2_methods_results.md", "part3_discussion.md", "part4_legends.md"]
text = "\n\n".join(open(os.path.join(EN, p), encoding="utf-8").read().rstrip() for p in parts)

## ---------- 英文表格 ----------
def read_csv(name):
    with open(os.path.join(ROOT, "稿", "tables", name), encoding="utf-8") as f:
        return list(csv.reader(f))

def md_table(header, rows, title, note):
    out = [f"**{title}**", "", "| " + " | ".join(header) + " |",
           "|" + "|".join(["---"] * len(header)) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(x) for x in r) + " |")
    out += ["", f"Note. {note}", ""]
    return "\n".join(out)

MAP = {"Female, n (%)": "Female, n (%)", "Male, n (%)": "Male, n (%)",
       "G2, n (%)": "Grade II, n (%)", "G3, n (%)": "Grade III, n (%)", "G4, n (%)": "Grade IV, n (%)",
       "WHO II, n (%)": "WHO grade II, n (%)", "WHO III, n (%)": "WHO grade III, n (%)",
       "WHO IV, n (%)": "WHO grade IV, n (%)",
       "Mutant, n (%)": "IDH mutant, n (%)", "WT, n (%)": "IDH wildtype, n (%)",
       "Wildtype, n (%)": "IDH wildtype, n (%)",
       "Methylated, n (%)": "MGMT methylated, n (%)", "Unmethylated, n (%)": "MGMT unmethylated, n (%)",
       "methylated, n (%)": "MGMT methylated, n (%)", "un-methylated, n (%)": "MGMT unmethylated, n (%)",
       "No. of patients": "No. of patients", "Age, median (IQR)": "Age, median (IQR), years",
       "Dead, n (%)": "Death, n (%)",
       "Follow-up, median (IQR), years": "Follow-up, median (IQR), years"}
a, b = read_csv("表1a_基线_TCGA.csv"), read_csv("表1b_基线_CGGA325.csv")
da = {r[0]: r for r in a[1:]}; db = {r[0]: r for r in b[1:]}
order = [r[0] for r in a[1:]] + [r[0] for r in b[1:] if r[0] not in [x[0] for x in a[1:]]]
rows = [[MAP.get(k, k), *(da[k][1:] if k in da else ["—"] * 3), *(db[k][1:] if k in db else ["—"] * 3)] for k in order]
t1 = md_table(["Characteristic", "TCGA high risk", "TCGA low risk", "P", "CGGA-325 high risk", "CGGA-325 low risk", "P"],
              rows, "Table 1. Baseline clinicopathological characteristics of the TCGA and CGGA-325 cohorts by risk group",
              "Continuous variables are shown as median (interquartile range) and compared by Wilcoxon rank sum tests; categorical variables are shown as counts (percentages) and compared by chi-square tests. MGMT denotes MGMT promoter methylation status. The TCGA cohort was divided into 399 high-risk and 400 low-risk patients and the CGGA-325 cohort into 162 high-risk and 163 low-risk patients by the median risk score.")

t2r = read_csv("表2_多因素Cox_15基因.csv")
t2 = md_table(t2r[0], t2r[1:], "Table 2. Fifteen-gene multivariable Cox prognostic model",
              "Coefficient is the regression coefficient of the multivariable Cox model, HR is the hazard ratio, 95% CI is the 95% confidence interval of the hazard ratio and P is from the Wald test. The C-index of the model in the TCGA training cohort was 0.824.")

t3r = read_csv("表3_独立预后分析.csv")
MOD = {"TCGA univariable": "TCGA, univariable", "TCGA multivariable": "TCGA, multivariable",
       "CGGA-325 multivariable (adjusted for WHO grade)": "CGGA-325, multivariable (adjusted for WHO grade)"}
VAR = {"riskScore": "Risk score", "age": "Age", "Sex": "Sex", "SexMale": "Male (vs female)",
       "Grade": "WHO grade", "GradeG3": "WHO grade G3 (vs G2)", "GradeG4": "WHO grade G4 (vs G2)",
       "IDH": "IDH status", "IDHWT": "IDH wildtype (vs mutant)",
       "GradeWHO II": "WHO grade II", "GradeWHO III": "WHO grade III", "GradeWHO IV": "WHO grade IV"}
t3rows = [[MOD.get(r[0], r[0]), VAR.get(r[1], r[1])] + r[2:] for r in t3r[1:]]
t3 = md_table(["Model", "Variable", "HR", "95% CI", "P"], t3rows,
              "Table 3. Univariable and multivariable independent prognostic analyses",
              "HR is the hazard ratio and 95% CI the 95% confidence interval. Covariates in the TCGA multivariable model were the risk score, age, sex, WHO grade and IDH status. WHO grade and IDH status were highly collinear in the CGGA-325 cohort, so WHO grade alone was used as the covariate.")

text = text.replace("<!--TABLES-->", t1 + "\n" + t2 + "\n" + t3)

## ---------- 引用编号 ----------
pool = {}
for f in ["参考文献_原始.json", "参考文献_精查.json", "参考文献_扩充.json"]:
    j = json.load(open(os.path.join(OUT, f), encoding="utf-8"))
    items = j if isinstance(j, list) else (j.get("landmarks", []) + sum(j.get("topics", {}).values(), []))
    for it in items:
        if it.get("pmid"):
            pool.setdefault(str(it["pmid"]), it)

order_ref = []
for p in re.findall(r"\[REF:(\d+)\]", text):
    if p not in order_ref:
        order_ref.append(p)

missing = [p for p in order_ref if p not in pool]
for p in missing:
    try:
        url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:{p}%20AND%20(SRC:MED)&format=json&resultType=core"
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=40) as r:
            res = json.load(r).get("resultList", {}).get("result", [])
        if res:
            x = res[0]
            auth = x.get("authorString", "").rstrip(".")
            ps = [y.strip() for y in auth.split(",") if y.strip()]
            if len(ps) > 3: auth = ", ".join(ps[:3]) + ", et al"
            ji = x.get("journalInfo", {})
            pool[p] = {"pmid": p, "doi": x.get("doi", ""), "title": x.get("title", "").strip(),
                       "journal": ji.get("journal", {}).get("title", ""), "year": x.get("pubYear", ""),
                       "volume": ji.get("volume", ""), "pages": x.get("pageInfo", ""), "authors": auth}
            print("fetched", p)
    except Exception as e:
        print("fail", p, e)

def vancouver(it):
    j = it.get("journal", "") or ""
    y = it.get("year", ""); v = it.get("volume", ""); pg = it.get("pages", "")
    loc = f"{y}" + (f";{v}" if v else "") + (f":{pg}" if pg else "")
    s = f"{it.get('authors','')}. {it.get('title','').strip().rstrip('.')}. {j}. {loc}."
    if it.get("doi"): s += f" doi:{it['doi']}."
    return s

num = {p: i + 1 for i, p in enumerate(order_ref)}
text = re.sub(r"\[REF:(\d+)\]", lambda m: f"[{num[m.group(1)]}]" if m.group(1) in num else m.group(0), text)
refmd = "# References\n\n" + "\n".join(f"{i}. {vancouver(pool[p])}" for p, i in sorted(num.items(), key=lambda kv: kv[1])) + "\n"
text = text.replace("<!--REFS-->", refmd)

final = os.path.join(OUT, "手稿_英文_v1.md")
open(final, "w", encoding="utf-8").write(text)
open(os.path.join(OUT, "参考文献_英文清单.md"), "w", encoding="utf-8").write(refmd)
print("references:", len(num), "| unresolved:", len(re.findall(r"\[REF:", text)))
print("saved:", final, len(text), "chars")
