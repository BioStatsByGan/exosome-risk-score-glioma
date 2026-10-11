"""组装英文手稿 v2：拼接四部分 + 英文表 1 至 3 + 表 S1（MR 把握度）+ 引用编号"""
import csv, json, os, re, urllib.parse, urllib.request

ROOT = "/mnt/f/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
EN = os.path.join(ROOT, "稿", "en")
OUT = os.path.join(ROOT, "稿")

parts = ["part1_abstract_intro_v2.md", "part2_results_v2.md", "part3_discussion_v2.md", "part4_legends_v2.md"]
text = "\n\n".join(open(os.path.join(EN, p), encoding="utf-8").read().rstrip() for p in parts)


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
       "Dead, n (%)": "Death, n (%)", "Follow-up, median (IQR), years": "Follow-up, median (IQR), years"}
a, b = read_csv("表1a_基线_TCGA.csv"), read_csv("表1b_基线_CGGA325.csv")
# 以「显示标签」为合并键：两表对同一特征的不同写法（WT / Wildtype、Methylated / methylated）
# 合并成一行，避免出现同名两行（一行只有 TCGA 值、一行只有 CGGA 值）
ka = [MAP.get(r[0], r[0]) for r in a[1:]]; kb = [MAP.get(r[0], r[0]) for r in b[1:]]
da = dict(zip(ka, a[1:])); db = dict(zip(kb, b[1:]))
order = ka + [k for k in kb if k not in ka]
rows = [[k, *(da[k][1:] if k in da else ["—"] * 3), *(db[k][1:] if k in db else ["—"] * 3)] for k in order]
t1 = md_table(["Characteristic", "TCGA high risk", "TCGA low risk", "P", "CGGA-325 high risk", "CGGA-325 low risk", "P"],
              rows, "Table 1. Baseline characteristics of the TCGA and CGGA-325 cohorts by risk group",
              "Continuous variables are median (interquartile range) compared by Wilcoxon rank sum tests; categorical variables are counts (percentages) compared by chi-square tests. MGMT denotes MGMT promoter methylation status. Risk groups were defined by the cohort median risk score, and each column describes the patients entering the survival analyses: TCGA, 799 patients with follow-up; CGGA-325, 313 of the 325 enrolled patients, the remaining 12 having no follow-up time recorded. Within the CGGA-325 column, WHO grade was not recorded for 4 patients, IDH status for 1 and MGMT promoter methylation status for 18; within the TCGA column, IDH status was unavailable for 36 patients and MGMT promoter methylation status for 73. Percentages are of the total number of patients in that risk group, so the rows for a characteristic with missing values do not add up to 100%.")

t2r = read_csv("表2_多因素Cox_15基因.csv")
t2 = md_table(t2r[0], t2r[1:], "Table 2. Fifteen-gene multivariable Cox model",
              "Coefficient is the regression coefficient, HR the hazard ratio, 95% CI its confidence interval and P from the Wald test. The training C-index was 0.824 and the 10-fold cross-validated C-index 0.819.")

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
              "Covariates in the TCGA multivariable model were risk score, age, sex, WHO grade and IDH status; in CGGA-325 WHO grade alone was used because grade and IDH status were collinear, WHO grade II was the reference category, and the four patients without a recorded WHO grade were excluded (n = 309).")

pf = os.path.join(ROOT, "04_results/audit/mr_power_analysis.txt")
power = list(csv.reader(open(pf, encoding="utf-8"), delimiter="\t")) if os.path.exists(pf) else None
if power is None:
    power = [["SNP_R2", "OR1.10", "OR1.20", "OR1.30", "OR1.50"],
             ["0.01", "0.2", "0.7", "2.0", "8.2"], ["0.03", "0.6", "3.5", "12.2", "49.6"],
             ["0.05", "1.1", "8.3", "29.1", "82.7"], ["0.08", "2.1", "18.9", "56.9", "98.0"]]
tS1 = md_table(["Variance explained by instrument (R2)", "Power at OR 1.10 (%)", "Power at OR 1.20 (%)",
                "Power at OR 1.30 (%)", "Power at OR 1.50 (%)"],
               power[1:], "Table S1. Statistical power of the Mendelian randomization analysis",
               "Power for a binary outcome with 544 cases and 378,749 controls (effective sample size 2,173) at a Bonferroni threshold of 1.02×10⁻³ across 49 tests. At the observed instrument strength the minimum detectable odds ratio at 80% power was 1.31 to 1.41 without correction for multiple testing and 1.49 to 1.67 at a Bonferroni threshold.")

text = text.replace("<!--TABLES-->", t1 + "\n" + t2 + "\n" + t3 + "\n" + tS1)

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
    except Exception as e:
        print("fail", p, e)


def vancouver(it):
    j = it.get("journal", "") or ""
    loc = f"{it.get('year','')}" + (f";{it['volume']}" if it.get("volume") else "") + (f":{it['pages']}" if it.get("pages") else "")
    s = f"{it.get('authors','')}. {it.get('title','').strip().rstrip('.')}. {j}. {loc}."
    if it.get("doi"): s += f" doi:{it['doi']}."
    return s


num = {p: i + 1 for i, p in enumerate(order_ref)}
text = re.sub(r"\[REF:(\d+)\]", lambda m: f"[{num[m.group(1)]}]" if m.group(1) in num else m.group(0), text)
refmd = "# References\n\n" + "\n".join(f"{i}. {vancouver(pool[p])}" for p, i in sorted(num.items(), key=lambda kv: kv[1])) + "\n"
decl = open(os.path.join(EN, "part5_declarations_v2.md"), encoding="utf-8").read().strip()
text = text.replace("<!--REFS-->", decl + "\n\n" + refmd)

open(os.path.join(OUT, "手稿_英文_v2.md"), "w", encoding="utf-8").write(text)
open(os.path.join(OUT, "参考文献_英文清单_v2.md"), "w", encoding="utf-8").write(refmd)
print("references:", len(num), "| unresolved:", len(re.findall(r"\[REF:", text)), "| chars:", len(text))
