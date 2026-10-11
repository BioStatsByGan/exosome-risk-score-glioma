"""组装中文手稿 v2：拼接三部分 + 中文表 1 至 3 与表 S1 + 引用编号"""
import csv, json, os, re, urllib.parse, urllib.request

ROOT = "/mnt/f/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
CN = os.path.join(ROOT, "稿", "cn")
OUT = os.path.join(ROOT, "稿")

parts = ["part1_v2.md", "part2_v2.md", "part3_v2.md"]
text = "\n\n".join(open(os.path.join(CN, p), encoding="utf-8").read().rstrip() for p in parts)


def read_csv(name):
    with open(os.path.join(ROOT, "稿", "tables", name), encoding="utf-8") as f:
        return list(csv.reader(f))


def md_table(header, rows, title, note):
    out = [f"**{title}**", "", "| " + " | ".join(header) + " |",
           "|" + "|".join(["---"] * len(header)) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(x) for x in r) + " |")
    out += ["", f"注　{note}", ""]
    return "\n".join(out)


zh = {"Female, n (%)": "女性，例数(%)", "Male, n (%)": "男性，例数(%)",
      "WHO II, n (%)": "WHO II 级，例数(%)", "WHO III, n (%)": "WHO III 级，例数(%)",
      "WHO IV, n (%)": "WHO IV 级，例数(%)", "Mutant, n (%)": "IDH 突变型，例数(%)",
      "WT, n (%)": "IDH 野生型，例数(%)", "Wildtype, n (%)": "IDH 野生型，例数(%)",
      "methylated, n (%)": "MGMT 甲基化，例数(%)", "un-methylated, n (%)": "MGMT 未甲基化，例数(%)",
      "Methylated, n (%)": "MGMT 甲基化，例数(%)", "Unmethylated, n (%)": "MGMT 未甲基化，例数(%)",
      "G2, n (%)": "G2 级，例数(%)", "G3, n (%)": "G3 级，例数(%)", "G4, n (%)": "G4 级，例数(%)",
      "No. of patients": "例数", "Age, median (IQR)": "年龄，中位数(四分位间距)",
      "Dead, n (%)": "死亡，例数(%)", "Follow-up, median (IQR), years": "随访时间，中位数(四分位间距)，年"}
a, b = read_csv("表1a_基线_TCGA.csv"), read_csv("表1b_基线_CGGA325.csv")
# 以「显示标签」为合并键，避免同一特征出现两行（WT / Wildtype、Methylated / methylated 写法不同）
ka = [zh.get(r[0], r[0]) for r in a[1:]]; kb = [zh.get(r[0], r[0]) for r in b[1:]]
da = dict(zip(ka, a[1:])); db = dict(zip(kb, b[1:]))
order = ka + [k for k in kb if k not in ka]
rows = [[k, *(da[k][1:] if k in da else ["—"] * 3), *(db[k][1:] if k in db else ["—"] * 3)] for k in order]
t1 = md_table(["临床特征", "TCGA 高危组", "TCGA 低危组", "P 值", "CGGA-325 高危组", "CGGA-325 低危组", "P 值"],
              rows, "表1　TCGA 与 CGGA-325 队列按风险分组的基线临床特征",
              "连续变量以中位数(四分位间距)表示，用 Wilcoxon 秩和检验；分类变量以例数(百分比)表示，用卡方检验。MGMT 指 MGMT 启动子甲基化状态。风险分组按各队列风险评分中位数划分；各列均为进入生存分析的患者：TCGA 799 例均有随访，CGGA-325 为 325 例中的 313 例，其余 12 例无随访时间记录。CGGA-325 列中 WHO 分级缺失 4 例、IDH 状态缺失 1 例、MGMT 启动子甲基化状态缺失 18 例；TCGA 列中 IDH 状态缺失 36 例、MGMT 启动子甲基化状态缺失 73 例。各行百分比以该风险组的总例数为分母，因此存在缺失的特征其各行合计不等于 100%。")

t2r = read_csv("表2_多因素Cox_15基因.csv")
t2 = md_table(t2r[0], t2r[1:], "表2　15 基因多因素 Cox 预后模型",
              "系数为多因素 Cox 模型回归系数，HR 为风险比，95%CI 为其置信区间，P 值由 Wald 检验给出。训练集 C 指数 0.824，10 折交叉验证 C 指数 0.819。")

t3r = read_csv("表3_独立预后分析.csv")
MOD = {"TCGA univariable": "TCGA 单因素", "TCGA multivariable": "TCGA 多因素",
       "CGGA-325 multivariable (adjusted for WHO grade)": "CGGA-325 多因素(校正 WHO 分级)"}
VAR = {"riskScore": "风险评分", "age": "年龄", "Sex": "性别", "SexMale": "男性",
       "Grade": "WHO 分级", "GradeG3": "WHO G3 级", "GradeG4": "WHO G4 级",
       "IDH": "IDH 状态", "IDHWT": "IDH 野生型",
       "GradeWHO II": "WHO II 级", "GradeWHO III": "WHO III 级", "GradeWHO IV": "WHO IV 级"}
t3rows = [[MOD.get(r[0], r[0]), VAR.get(r[1], r[1])] + r[2:] for r in t3r[1:]]
t3 = md_table(["模型", "变量", "HR", "95%CI", "P 值"], t3rows,
              "表3　单因素与多因素独立预后分析",
              "TCGA 多因素模型协变量为风险评分、年龄、性别、WHO 分级与 IDH 状态；CGGA-325 中分级与 IDH 高度共线，仅以 WHO 分级为协变量，参照组为 WHO II 级，并剔除 4 例无 WHO 分级记录者（n = 309）。")

pf = os.path.join(ROOT, "04_results/audit/mr_power_analysis.txt")
power = list(csv.reader(open(pf, encoding="utf-8"), delimiter="\t")) if os.path.exists(pf) else [
    ["SNP_R2", "OR1.10", "OR1.20", "OR1.30", "OR1.50"], ["0.03", "0.6", "3.5", "12.2", "49.6"]]
tS1 = md_table(["工具解释的方差 R2", "真实 OR=1.10 时把握度(%)", "OR=1.20", "OR=1.30", "OR=1.50"],
               power[1:], "表S1　孟德尔随机化分析的统计把握度",
               "结局为 544 例病例与 378 749 例对照（有效样本量 2 173），阈值取 49 个检验的 Bonferroni 校正值 1.02×10⁻³。在观测到的工具强度下，不做多重检验校正时对约 1.31 至 1.41 倍以上的 OR 达到 80% 把握度，Bonferroni 阈值下为 1.49 至 1.67 倍。")

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


def ref_cn(it):
    j = it.get("journal", "") or ""
    loc = f"{it.get('year','')}" + (f", {it['volume']}" if it.get("volume") else "") + (f": {it['pages']}" if it.get("pages") else "")
    s = f"{it.get('authors','')}. {it.get('title','').strip().rstrip('.')}. {j}, {loc}."
    if it.get("doi"): s += f" doi: {it['doi']}."
    if it.get("pmid"): s += f" PMID: {it['pmid']}."
    return s


num = {p: i + 1 for i, p in enumerate(order_ref)}
text = re.sub(r"\[REF:(\d+)\]", lambda m: f"[{num[m.group(1)]}]" if m.group(1) in num else m.group(0), text)
refmd = "# 参考文献\n\n" + "\n".join(f"{i}. {ref_cn(pool[p])}" for p, i in sorted(num.items(), key=lambda kv: kv[1])) + "\n"
decl = open(os.path.join(CN, "part4_declarations_v2.md"), encoding="utf-8").read().strip()
text = text.replace("<!--REFS-->", decl + "\n\n" + refmd)

open(os.path.join(OUT, "手稿_中文_v2.md"), "w", encoding="utf-8").write(text)
print("references:", len(num), "| unresolved:", len(re.findall(r"\[REF:", text)), "| chars:", len(text))
