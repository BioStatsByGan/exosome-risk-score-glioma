"""把 CSV 结果转成三线表形式的 Markdown 管道表，写入正文的表章节（表1 至表3）。"""
import csv, os, re

ROOT = "/mnt/f/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
T = os.path.join(ROOT, "稿", "tables")
MD = os.path.join(ROOT, "稿", "手稿_中文_v1.md")

def read_csv(name):
    with open(os.path.join(T, name), encoding="utf-8") as f:
        return list(csv.reader(f))

def md_table(header, rows, title, note):
    out = [f"**{title}**", "", "| " + " | ".join(header) + " |",
           "|" + "|".join(["---"] * len(header)) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(x) if str(x) != "NA" else "NA" for x in r) + " |")
    out += ["", f"注　{note}", ""]
    return "\n".join(out)

# 表1：TCGA 与 CGGA-325 基线并排
a = read_csv("表1a_基线_TCGA.csv"); b = read_csv("表1b_基线_CGGA325.csv")
# CGGA 的性别/分级/IDH/MGMT 标签本地化
zh = {"Female, n (%)": "女性，例数(%)", "Male, n (%)": "男性，例数(%)",
      "WHO II, n (%)": "WHO II 级，例数(%)", "WHO III, n (%)": "WHO III 级，例数(%)",
      "WHO IV, n (%)": "WHO IV 级，例数(%)", "Mutant, n (%)": "突变型，例数(%)",
      "WT, n (%)": "野生型，例数(%)", "Wildtype, n (%)": "野生型，例数(%)",
      "methylated, n (%)": "甲基化，例数(%)", "un-methylated, n (%)": "未甲基化，例数(%)",
      "Methylated, n (%)": "甲基化，例数(%)", "Unmethylated, n (%)": "未甲基化，例数(%)",
      "G2, n (%)": "G2 级，例数(%)", "G3, n (%)": "G3 级，例数(%)", "G4, n (%)": "G4 级，例数(%)",
      "No. of patients": "例数", "Age, median (IQR)": "年龄，中位数(四分位间距)",
      "Dead, n (%)": "死亡，例数(%)", "Follow-up, median (IQR), years": "随访时间，中位数(四分位间距)，年"}
rows = []
da = {r[0]: r for r in a[1:]}; db = {r[0]: r for r in b[1:]}
order = [r[0] for r in a[1:]] + [k for k in [r[0] for r in b[1:]] if k not in [r[0] for r in a[1:]]]
for k in order:
    ra = da.get(k); rb = db.get(k)
    label = zh.get(k, k)
    rows.append([label,
                 ra[1] if ra else "—", ra[2] if ra else "—", ra[3] if ra else "—",
                 rb[1] if rb else "—", rb[2] if rb else "—", rb[3] if rb else "—"])
t1 = md_table(["临床特征", "TCGA 高危组", "TCGA 低危组", "P 值", "CGGA-325 高危组", "CGGA-325 低危组", "P 值"],
              rows, "表1　TCGA 与 CGGA-325 队列按风险分组的基线临床特征",
              "连续变量以中位数(四分位间距)表示，用 Wilcoxon 秩和检验；分类变量以例数(百分比)表示，用卡方检验。MGMT 指 MGMT 启动子甲基化状态。TCGA 队列按风险评分中位数分为高危组 399 例与低危组 400 例，CGGA-325 队列分为高危组 162 例与低危组 163 例。")

# 表2：15 基因多因素 Cox
t2rows = read_csv("表2_多因素Cox_15基因.csv")
t2 = md_table(t2rows[0], t2rows[1:], "表2　15 基因多因素 Cox 预后模型",
              "系数为多因素 Cox 模型的回归系数，HR 为风险比，95%CI 为风险比的 95% 置信区间，P 值由 Wald 检验给出。模型在 TCGA 训练队列的 C 指数为 0.824。")

# 表3：独立预后分析
t3rows = read_csv("表3_独立预后分析.csv")
mod = {"TCGA univariable": "TCGA 单因素", "TCGA multivariable": "TCGA 多因素",
       "CGGA-325 multivariable (adjusted for WHO grade)": "CGGA-325 多因素(校正 WHO 分级)"}
vn = {"riskScore": "风险评分", "age": "年龄", "Sex": "性别", "SexMale": "男性",
      "Grade": "WHO 分级", "GradeG3": "WHO G3 级", "GradeG4": "WHO G4 级",
      "IDH": "IDH 状态", "IDHWT": "IDH 野生型",
      "GradeWHO II": "WHO II 级", "GradeWHO III": "WHO III 级", "GradeWHO IV": "WHO IV 级"}
t3rows = [[mod.get(r[0], r[0]), vn.get(r[1], r[1])] + r[2:] for r in t3rows[1:]]
t3 = md_table(["模型", "变量", "HR", "95%CI", "P 值"], t3rows,
              "表3　单因素与多因素独立预后分析",
              "HR 为风险比，95%CI 为 95% 置信区间。TCGA 多因素模型协变量包括风险评分、年龄、性别、WHO 分级与 IDH 状态。CGGA-325 队列中 WHO 分级与 IDH 状态高度共线，故仅以 WHO 分级作为协变量。")

# 替换正文表章节
s = open(MD, encoding="utf-8").read()
new_section = "## 表\n\n" + t1 + "\n" + t2 + "\n" + t3
s = re.sub(r"## 表\n[\s\S]*$", new_section, s)
open(MD, "w", encoding="utf-8").write(s)
print("表章节已写入正文，长度", len(s))
