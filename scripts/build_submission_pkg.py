"""组装投稿包：正文 + 英文命名的图件 + 表格 + 清单"""
import os, shutil
from docx import Document

ROOT = "/mnt/f/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
PKG = os.path.join(ROOT, "投稿包")
for sub in ["Figures", "Tables"]:
    os.makedirs(os.path.join(PKG, sub), exist_ok=True)

# 1) 正文
shutil.copy(os.path.join(ROOT, "稿/Manuscript_English_v2.docx"),
            os.path.join(PKG, "Manuscript_English.docx"))

# 2) 主图 + 补充图（英文命名）
main = {"图1_差异与富集": "Figure1", "图2_模型构建与SHAP": "Figure2",
        "图3_训练队列风险分层": "Figure3", "图4_外部验证": "Figure4",
        "图5_独立预后与临床": "Figure5", "图6_亚组与通路": "Figure6",
        "图7_免疫景观": "Figure7", "图8_TMB与药敏": "Figure8",
        "图9_单细胞": "Figure9",
        "图S1_SHAP个体解释": "FigureS1", "图S2_TIDE免疫治疗预测": "FigureS2",
        "图S3_MR与共定位": "FigureS3", "图S4_药敏箱线": "FigureS4",
        "图S5_中性粒方法比较": "FigureS5"}
src_dirs = [os.path.join(ROOT, "稿/figures"), os.path.join(ROOT, "稿/figures_v2")]
n_fig = 0
for cn, en in main.items():
    for d in src_dirs:
        p = os.path.join(d, cn + ".pdf")
        if os.path.exists(p):
            shutil.copy(p, os.path.join(PKG, "Figures", en + ".pdf")); n_fig += 1; break
    else:
        print("  [缺图]", cn)

# 3) 表格：从 docx 中提取 4 张表 → CSV，并复制已有主表 CSV
d = Document(os.path.join(PKG, "Manuscript_English.docx"))
import csv
names = ["Table1_baseline", "Table2_multicox", "Table3_independence", "TableS1_MR_power"]
for i, tb in enumerate(d.tables):
    nm = names[i] if i < len(names) else f"Table_extra{i}"
    with open(os.path.join(PKG, "Tables", nm + ".csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        for row in tb.rows:
            w.writerow([c.text.strip() for c in row.cells])

# 补充表 S2-S4
for src, dst in [("稿/tables/表S2_分层验证.csv", "TableS2_subgroup.csv"),
                 ("稿/tables/表S3_Hallmark通路.csv", "TableS3_hallmark.csv"),
                 ("稿/tables/表S4_MCPcounter.csv", "TableS4_MCPcounter.csv"),
                 ("稿/tables/表S5_中性粒六方法比较.csv", "TableS5_neutrophil_methods.csv"),
                 ("稿/tables/表S5b_方法相关矩阵.csv", "TableS5b_method_correlation.csv")]:
    sp = os.path.join(ROOT, src)
    if os.path.exists(sp):
        shutil.copy(sp, os.path.join(PKG, "Tables", dst)); print("  补充表已拷入", dst)
    else:
        print("  [缺补充表]", src)

print(f"图件 {n_fig}/13 张，表 {len(d.tables)} 张")
print("\n投稿包内容：")
for r, ds, fs in os.walk(PKG):
    for fn in sorted(fs):
        p = os.path.join(r, fn)
        print(f"  {os.path.relpath(p, PKG):<48} {os.path.getsize(p)/1024:8.1f} KB")
