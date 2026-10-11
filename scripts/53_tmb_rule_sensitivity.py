# -*- coding: utf-8 -*-
"""步骤 53：TMB 文件选择规则的敏感性分析（回应外部审查第 9 条）
现有规则：同一患者多份 MAF 时取「变异总数最多」的一份。该规则可能偏向高突变样本。
本脚本：① 报告每个患者的文件数分布与所选文件的样本类型构成；
        ② 用替代规则（最少文件、中位、仅原发 01 样本、跨文件合并去重）分别重算 TMB；
        ③ 每种规则下重做「高危 vs 低危」比较与「校正分级后」的关联，检验结论是否依赖规则。
输出：04_results/audit/tmb_rule_sensitivity.txt
      稿/tables/表S15_TMB文件规则敏感性.csv
"""
import glob, gzip, os, statistics
import pandas as pd
import numpy as np
from scipy import stats

ROOT = r"F:\0000_Bioinfo\222_tcgaSHAP_胶质瘤"
MAF = os.path.join(ROOT, "00_data/maf/extracted")
OUTD = os.path.join(ROOT, "04_results/audit")
os.makedirs(OUTD, exist_ok=True)

NS = {"Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Del",
      "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site", "Translation_Start_Site"}
EXOME_MB = 38.0

files = sorted(glob.glob(os.path.join(MAF, "**", "*.maf.gz"), recursive=True)) + \
        sorted(glob.glob(os.path.join(MAF, "**", "*.maf"), recursive=True))

recs, per_pat = [], {}
for f in files:
    try:
        op = gzip.open if f.endswith(".gz") else open
        with op(f, "rt", errors="replace") as fh:
            df = pd.read_csv(fh, sep="\t", comment="#", low_memory=False)
        if "Tumor_Sample_Barcode" not in df.columns:
            continue
        b = str(df["Tumor_Sample_Barcode"].iloc[0])
        n_ns = int(df["Variant_Classification"].isin(NS).sum()) if "Variant_Classification" in df.columns else 0
        pat = b[:12]
        st = b[13:15]
        KEYC = ("Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Allele")
        if all(c in df.columns for c in KEYC):
            keys = set(zip(df["Chromosome"].astype(str), df["Start_Position"], df["End_Position"],
                           df["Reference_Allele"].astype(str), df["Tumor_Allele"].astype(str)))
        else:
            keys = None
        recs.append({"file": os.path.basename(os.path.dirname(f)), "patient": pat, "sample_type": st,
                     "n_nonsyn": n_ns, "n_total": len(df)})
        per_pat.setdefault(pat, []).append((n_ns, len(df), st, keys))
    except Exception as e:
        print("skip", os.path.basename(f), type(e).__name__)

f_ = open(os.path.join(OUTD, "tmb_rule_sensitivity.txt"), "w", encoding="utf-8")
def say(s=""):
    print(s); f_.write(s + "\n")

recs = pd.DataFrame(recs)
tum = recs[recs["sample_type"].isin(["01", "02"])]
say("== TMB 文件选择规则敏感性（%s）==" % pd.Timestamp.now().strftime("%Y-%m-%d %H:%M"))
say("MAF 文件 %d 个；肿瘤样本（01/02）文件 %d 个；患者 %d 例" % (len(recs), len(tum), tum["patient"].nunique()))
cnt = tum.groupby("patient").size()
say("每个患者的肿瘤 MAF 文件数：1 份 %d 例、2 份 %d 例、≥3 份 %d 例（最多 %d 份）"
    % (sum(cnt == 1), sum(cnt == 2), sum(cnt >= 3), cnt.max()))
say("样本类型构成（全部肿瘤文件）：01（原发）%d，02（复发）%d"
    % (sum(tum["sample_type"] == "01"), sum(tum["sample_type"] == "02")))

## 各种规则
def rule_max(d):
    return d.sort_values("n_total", ascending=False).iloc[0]
def rule_min(d):
    return d.sort_values("n_total", ascending=True).iloc[0]
def rule_med(d):
    return d.sort_values("n_total").iloc[len(d) // 2]
def rule_primary(d):
    dd = d[d["sample_type"] == "01"]
    return (dd if len(dd) else d).sort_values("n_total", ascending=False).iloc[0]

rows = []
sel_type = {}
for name, fn in [("max_variants (current)", rule_max), ("min_variants", rule_min),
                 ("median_file", rule_med), ("primary_sample_only", rule_primary)]:
    sel = tum.groupby("patient", group_keys=False).apply(fn)
    sel_type[name] = sel["sample_type"].value_counts().to_dict()
    rows.append({"rule": name, "n_patients": len(sel),
                 "n_01": int((sel["sample_type"] == "01").sum()),
                 "n_02": int((sel["sample_type"] == "02").sum()),
                 "median_TMB": round(float((sel["n_nonsyn"] / EXOME_MB).median()), 3)})

# 合并去重规则：仅对「有 >1 份文件」的患者需要合并（其余患者与单文件规则相同）
pool = {}
for pat, lst in per_pat.items():
    if len(lst) == 1:
        pool[pat] = lst[0][0]
        continue
    ks = set()
    ok = True
    for _, _, _, keys in lst:
        if keys is None:
            ok = False; break
        ks |= keys
    pool[pat] = len(ks) if ok else max(x[0] for x in lst)
rows.append({"rule": "pooled_unique_variants", "n_patients": len(pool), "n_01": np.nan, "n_02": np.nan,
             "median_TMB": round(float(np.median([v / EXOME_MB for v in pool.values()])), 3)})

say("\n各规则下所选文件类型与 TMB 中位数：")
for r in rows:
    say("  %-24s 患者 %4d | 01 %s / 02 %s | 中位 TMB %.3f"
        % (r["rule"], r["n_patients"], r["n_01"], r["n_02"], r["median_TMB"]))

## 与风险分组、分级的关系（用当前的 risk.TCGA.txt 与临床表）
risk = pd.read_csv(os.path.join(ROOT, "04_results/20.SHAP/risk.TCGA.txt"), sep="\t")
cli = pd.read_csv(os.path.join(ROOT, "02_inputs/clinical_glioma.txt"), sep="\t")
clin = pd.merge(risk[["id", "riskScore", "Risk"]], cli[["patient", "grade"]], left_on="id", right_on="patient")

say("\n各规则下的「高危 vs 低危」TMB 比较（log2(TMB+1)，Wilcoxon；并给出校正 WHO 分级后的线性模型系数）：")
out = []
for r in rows:
    if r["rule"] == "pooled_unique_variants":
        s = pd.Series(pool, name="n_ns").rename_axis("patient")
    else:
        fn = {"max_variants (current)": rule_max, "min_variants": rule_min,
              "median_file": rule_med, "primary_sample_only": rule_primary}[r["rule"]]
        sel = tum.groupby("patient", group_keys=False).apply(fn)
        s = sel.set_index("patient")["n_nonsyn"]
    d = clin.merge(pd.DataFrame({"patient": s.index, "n_ns": s.values}), on="patient")
    d["logTMB"] = np.log2(d["n_ns"] / EXOME_MB + 1)
    hi = d.loc[d["Risk"].str.contains("high", case=False, na=False), "logTMB"]
    lo = d.loc[d["Risk"].str.contains("low", case=False, na=False), "logTMB"]
    p = stats.mannwhitneyu(hi, lo).pvalue if len(hi) and len(lo) else np.nan
    dd = d[d["grade"].isin(["G2", "G3", "G4"]) & d["Risk"].notna()].copy()
    dd["Grade"] = pd.Categorical(dd["grade"], categories=["G2", "G3", "G4"])
    dd["RiskC"] = pd.Categorical(dd["Risk"].str.lower(), categories=["low", "high"])
    beta = pv = np.nan
    if dd["RiskC"].nunique() == 2 and dd["Grade"].nunique() >= 2:
        import statsmodels.formula.api as smf
        m = smf.ols("logTMB ~ C(RiskC) + C(Grade)", data=dd).fit()
        beta = m.params.get("C(RiskC)[T.high]", np.nan)
        pv = m.pvalues.get("C(RiskC)[T.high]", np.nan)
    say("  %-24s n=%d | 中位 log2TMB 高危 %.2f vs 低危 %.2f | P=%.2g | 校正分级后 β(高危)=%.3f，P=%.3f"
        % (r["rule"], len(d), hi.median(), lo.median(), p, beta, pv))
    out.append({"rule": r["rule"], "n_patients": len(d), "median_log2TMB_high": round(float(hi.median()), 3),
                "median_log2TMB_low": round(float(lo.median()), 3), "P_high_vs_low": "%.2g" % p,
                "beta_highrisk_grade_adj": round(float(beta), 3), "P_grade_adj": "%.3f" % pv})

tab = pd.DataFrame(out)
tab.to_csv(os.path.join(ROOT, "稿/tables/表S15_TMB文件规则敏感性.csv"), index=False)
say("\n表 S15 已写入 稿/tables/表S15_TMB文件规则敏感性.csv")
f_.close()
print("written:", os.path.join(OUTD, "tmb_rule_sensitivity.txt"))
