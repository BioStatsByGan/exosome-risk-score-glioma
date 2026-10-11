"""审计：跨组织复现 MR —— 全血 eQTLGen（本地 IEU 预处理格式）工具变量 vs FinnGen 胶质瘤
思路：与脑组织 GTEx 结果做方向一致性比对，检验阴性是否由组织特异性造成。
输出：04_results/audit/cross_tissue_mr.txt
"""
import os, json, math

import pandas as pd

ROOT = "/mnt/f/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
AUD = os.path.join(ROOT, "04_results/audit")
MR = os.path.join(ROOT, "04_results/MR")

blood = pd.read_csv(os.path.join(AUD, "eqtlgen_blood_genes.txt"), sep="\t")
blood = blood[blood["pval.exposure"].notna() & (blood["pval.exposure"] < 5e-8)]
print("eQTLGen 显著 cis-eQTL 行:", len(blood), "| 基因:", blood["gene"].nunique())

# 每个基因取最显著 lead SNP
lead = blood.sort_values("pval.exposure").groupby("gene").head(1)
lead = lead[["gene", "SNP", "effect_allele.exposure", "other_allele.exposure",
             "beta.exposure", "se.exposure", "pval.exposure", "eaf.exposure"]]
lead.columns = ["gene", "SNP", "EA", "OA", "beta_exp", "se_exp", "p_exp", "eaf_exp"]
print("lead SNP 数:", len(lead))
print(lead.to_string(index=False))

# 从 FinnGen 主结局按 rsID 抽这些 SNP
want = set(lead["SNP"])
metaf = json.load(open("/mnt/e/R12.xlsx", "rb")) if False else None
N_CASE, N_CTRL = 544, 378749
src = os.path.join(ROOT, "00_data/finngen_R12_C3_GBM_ASTROCYTOMA_EXALLC.gz")
keep = []
for chunk in pd.read_csv(src, sep="\t", chunksize=500000, low_memory=False, dtype={"#chrom": str}):
    hit = chunk[chunk["rsids"].astype(str).isin(want)]
    if len(hit):
        keep.append(hit)
out = pd.concat(keep, ignore_index=True) if keep else pd.DataFrame()
print("FinnGen 命中 SNP:", len(out))
out.to_csv(os.path.join(AUD, "outcome_by_rsid_blood_snps.txt"), sep="\t", index=False)

# Wald ratio
rows = []
for _, r in lead.iterrows():
    o = out[out["rsids"] == r["SNP"]]
    if not len(o):
        rows.append({"gene": r["gene"], "SNP": r["SNP"], "nSNP": 0, "note": "结局缺该 SNP"})
        continue
    o = o.iloc[0]
    ea_o, oa_o = str(o["alt"]).upper(), str(o["ref"]).upper()
    b_out, se_out = float(o["beta"]), float(o["sebeta"])
    # 等位基因方向对齐
    if r["EA"].upper() == ea_o:
        b_exp, se_exp = r["beta_exp"], r["se_exp"]
    elif r["EA"].upper() == oa_o:
        b_exp, se_exp, b_out = -r["beta_exp"], r["se_exp"], -b_out
    else:
        rows.append({"gene": r["gene"], "SNP": r["SNP"], "nSNP": 0, "note": f"等位基因不匹配 {r['EA']}/{r['OA']} vs {ea_o}/{oa_o}"})
        continue
    beta_mr = b_out / b_exp
    se_mr = abs(se_out / b_exp)
    z = beta_mr / se_mr
    from scipy import stats
    p = 2 * stats.norm.sf(abs(z))
    rows.append({"gene": r["gene"], "SNP": r["SNP"], "nSNP": 1,
                 "beta_exp_blood": b_exp, "p_exp": r["p_exp"], "F": (b_exp / se_exp) ** 2,
                 "OR": math.exp(beta_mr), "OR_low": math.exp(beta_mr - 1.96 * se_mr),
                 "OR_high": math.exp(beta_mr + 1.96 * se_mr), "p": p, "note": ""})

res = pd.DataFrame(rows)
res.to_csv(os.path.join(AUD, "cross_tissue_mr_blood.txt"), sep="\t", index=False)
print("\n===== 全血 eQTLGen → FinnGen 胶质瘤（Wald ratio）=====")
print(res.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
sig = res[(res["nSNP"] == 1) & (res["p"] < 0.05)]
print("\n名义显著基因数:", len(sig), "/", (res["nSNP"] == 1).sum())
