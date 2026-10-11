"""步骤 33：TIDE（TIDEpy 本地版，替代已退役的 TIDE 网页服务器）
输入：02_inputs/TCGA_glioma_tumor_log2tpm.csv（log2 TPM，还原成线性值后交给 TIDEpy 自己做 log2+中心化）
输出：04_results/33.TIDE/TIDE.csv
"""
import os

import pandas as pd
from tidepy.pred import TIDE

ROOT = "/mnt/f/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
OUT = os.path.join(ROOT, "04_results", "33.TIDE")
os.makedirs(OUT, exist_ok=True)

risk = pd.read_csv(os.path.join(ROOT, "04_results/20.SHAP/risk.TCGA.txt"), sep="\t", index_col=0)
expr = pd.read_csv(os.path.join(ROOT, "02_inputs/TCGA_glioma_tumor_log2tpm.csv"), index_col=0)
cols = [c for c in expr.columns if c in risk.index]
expr = expr[cols]
print("样本数:", len(cols), "基因数:", expr.shape[0])

linear = 2 ** expr - 1          # 还原为线性 TPM，TIDEpy 会自行做 log2(x+1) 与中心化
linear = linear.groupby(level=0).mean()
res = TIDE(linear, cancer="Other")
res.index.name = "SampleID"
res.to_csv(os.path.join(OUT, "TIDE.csv"))
print("TIDE 完成，样本数:", res.shape[0])
print(res.head(3).iloc[:, :6].to_string())
print("响应者比例:", round(100 * (1 - res["No benefits"].mean()), 1), "% | 无获益比例:", round(100 * res["No benefits"].mean(), 1), "%")
