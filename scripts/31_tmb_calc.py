"""步骤 31 输入：从 GDC MAF 计算 TMB（每 Mb 非同义突变数，外显子组捕获区按 38 Mb 计）
输出：04_results/31.riskTMB/TMB.txt
"""
import glob, gzip, json, os

import pandas as pd

ROOT = r"F:\0000_Bioinfo\222_tcgaSHAP_胶质瘤"          # Windows 下直接可用
if not os.path.isdir(ROOT):                              # 兼容 WSL 环境
    ROOT = "/mnt/f/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
MAF = os.path.join(ROOT, "00_data/maf/extracted")
OUT = os.path.join(ROOT, "04_results/31.riskTMB")
os.makedirs(OUT, exist_ok=True)

NS = {"Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
      "In_Frame_Del", "In_Frame_Ins", "Splice_Site", "Translation_Start_Site"}
EXOME_MB = 38.0

files = sorted(glob.glob(os.path.join(MAF, "**", "*.maf.gz"), recursive=True)) + \
        sorted(glob.glob(os.path.join(MAF, "**", "*.maf"), recursive=True))
print("MAF 文件数:", len(files))
rows, n_used = [], 0
for f in files:
    try:
        op = gzip.open if f.endswith(".gz") else open
        with op(f, "rt", errors="replace") as fh:
            df = pd.read_csv(fh, sep="\t", comment="#", low_memory=False)
        if "Tumor_Sample_Barcode" not in df.columns or "Variant_Classification" not in df.columns:
            continue
        b = df["Tumor_Sample_Barcode"].astype(str).iloc[0]
        n = int(df["Variant_Classification"].isin(NS).sum())
        rows.append({"Tumor_Sample_Barcode": b, "n_nonsyn": n, "n_total": len(df)})
        n_used += 1
    except Exception as e:
        print("跳过", os.path.basename(f), type(e).__name__)

tmb = pd.DataFrame(rows)
print("可用 MAF:", n_used, "| 样本数:", tmb.shape[0])
tmb["patient"] = tmb["Tumor_Sample_Barcode"].str[:12]
tmb["sample_type"] = tmb["Tumor_Sample_Barcode"].str[13:15]
tmb = tmb[tmb["sample_type"].isin(["01", "02"])]            # 原发与复发肿瘤
tmb = tmb.sort_values("n_total", ascending=False).drop_duplicates("patient")
tmb["TMB"] = tmb["n_nonsyn"] / EXOME_MB
tmb[["patient", "TMB", "n_nonsyn", "n_total"]].to_csv(os.path.join(OUT, "TMB.txt"), sep="\t", index=False)
print(tmb["TMB"].describe().round(2).to_string())
print("中位 TMB:", round(tmb["TMB"].median(), 2), "突变/Mb")
print("saved:", os.path.join(OUT, "TMB.txt"))
