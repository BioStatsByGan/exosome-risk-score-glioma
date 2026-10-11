"""审计用：从本地 eQTLGen 全血 cis-eQTL（IEU 预处理格式，约 24 GB）抽取 15 个模型基因的行，
用于与脑组织 GTEx eQTL 做跨组织复现（winner's curse 检验）。
输出：04_results/audit/eqtlgen_blood_genes.txt
"""
import glob, os

import pandas as pd

SRC = "/mnt/i/05-eqtlgen_IEU_Source_Data/06_cis_colocation_data"
OUT = "/mnt/f/0000_Bioinfo/222_tcgaSHAP_胶质瘤/04_results/audit"
os.makedirs(OUT, exist_ok=True)
GENES = {"PPIA", "UTP20", "RAB27A", "MSN", "EZH2", "GNG12", "ADGRE5", "LRRK2",
         "AQP1", "ACY3", "MYC", "EZR", "CCT3", "SNX18", "EIF2S1"}

files = sorted(glob.glob(os.path.join(SRC, "eqtlgen_ieu_1mbcis_[0-9]*.csv")))
print("待扫描文件:", [os.path.basename(f) for f in files], flush=True)
kept = []
for f in files:
    n = 0
    for chunk in pd.read_csv(f, chunksize=500000, low_memory=False):
        n += len(chunk)
        gcol = "gene" if "gene" in chunk.columns else None
        if gcol is None:
            print("缺少 gene 列:", chunk.columns.tolist()[:6]); break
        g = chunk[gcol].astype(str).str.strip()
        hit = chunk[g.isin(GENES)]
        if len(hit):
            kept.append(hit)
    print(f"{os.path.basename(f)}: 扫描 {n} 行，累计命中 {sum(len(k) for k in kept)}", flush=True)

res = pd.concat(kept, ignore_index=True) if kept else pd.DataFrame()
res.to_csv(os.path.join(OUT, "eqtlgen_blood_genes.txt"), sep="\t", index=False)
print("命中总数:", len(res))
if len(res):
    print(res.groupby("gene").size().to_string())
