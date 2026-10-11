"""MR 层第 2 步（修订）：按 chr:pos 匹配从 FinnGen R12 抽结局统计量
输出：04_results/MR/outcome_<tag>.txt（含 chrom, pos, rsids, alt, ref, beta, sebeta, pval, af_alt）
"""
import glob, os

import pandas as pd

ROOT = "/mnt/f/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
MR = os.path.join(ROOT, "04_results/MR")
DATA = os.path.join(ROOT, "00_data")

frames = []
for f in sorted(glob.glob(os.path.join(MR, "eqtl_brain_Brain_*.txt"))):
    tissue = os.path.basename(f).replace("eqtl_brain_", "").replace(".txt", "")
    d = pd.read_csv(f, sep="\t")
    d["tissue"] = tissue
    frames.append(d)
eqtl = pd.concat(frames, ignore_index=True)
eqtl.to_csv(os.path.join(MR, "instruments_all.txt"), sep="\t", index=False)
ch = pd.to_numeric(eqtl["chr"], errors="coerce")
ps = pd.to_numeric(eqtl["pos"], errors="coerce")
loci = {(int(c), int(p)) for c, p in zip(ch.dropna(), ps[ch.notna()])}   # 只保留常染色体位点
print("eQTL 行:", len(eqtl), "| 唯一 chr:pos 位点:", len(loci),
      "| 有 rsID 的行:", eqtl["rsids2"].notna().sum())

import sys
TAGS = {"GBM_ASTROCYTOMA": "finngen_R12_C3_GBM_ASTROCYTOMA_EXALLC.gz",
        "BRAIN_EXALLC": "finngen_R12_C3_BRAIN_EXALLC.gz"}
want = sys.argv[1:] or list(TAGS)
for tag in want:
    fn = TAGS[tag]
    src = os.path.join(DATA, fn)
    out = os.path.join(MR, f"outcome_{tag}.txt")
    keep, n = [], 0
    for chunk in pd.read_csv(src, sep="\t", chunksize=500000, low_memory=False, dtype={"#chrom": str}):
        n += len(chunk)
        if n % 2000000 == 0:
            print(f"  {tag}: 已扫描 {n/1e6:.0f}M 行，命中 {sum(len(k) for k in keep)}", flush=True)
        ch = chunk["#chrom"].astype(str).str.replace("chr", "", regex=False)
        ok = pd.to_numeric(ch, errors="coerce")
        pos = pd.to_numeric(chunk["pos"], errors="coerce")
        m = [(int(a), int(b)) in loci for a, b in zip(ok.fillna(-1), pos.fillna(-1))]
        hit = chunk[pd.Series(m, index=chunk.index)]
        if len(hit):
            keep.append(hit)
    res = pd.concat(keep, ignore_index=True) if keep else pd.DataFrame()
    res.to_csv(out, sep="\t", index=False)
    print(f"{tag}: 扫描 {n} 行，命中 {len(res)} 个位点 -> {os.path.basename(out)}")
