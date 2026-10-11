"""修正工具变量文件的染色体编号（课程包的 signif_pairs 里部分基因 chr 列错写成 1，坐标本身正确）
做法：用 Ensembl REST 按 ENSG 查出真实染色体与基因区间，校正 chr 并做区间自检；
      随后按 chr:pos 重新抽取 FinnGen 结局。
输出：04_results/MR/instruments_fixed.txt、outcome_<tag>.txt
"""
import json, os, sys, urllib.request

import pandas as pd

ROOT = "/mnt/f/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
MR = os.path.join(ROOT, "04_results/MR")
DATA = os.path.join(ROOT, "00_data")

inst = pd.read_csv(os.path.join(MR, "instruments_all.txt"), sep="\t")
ensg = sorted(inst["ensembl_gene_id"].dropna().unique())
req = urllib.request.Request("https://rest.ensembl.org/lookup/id?expand=0",
                             data=json.dumps({"ids": ensg}).encode(),
                             headers={"Content-Type": "application/json", "User-Agent": "Mozilla/5.0"})
with urllib.request.urlopen(req, timeout=60) as r:
    anno = json.load(r)

rows = []
for _, x in inst.iterrows():
    a = anno.get(x["ensembl_gene_id"], {})
    true_chr = str(a.get("seq_region_name", x["chr"]))
    start, end = a.get("start", 0), a.get("end", 0)
    inwin = (start - 1_500_000) <= x["pos"] <= (end + 1_500_000)
    rows.append({"ensembl_gene_id": x["ensembl_gene_id"], "external_gene_name": x["external_gene_name"],
                 "ensembl_symbol": a.get("display_name", ""), "A1": x["A1"], "A2": x["A2"],
                 "chr_file": x["chr"], "chr": true_chr, "pos": int(x["pos"]), "eaf": x["eaf"],
                 "beta": x["beta"], "se": x["se"], "p": x["p"], "n": x["n"], "rsids2": x["rsids2"],
                 "tissue": x["tissue"], "in_cis_window": inwin})
fx = pd.DataFrame(rows)
fx.to_csv(os.path.join(MR, "instruments_fixed.txt"), sep="\t", index=False)

print("染色体校正情况：")
chk = fx.groupby(["external_gene_name", "ensembl_gene_id", "chr_file", "chr"]).size().reset_index(name="n")
print(chk.to_string(index=False))
print("\n坐标落在基因 ±1.5 Mb 窗口内的比例: %.1f%%" % (100 * fx["in_cis_window"].mean()))
print("染色体被修正的行数:", int((fx["chr_file"].astype(str) != fx["chr"]).sum()), "/", len(fx))

TAGS = {"GBM_ASTROCYTOMA": "finngen_R12_C3_GBM_ASTROCYTOMA_EXALLC.gz",
        "BRAIN_EXALLC": "finngen_R12_C3_BRAIN_EXALLC.gz"}
want = sys.argv[1:] or list(TAGS)
ch = pd.to_numeric(fx["chr"], errors="coerce")
ps = pd.to_numeric(fx["pos"], errors="coerce")
loci = {(int(c), int(p)) for c, p in zip(ch.dropna(), ps[ch.notna()])}
print("\n校正后唯一 chr:pos 位点数:", len(loci))

for tag in want:
    src = os.path.join(DATA, TAGS[tag])
    out = os.path.join(MR, f"outcome_{tag}.txt")
    keep, n = [], 0
    for chunk in pd.read_csv(src, sep="\t", chunksize=500000, low_memory=False, dtype={"#chrom": str}):
        n += len(chunk)
        c2 = pd.to_numeric(chunk["#chrom"].astype(str).str.replace("chr", "", regex=False), errors="coerce")
        p2 = pd.to_numeric(chunk["pos"], errors="coerce")
        m = pd.Series([(int(a), int(b)) in loci for a, b in zip(c2.fillna(-1), p2.fillna(-1))], index=chunk.index)
        hit = chunk[m]
        if len(hit):
            keep.append(hit)
        if n % 2000000 == 0:
            print(f"  {tag}: 扫描 {n/1e6:.0f}M 行，命中 {sum(len(k) for k in keep)}", flush=True)
    res = pd.concat(keep, ignore_index=True) if keep else pd.DataFrame()
    res.to_csv(out, sep="\t", index=False)
    print(f"{tag}: 命中 {len(res)} 个位点 -> {os.path.basename(out)}")
