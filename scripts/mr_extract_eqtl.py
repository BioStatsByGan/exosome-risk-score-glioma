"""MR 层第 1 步：从本地 GTEx v8 脑组织 signif_pairs 中抽取 15 个模型基因的顺式 eQTL
注意：文件名为 .gz 但实为纯文本，用普通文本流读取。
输出：04_results/MR/eqtl_brain_<组织>.txt（并在 eqtl_brain_summary.txt 汇总）
"""
import glob, os, time

D = "/mnt/i/EQTL/eqtldata"
OUT = "/mnt/f/0000_Bioinfo/222_tcgaSHAP_胶质瘤/04_results/MR"
os.makedirs(OUT, exist_ok=True)
GENES = {"PPIA", "UTP20", "RAB27A", "MSN", "EZH2", "GNG12", "ADGRE5", "LRRK2",
         "AQP1", "ACY3", "MYC", "EZR", "CCT3", "SNX18", "EIF2S1"}
summary = []
files = sorted(glob.glob(os.path.join(D, "Brain_*.signif_pairs.txt.gz")))
print("脑组织文件:", len(files))
for f in files:
    tissue = os.path.basename(f).replace(".v8.EUR.signif_pairs.txt.gz", "")
    t0 = time.time()
    header, rows = None, []
    with open(f, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if header is None:
                header = line.rstrip("\n")
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2 and parts[1] in GENES:
                rows.append(line.rstrip("\n"))
    out = os.path.join(OUT, f"eqtl_brain_{tissue}.txt")
    with open(out, "w", encoding="utf-8") as o:
        o.write(header + "\n")
        o.write("\n".join(rows) + ("\n" if rows else ""))
    genes_hit = sorted({r.split("\t")[1] for r in rows})
    summary.append((tissue, len(rows), len(genes_hit), ",".join(genes_hit)))
    print(f"{tissue:36s} {len(rows):6d} SNP  命中基因 {len(genes_hit):2d}  ({time.time()-t0:.0f}s)")

with open(os.path.join(OUT, "eqtl_brain_summary.txt"), "w", encoding="utf-8") as o:
    o.write("tissue\tn_snp\tn_gene\tgenes\n")
    for t, n, ng, g in summary:
        o.write(f"{t}\t{n}\t{ng}\t{g}\n")
tot_genes = sorted({g for _, _, _, gs in summary for g in gs.split(",") if g})
print("\n有顺式 eQTL 的模型基因:", len(tot_genes), tot_genes)
print("无顺式 eQTL 的模型基因:", sorted(GENES - set(tot_genes)))
print("DONE")
