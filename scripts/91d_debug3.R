options(stringsAsFactors = FALSE, warn = -1)
suppressMessages(library(data.table))
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
inst <- fread("04_results/MR/instruments_all.txt")
o <- fread("04_results/MR/outcome_GBM_ASTROCYTOMA.txt"); setnames(o, "#chrom", "chrom", skip_absent = TRUE)
cat("inst$chr class:", class(inst$chr), "| o$chrom class:", class(o$chrom), "\n")
cat("inst 前 3 chr/pos:", paste(head(inst$chr, 3), head(inst$pos, 3), collapse = " "), "\n")
cat("out  前 3 chrom/pos:", paste(head(o$chrom, 3), head(o$pos, 3), collapse = " "), "\n")
o[, chrom_i := as.integer(gsub("chr", "", chrom))]; o[, pos_i := as.integer(pos)]
inst[, chr_i := as.integer(chr)]; inst[, pos_i := as.integer(pos)]
cat("as.integer 后交集位点:", nrow(merge(unique(inst[, .(chr_i, pos_i)]), unique(o[, .(chrom_i, pos_i)]),
                                            by.x = c("chr_i", "pos_i"), by.y = c("chrom_i", "pos_i"))), "\n")
g <- "RAB27A"; t <- "Brain_Cortex"
d <- inst[external_gene_name == g & tissue == t][order(p)][1]
cat("\n", g, t, "lead:", d$chr_i, d$pos_i, " p=", d$p, "\n")
cat("该位点在结局文件中:", nrow(o[chrom_i == d$chr_i & pos_i == d$pos_i]), "\n")
## 检查是否坐标系不同：把该基因在 eQTL 里的 pos 与 FinnGen 中 RAB27A 附近的位点比对
cat("该基因在 eQTL 的 pos 范围:", paste(range(inst[external_gene_name == g]$pos), collapse = "-"), "\n")
near <- o[nearest_genes %like% g]
cat("FinnGen 中 nearest_genes 含", g, "的位点数:", nrow(near), "\n")
if (nrow(near)) print(near[1:3, .(chrom, pos, rsids, nearest_genes, pval)])
