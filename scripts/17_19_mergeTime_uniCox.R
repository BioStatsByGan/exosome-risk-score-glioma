## 步骤 17 + 19：TCGA 表达与生存合并（mergeTime）+ 单因素 Cox
## 输入：02_inputs/TCGA_glioma_tumor_log2tpm.csv（基因×样本）、02_inputs/clinical_glioma.txt、
##       04_results/13.diff/glioma_exosome_diff.txt（差异外泌体基因 422 个）
## 输出：02_inputs/TCGA.expTime.txt、04_results/19.uniCox/TCGA.uniCox.txt、TCGA.uniSigExp.txt
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(survival); library(ggplot2)})

ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
setwd(ROOT); dir.create("04_results/19.uniCox", recursive = TRUE, showWarnings = FALSE)

deg <- fread("04_results/13.diff/glioma_exosome_diff.txt")
genes <- unique(deg$gene)
cat("差异外泌体基因:", length(genes), "\n")

ex <- fread("02_inputs/TCGA_glioma_tumor_log2tpm.csv", showProgress = FALSE)
gn <- ex[[1]]; mat <- as.matrix(ex[, -1]); rownames(mat) <- gn; rm(ex)
mat <- mat[rownames(mat) %in% genes, , drop = FALSE]
cat("表达矩阵命中基因:", nrow(mat), " 样本:", ncol(mat), "\n")

cli <- fread("02_inputs/clinical_glioma.txt")
cli <- cli[!is.na(futime) & !is.na(fustat) & futime > 0]
same <- intersect(colnames(mat), cli$patient)
cat("生存信息可匹配样本:", length(same), "\n")
cli <- as.data.frame(cli); rownames(cli) <- cli$patient; cli <- cli[same, ]
mat <- mat[, same, drop = FALSE]

rt <- data.frame(futime = cli$futime, fustat = cli$fustat, t(mat))
out <- cbind(id = rownames(rt), rt)
fwrite(out, "02_inputs/TCGA.expTime.txt", sep = "\t", quote = FALSE)
cat("saved: 02_inputs/TCGA.expTime.txt\n")

## ---- 单因素 Cox ----
rt <- read.table("02_inputs/TCGA.expTime.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
pfilter <- 0.05
outTab <- data.frame(); sigGenes <- c("futime", "fustat")
for (g in colnames(rt)[3:ncol(rt)]) {
  cox <- coxph(Surv(futime, fustat) ~ rt[, g], data = rt)
  s <- summary(cox)
  if (s$coefficients[, "Pr(>|z|)"] < pfilter) {
    sigGenes <- c(sigGenes, g)
    outTab <- rbind(outTab, data.frame(
      id = g, HR = s$conf.int[, "exp(coef)"], HR.95L = s$conf.int[, "lower .95"],
      HR.95H = s$conf.int[, "upper .95"], pvalue = s$coefficients[, "Pr(>|z|)"]))
  }
}
outTab <- outTab[order(outTab$HR, decreasing = TRUE), ]
fwrite(outTab, "04_results/19.uniCox/TCGA.uniCox.txt", sep = "\t", quote = FALSE)
cat("单因素 Cox 显著基因(p<0.05):", nrow(outTab), "\n")
cat("  保护性 HR<1:", sum(outTab$HR < 1), " 风险性 HR>1:", sum(outTab$HR > 1), "\n")
cat("  HR 范围:", sprintf("%.2f - %.2f", min(outTab$HR), max(outTab$HR)), "\n")

uni <- rt[, sigGenes]
fwrite(cbind(id = rownames(uni), uni), "04_results/19.uniCox/TCGA.uniSigExp.txt", sep = "\t", quote = FALSE)
cat("saved: TCGA.uniSigExp.txt\n")
cat("STEP17-19 DONE\n")
