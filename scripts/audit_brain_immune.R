## 审计补救 B1：脑特异免疫细胞集（含小胶质）的 ssGSEA 敏感性分析
## 目的：LM22 不含小胶质，用含小胶质的基因集复核风险分层的髓系信号是否仍在
## 输出：04_results/audit/brain_immune_ssgsea_sensitivity.txt（分级、分级+IDH 两种校正）
##       04_results/audit/brain_aware_neutrophil_scores.txt（患者级中性粒署名评分，供 39 脚本/图S5 用）
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(GSVA); library(limma)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)

SETS <- list(
  Microglia = c("P2RY12", "TMEM119", "CX3CR1", "SALL1", "GPR34", "SIGLEC8", "OLFML3", "BHLHE41"),
  `Macrophage M1` = c("CD86", "NOS2", "IL1B", "TNF", "CXCL9", "CXCL10", "SOCS3", "IDO1"),
  `Macrophage M2` = c("CD163", "MRC1", "MSR1", "IL10", "TGFB1", "CCL22", "F13A1", "STAB1"),
  Monocyte = c("S100A8", "S100A9", "FCN1", "VCAN", "CD14", "LYZ", "FCGR3A", "CSF1R"),
  `T cell` = c("CD3D", "CD3E", "CD2", "IL7R", "CD8A", "CD4", "CCL5", "TRAC"),
  `B cell` = c("MS4A1", "CD79A", "CD79B", "CD19", "TCL1A", "BANK1"),
  Neutrophil = c("FCGR3B", "CSF3R", "CXCR2", "S100A12", "ELANE", "MPO"),
  Endothelial = c("PECAM1", "VWF", "CLDN5", "FLT1", "CDH5", "EGFL7"))

risk <- read.table("04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t",
                   check.names = FALSE, row.names = 1)
ex <- fread("02_inputs/TCGA_glioma_tumor_log2tpm.csv", showProgress = FALSE)
gn <- ex[[1]]; mat <- as.matrix(ex[, -1]); rownames(mat) <- gn; rm(ex)
mat <- mat[, intersect(colnames(mat), rownames(risk)), drop = FALSE]
risk <- risk[colnames(mat), ]
cat("矩阵:", nrow(mat), "x", ncol(mat), "\n")

for (k in names(SETS)) cat(sprintf("%-16s 命中 %d/%d\n", k, sum(SETS[[k]] %in% rownames(mat)), length(SETS[[k]])))

ss <- gsva(ssgseaParam(mat, SETS, normalize = TRUE), verbose = FALSE)
cli <- read.table("02_inputs/clinical_glioma.txt", header = TRUE, sep = "\t", check.names = FALSE)
cli <- cli[match(colnames(ss), cli$patient), ]
d <- data.frame(t(ss), check.names = FALSE)
d$Risk <- risk$Risk; d$Grade <- factor(cli$grade); d$IDH <- factor(cli$idh)

## 同时给出「仅校正分级」与「校正分级+IDH」两种口径，避免同一数值被两处按不同口径引用
res <- do.call(rbind, lapply(names(SETS), function(k) {
  m1 <- summary(lm(as.formula(sprintf("`%s` ~ Risk", k)), data = d))$coefficients
  m2 <- summary(lm(as.formula(sprintf("`%s` ~ Risk + Grade", k)), data = d))$coefficients
  m3 <- summary(lm(as.formula(sprintf("`%s` ~ Risk + Grade + IDH", k)), data = d))$coefficients
  data.frame(CellSet = k, mean_high = mean(d[[k]][d$Risk == "high"]), mean_low = mean(d[[k]][d$Risk == "low"]),
             p_unadj = m1["Risklow", 4], p_adj_grade = m2["Risklow", 4], p_adj_grade_idh = m3["Risklow", 4])
}))
res <- res[order(res$p_unadj), ]
print(res, row.names = FALSE, digits = 3)
fwrite(res, "04_results/audit/brain_immune_ssgsea_sensitivity.txt", sep = "\t", quote = FALSE)
## 患者级中性粒脑特异署名评分（图 S5 第 6 个估计用）
fwrite(data.frame(patient = colnames(ss), BrainAware_Neutrophil = as.numeric(ss["Neutrophil", ])),
       "04_results/audit/brain_aware_neutrophil_scores.txt", sep = "\t", quote = FALSE)
cat("\n判读：小胶质（LM22 缺失的那一类）在高低危组间的差异（校正分级后）:", 
    ifelse(res$p_adj_grade[res$CellSet == "Microglia"] < 0.05, "仍在", "不显著"), "\n")
cat("AUDIT-B1 DONE\n")
