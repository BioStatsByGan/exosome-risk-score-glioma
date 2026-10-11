## 审计（MR 把握度 + 单细胞模块评分零分布）
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(Seurat)})
OUT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤/04_results/audit"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cat("========== A. MR 把握度分析（FinnGen C3_GBM_ASTROCYTOMA：544 例 / 378749 对照）==========\n")
ncase <- 544; nctrl <- 378749
neff <- 4 * ncase * nctrl / (ncase + nctrl)
cat(sprintf("有效样本量 n_eff = %.0f\n", neff))
alpha_bonf <- 0.05 / 49
cat(sprintf("Bonferroni 阈值 alpha = %.2e（49 个检验）\n\n", alpha_bonf))
snp_r2 <- c(0.01, 0.03, 0.05, 0.08)
ors <- c(1.10, 1.20, 1.30, 1.50)
tab <- do.call(rbind, lapply(snp_r2, function(r2) {
  data.frame(SNP_R2 = r2,
             OR1.10 = round(100 * pchisq(qchisq(1 - alpha_bonf, 1), 1, ncp = neff * r2 * log(1.10)^2, lower.tail = FALSE), 1),
             OR1.20 = round(100 * pchisq(qchisq(1 - alpha_bonf, 1), 1, ncp = neff * r2 * log(1.20)^2, lower.tail = FALSE), 1),
             OR1.30 = round(100 * pchisq(qchisq(1 - alpha_bonf, 1), 1, ncp = neff * r2 * log(1.30)^2, lower.tail = FALSE), 1),
             OR1.50 = round(100 * pchisq(qchisq(1 - alpha_bonf, 1), 1, ncp = neff * r2 * log(1.50)^2, lower.tail = FALSE), 1))
}))
cat("把握度 (%)，行=工具解释的暴露方差，列=真实 OR：\n")
print(tab, row.names = FALSE)
fwrite(tab, file.path(OUT, "mr_power_analysis.txt"), sep = "\t", quote = FALSE)
## 反解：80% 把握度可检出的最小 OR
for (r2 in snp_r2) {
  f <- function(b) pchisq(qchisq(1 - alpha_bonf, 1), 1, ncp = neff * r2 * b^2, lower.tail = FALSE) - 0.80
  or_min <- tryCatch(exp(uniroot(f, c(0.001, 2))$root), error = function(e) NA)
  cat(sprintf("R2 = %.2f 时，80%% 把握度可检出的最小 OR = %.2f\n", r2, or_min))
}

cat("\n========== B. 单细胞风险模块评分：随机基因集零分布 ==========\n")
obj <- readRDS("F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤/04_results/SC/cgga_scRNA_seurat.rds")
GENES <- c("PPIA", "UTP20", "RAB27A", "MSN", "EZH2", "GNG12", "ADGRE5", "LRRK2",
           "AQP1", "ACY3", "MYC", "EZR", "CCT3", "SNX18", "EIF2S1")
GENES <- intersect(GENES, rownames(obj))
## 缓存对象保存于注释之前，需从已保存的注释表补回 celltype
if (!"celltype" %in% colnames(obj@meta.data)) {
  ann <- fread("F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤/04_results/SC/cell_annotation.txt")
  ct <- setNames(ann$celltype, ann$cell)
  obj$celltype <- factor(unname(ct[colnames(obj)]))
  cat("从注释表补回 celltype，缺失:", sum(is.na(obj$celltype)), "\n")
}
cat("可用模型基因:", length(GENES), "\n")
obj <- AddModuleScore(obj, features = list(GENES), name = "Obs", seed = 1)
obs <- tapply(obj$Obs1, obj$celltype, mean)
ct_ord <- names(obs)

## 随机基因集：从表达量分布相近的基因里抽 15 个，做 200 次
expr_mean <- rowMeans(GetAssayData(obj, layer = "data"))
target <- mean(expr_mean[GENES])
pool <- names(sort(abs(expr_mean - target)))[1:1000]
set.seed(7)
mat <- matrix(NA, 200, length(ct_ord), dimnames = list(NULL, ct_ord))
for (i in 1:200) {
  g <- sample(pool, length(GENES))
  o2 <- AddModuleScore(obj, features = list(g), name = paste0("Rnd", i), seed = i)
  mat[i, ] <- tapply(o2@meta.data[[paste0("Rnd", i, "1")]], o2$celltype, mean)
}
cat("\n观测值（各细胞类型风险模块均值）:\n"); print(round(obs, 3))
cat("\n随机基因集零分布（均值 ± SD）:\n")
print(round(data.frame(mean = colMeans(mat), sd = apply(mat, 2, sd)), 3))
## 观测值在零分布中的分位（按细胞类型）
pct <- sapply(ct_ord, function(c) mean(mat[, c] < obs[c]))
cat("\n观测值在零分布的百分位（越接近 1 说明该细胞类型评分偏高不是随机技术效应）:\n")
print(round(sort(pct, decreasing = TRUE), 3))
fwrite(data.frame(celltype = ct_ord, observed = as.numeric(obs[ct_ord]),
                  null_mean = colMeans(mat), null_sd = apply(mat, 2, sd),
                  percentile = as.numeric(pct[ct_ord])),
       file.path(OUT, "scrna_module_null_test.txt"), sep = "\t", quote = FALSE)
cat("\n随机集的细胞类型排序（前 3 名出现次数）:\n")
top3 <- apply(mat, 1, function(x) paste(names(sort(x, decreasing = TRUE))[1:3], collapse = ">"))
print(head(sort(table(top3), decreasing = TRUE), 5))
cat("AUDIT-SC DONE\n")
