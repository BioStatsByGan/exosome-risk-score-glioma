## 步骤 27：CIBERSORT 免疫浸润（LM22，perm=1000）
## 输入：02_inputs/TCGA_glioma_tumor_log2tpm.csv（log2 TPM，回归到线性值后送入）、00_data/ref.txt（LM22）
## 输出：04_results/27.CIBERSORT/CIBERSORT-Results.txt
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(limma); library(e1071); library(parallel)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
dir.create("04_results/27.CIBERSORT", recursive = TRUE, showWarnings = FALSE)
setwd("04_results/27.CIBERSORT")

## 载入 CIBERSORT 函数（课程原始脚本内含 CoreAlg / doPerm / CIBERSORT）
source("../../01_模板_原始脚本/27.CIBERSORT/tcgaSHAP27.CIBERSORT.R", local = FALSE)
if (!exists("CIBERSORT")) stop("未找到 CIBERSORT 函数")

risk <- read.table("../../04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t",
                   check.names = FALSE, row.names = 1)
ex <- fread("../../02_inputs/TCGA_glioma_tumor_log2tpm.csv", showProgress = FALSE)
gn <- ex[[1]]; mat <- as.matrix(ex[, -1]); rownames(mat) <- gn; rm(ex)
mat <- mat[, intersect(colnames(mat), rownames(risk)), drop = FALSE]
## CIBERSORT 要求线性表达值，log2(TPM+1) 反变换
mat <- 2^mat - 1
mat <- avereps(mat)
mat <- mat[rowMeans(mat) > 0, ]
cat("送入 CIBERSORT:", nrow(mat), "基因 x", ncol(mat), "样本\n")
out <- rbind(ID = colnames(mat), mat)
write.table(out, file = "uniq.symbol.txt", sep = "\t", quote = FALSE, col.names = FALSE)

res <- CIBERSORT("../../00_data/ref.txt", "uniq.symbol.txt", perm = 1000, QN = FALSE)
write.table(res, file = "CIBERSORT-Results.txt", sep = "\t", quote = FALSE, col.names = NA)
cat("CIBERSORT 完成，样本数:", nrow(res), "\n")
cat("P-value<0.05 的样本:", sum(res[, "P-value"] < 0.05), "\n")
cat("STEP27 DONE\n")
