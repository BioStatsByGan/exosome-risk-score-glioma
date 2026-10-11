## 步骤 34：oncoPredict 药敏预测（GDSC2 训练集）
## 输入：02_inputs/TCGA_glioma_tumor_log2tpm.csv、00_data/gdsc/GDSC2_Expr.rds、GDSC2_Res.rds
## 输出：04_results/34.oncoPredict/calcPhenotype_Output/DrugPredictions.csv
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(limma); library(oncoPredict); library(parallel)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
OUT <- "04_results/34.oncoPredict"; dir.create(file.path(OUT, "calcPhenotype_Output"), recursive = TRUE, showWarnings = FALSE)
setwd(OUT)

risk <- read.table("../../04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t",
                   check.names = FALSE, row.names = 1)
ex <- fread("../../02_inputs/TCGA_glioma_tumor_log2tpm.csv", showProgress = FALSE)
gn <- ex[[1]]; data <- as.matrix(ex[, -1]); rownames(data) <- gn; rm(ex)
data <- data[, intersect(colnames(data), rownames(risk)), drop = FALSE]
data <- avereps(data)
data <- data[rowMeans(data) > 0.5, ]
cat("测试集:", nrow(data), "基因 x", ncol(data), "样本\n")

GDSC2_Expr <- readRDS("../../00_data/gdsc/GDSC2_Expr.rds")
GDSC2_Res <- readRDS("../../00_data/gdsc/GDSC2_Res.rds")
GDSC2_Res <- exp(GDSC2_Res)
cat("训练集:", nrow(GDSC2_Expr), "基因 x", ncol(GDSC2_Expr), "细胞系 | 药物数:", ncol(GDSC2_Res), "\n")

calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = data,
              batchCorrect = "eb",
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10,
              printOutput = FALSE,
              removeLowVaringGenesFrom = "rawData")

f <- list.files(".", pattern = "DrugPredictions.*\\.csv$", full.names = TRUE)
cat("输出文件:", f, "\n")
if (length(f)) file.rename(f[1], file.path("calcPhenotype_Output", "DrugPredictions.csv"))
cat("STEP34 DONE\n")
