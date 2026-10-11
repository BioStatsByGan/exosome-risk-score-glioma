## 步骤 29：ssGSEA 免疫功能评分（28 免疫细胞 + 免疫相关功能）高低危组比较
## 输入：02_inputs/TCGA_glioma_tumor_log2tpm.csv、04_results/20.SHAP/risk.TCGA.txt、immune.gmt
## 输出：04_results/29.immFunction/ssgsea_scores.txt、ssgsea_diff.txt、03_figures/Fig6A_ssgsea_immune.pdf
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({
  library(data.table); library(GSVA); library(GSEABase); library(ggplot2)
  library(reshape2); library(ggpubr); library(patchwork)
})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
source("scripts/plot_style.R")
dir.create("04_results/29.immFunction", recursive = TRUE, showWarnings = FALSE)

risk <- read.table("04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t",
                   check.names = FALSE, row.names = 1)
ex <- fread("02_inputs/TCGA_glioma_tumor_log2tpm.csv", showProgress = FALSE)
gn <- ex[[1]]; mat <- as.matrix(ex[, -1]); rownames(mat) <- gn; rm(ex)
mat <- mat[, intersect(colnames(mat), rownames(risk)), drop = FALSE]
risk <- risk[colnames(mat), ]
cat("表达矩阵:", nrow(mat), "基因 x", ncol(mat), "样本\n")

cache <- "04_results/29.immFunction/ssgsea_scores.txt"
if (file.exists(cache)) {
  cat("读取已算好的 ssGSEA 缓存\n")
  ss <- as.matrix(read.table(cache, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1))
} else {
  gmt <- getGmt("01_模板_原始脚本/29.immFunction/immune.gmt")
  ss <- gsva(ssgseaParam(mat, gmt, normalize = TRUE), verbose = FALSE)
  fwrite(cbind(id = rownames(ss), as.data.frame(ss)), cache, sep = "\t", quote = FALSE)
}
cat("ssGSEA 结果:", nrow(ss), "个免疫特征 x", ncol(ss), "样本\n")

## 组间比较
grp <- setNames(as.character(risk$Risk), rownames(risk))   # 必须带名字，否则 grp[dl$sample] 全 NA
res <- do.call(rbind, lapply(rownames(ss), function(x) {
  v <- as.numeric(ss[x, ])
  p <- wilcox.test(v[grp == "high"], v[grp == "low"])$p.value
  data.frame(Feature = x, mean_high = mean(v[grp == "high"]), mean_low = mean(v[grp == "low"]),
             diff = mean(v[grp == "high"]) - mean(v[grp == "low"]), p = p)
}))
res$FDR <- p.adjust(res$p, "BH")
res <- res[order(res$p), ]
fwrite(res, "04_results/29.immFunction/ssgsea_diff.txt", sep = "\t", quote = FALSE)
cat("显著(FDR<0.05)免疫特征:", sum(res$FDR < 0.05), "/", nrow(res), "\n")
print(head(res, 8), row.names = FALSE, digits = 3)

## 箱线图（按 p 值排序，取全部 29 个特征，分两页）
dl <- data.frame(sample = rep(colnames(ss), each = nrow(ss)),
                 Feature = rep(rownames(ss), times = ncol(ss)),
                 score = as.vector(ss))
dl$Risk <- grp[dl$sample]
dl$Risk <- factor(dl$Risk, levels = c("low", "high"), labels = c("Low risk", "High risk"))
dl$Feature <- factor(dl$Feature, levels = res$Feature)
lab <- setNames(sprintf("%s\nP = %s", res$Feature,
                        ifelse(res$p < 0.001, "<0.001", sprintf("%.3f", res$p))), res$Feature)

mk <- function(feats) {
  d <- dl[dl$Feature %in% feats, ]
  ggplot(d, aes(Feature, score, fill = Risk)) +
    geom_boxplot(width = 0.66, outlier.size = 0.35, linewidth = 0.3, colour = "grey25") +
    scale_fill_manual(values = c("Low risk" = "#0072B2", "High risk" = "#D55E00"), name = "") +
    scale_x_discrete(labels = lab) +
    labs(x = NULL, y = "ssGSEA score") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          legend.position = "top")
}
half <- ceiling(nrow(res) / 2)
save_pdf(mk(res$Feature[1:half]), "03_figures/Fig6A_ssgsea_immune.pdf", 10, 5.2)
save_pdf(mk(res$Feature[(half + 1):nrow(res)]), "03_figures/Fig6B_ssgsea_immune2.pdf", 10, 5.2)
cat("STEP29 DONE\n")
