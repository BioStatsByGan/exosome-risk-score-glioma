## 步骤 28：免疫细胞浸润（CIBERSORT 22 细胞）高低危组比较
## 输入：04_results/27.CIBERSORT/CIBERSORT-Results.txt、04_results/20.SHAP/risk.TCGA.txt
## 输出：04_results/28.immDiff/imm_cell_diff.txt、03_figures/Fig6C_imm_cells.pdf、Fig6D_imm_corr.pdf
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(ggplot2); library(ggpubr); library(reshape2); library(data.table); library(pheatmap)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT); source("scripts/plot_style.R")
dir.create("04_results/28.immDiff", recursive = TRUE, showWarnings = FALSE)

imm <- read.table("04_results/27.CIBERSORT/CIBERSORT-Results.txt", header = TRUE, sep = "\t",
                  check.names = FALSE, row.names = 1)
cat("CIBERSORT 结果:", nrow(imm), "样本 | 列:", ncol(imm), "\n")
cat("P-value < 0.05 的样本:", sum(imm$`P-value` < 0.05), "（本分析保留全部样本，另报告该比例）\n")
cells <- setdiff(colnames(imm), c("P-value", "Correlation", "RMSE"))
imm <- imm[, c(cells, "P-value", "Correlation", "RMSE")]
# 样本 ID 归一化到患者
rownames(imm) <- gsub("^(.*?)-(.*?)-(.*?)-.*", "\\1-\\2-\\3", rownames(imm))
imm <- imm[!duplicated(rownames(imm)), ]

risk <- read.table("04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t",
                   check.names = FALSE, row.names = 1)
same <- intersect(rownames(imm), rownames(risk))
cat("配对样本:", length(same), "\n")
imm <- imm[same, , drop = FALSE]; risk <- risk[same, ]

res <- do.call(rbind, lapply(cells, function(cc) {
  hi <- as.numeric(imm[risk$Risk == "high", cc]); lo <- as.numeric(imm[risk$Risk == "low", cc])
  data.frame(Cell = cc, mean_high = mean(hi), mean_low = mean(lo),
             p = wilcox.test(hi, lo)$p.value)
}))
res$FDR <- p.adjust(res$p, "BH")
res <- res[order(res$p), ]
fwrite(res, "04_results/28.immDiff/imm_cell_diff.txt", sep = "\t", quote = FALSE)
cat("组间差异细胞类型（FDR<0.05）:", sum(res$FDR < 0.05), "/", nrow(res), "\n")
print(head(res, 8), row.names = FALSE, digits = 3)

d <- data.frame(sample = rep(rownames(imm), length(cells)),
                Cell = rep(cells, each = nrow(imm)),
                Fraction = as.vector(as.matrix(imm[, cells])))
d$Risk <- factor(rep(risk$Risk, length(cells)), levels = c("low", "high"), labels = c("Low risk", "High risk"))
d$Cell <- factor(d$Cell, levels = res$Cell)
lab <- setNames(sprintf("%s\nP = %s", res$Cell,
                        ifelse(res$p < 0.001, "<0.001", sprintf("%.3f", res$p))), res$Cell)
p1 <- ggplot(d, aes(Cell, Fraction, fill = Risk)) +
  geom_boxplot(width = 0.66, outlier.size = 0.25, linewidth = 0.28, colour = "grey25") +
  scale_fill_manual(values = c("Low risk" = "#0072B2", "High risk" = "#D55E00"), name = "") +
  scale_x_discrete(labels = lab) +
  labs(x = NULL, y = "Estimated fraction") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6.2), legend.position = "top")
save_pdf(p1, "03_figures/Fig6C_imm_cells.pdf", 11, 5.6)

## 风险评分与显著细胞类型的相关性热图
sig <- res$Cell[res$FDR < 0.05]
cm <- cor(as.matrix(imm[, sig]), risk$riskScore, method = "spearman")
rownames(cm) <- sig
save_pdf(NULL, "03_figures/tmp.pdf", 1, 1)
cairo_pdf("03_figures/Fig6D_imm_corr.pdf", width = 4.6, height = max(3.0, 0.22 * length(sig) + 1.4))
pheatmap(t(cm), cluster_cols = TRUE, cluster_rows = FALSE, display_numbers = TRUE, number_format = "%.2f",
         number_color = "black", fontsize_number = 7, fontsize_row = 8,
         color = colorRampPalette(c("#0072B2", "white", "#D55E00"))(50),
         border_color = "grey85", angle_col = 45, main = NA)
dev.off()
file.remove("03_figures/tmp.pdf")
cat("saved: Fig6D_imm_corr.pdf\n")
cat("STEP28 DONE\n")
