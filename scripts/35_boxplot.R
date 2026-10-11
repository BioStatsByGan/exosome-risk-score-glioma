## 步骤 35：药敏预测结果的组间比较与箱线图（含血脑屏障穿透药重点标注）
## 输入：04_results/34.oncoPredict/calcPhenotype_Output/DrugPredictions.csv、risk.TCGA.txt
## 输出：04_results/35.drugSensitivity/{drug_diff.txt}、03_figures/Fig11A_drug_box.pdf、Fig11B_drug_volcano.pdf
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(ggplot2); library(ggpubr); library(reshape2)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT); source("scripts/plot_style.R")
dir.create("04_results/35.drugSensitivity", recursive = TRUE, showWarnings = FALSE)

dp <- fread("04_results/34.oncoPredict/calcPhenotype_Output/DrugPredictions.csv", check.names = FALSE)
ids <- dp[[1]]; dp <- as.data.frame(dp[, -1]); rownames(dp) <- ids
risk <- read.table("04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t",
                   check.names = FALSE, row.names = 1)
same <- intersect(rownames(dp), rownames(risk))
dp <- dp[same, , drop = FALSE]; risk <- risk[same, ]
cat("药敏矩阵:", nrow(dp), "样本 x", ncol(dp), "药物\n")

drugs <- colnames(dp)
res <- do.call(rbind, lapply(drugs, function(d) {
  hi <- as.numeric(dp[risk$Risk == "high", d]); lo <- as.numeric(dp[risk$Risk == "low", d])
  data.frame(Drug = d, IC50_high = mean(hi), IC50_low = mean(lo),
             log2FC = mean(hi) - mean(lo), p = wilcox.test(hi, lo)$p.value)
}))
res$FDR <- p.adjust(res$p, "BH")
res <- res[order(res$p), ]
fwrite(res, "04_results/35.drugSensitivity/drug_diff.txt", sep = "\t", quote = FALSE)
cat("显著药物（FDR<0.05）:", sum(res$FDR < 0.05), "/", nrow(res),
    "| 高危组更敏感（IC50 更低）:", sum(res$FDR < 0.05 & res$log2FC < 0),
    "| 高危组更耐药:", sum(res$FDR < 0.05 & res$log2FC > 0), "\n")
print(head(res, 10), row.names = FALSE, digits = 3)

## 一线/血脑屏障穿透药物是否在列
key <- grep("Temozolomide|Lomustine|Carmustine|Cisplatin|Vincristine|Cytarabine|Doxorubicin|Etoposide|Melphalan|Bortezomib|Irinotecan", drugs,
            value = TRUE, ignore.case = TRUE)
cat("\n重点药物（化疗/穿透 BBB）:\n")
print(res[res$Drug %in% key, ], row.names = FALSE, digits = 3)

## 图 11A：显著药物的 IC50 箱线图（前 20）
top <- head(res[res$FDR < 0.05, ], 20)
if (nrow(top)) {
  d <- reshape2::melt(as.matrix(dp[, top$Drug, drop = FALSE]))
  colnames(d) <- c("sample", "Drug", "IC50")
  d$Risk <- factor(rep(risk$Risk, ncol(dp[, top$Drug, drop = FALSE])), levels = c("low", "high"),
                   labels = c("Low risk", "High risk"))
  d$Drug <- factor(d$Drug, levels = top$Drug)
  ## 选择依据是 FDR（上一行 top <- ...），正文报告的是名义 P → 图内标注名义 P，保持口径一致
  lab <- setNames(sprintf("%s\nP = %s", top$Drug,
                          ifelse(top$p < 0.001, sprintf("%.1e", top$p), sprintf("%.3f", top$p))), top$Drug)
  p1 <- ggplot(d, aes(Drug, IC50, fill = Risk)) +
    geom_boxplot(width = 0.66, outlier.size = 0.25, linewidth = 0.28, colour = "grey25") +
    scale_fill_manual(values = c("Low risk" = "#0072B2", "High risk" = "#D55E00"), name = "") +
    scale_x_discrete(labels = lab) +
    labs(x = NULL, y = "Predicted IC50") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6.2), legend.position = "top")
  save_pdf(p1, "03_figures/Fig11A_drug_box.pdf", 11, 5.4)
}

## 图 11B：火山图（IC50 差异 vs 显著性）
res$sig <- ifelse(res$FDR < 0.05 & res$log2FC < 0, "Sensitive in high risk",
                  ifelse(res$FDR < 0.05 & res$log2FC > 0, "Resistant in high risk", "NS"))
lab <- rbind(head(res[res$sig != "NS", ][order(-abs(res[res$sig != "NS", ]$log2FC)), ], 10),
             res[res$Drug %in% key, ])
lab <- lab[!duplicated(lab$Drug), ]
p2 <- ggplot(res, aes(log2FC, -log10(p))) +
  geom_point(data = res[res$sig == "NS", ], colour = "grey80", size = 1) +
  geom_point(data = res[res$sig != "NS", ], aes(colour = sig), size = 1.6, alpha = 0.9) +
  scale_colour_manual(values = c("Sensitive in high risk" = "#0072B2", "Resistant in high risk" = "#D55E00", "NS" = "grey80"), name = "") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_hline(yintercept = -log10(max(res$p[res$FDR < 0.05])), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  ggrepel::geom_text_repel(data = lab, aes(label = Drug), size = 2.4, max.overlaps = 30,
                           segment.size = 0.2, min.segment.length = 0) +
  labs(x = expression(log[2]~"IC50 (high risk) - log"[2]~"IC50 (low risk)"), y = expression(-log[10]~"P")) +
  theme(legend.position = "top")
save_pdf(p2, "03_figures/Fig11B_drug_volcano.pdf", 8.6, 5.4)
cat("STEP35 DONE\n")
