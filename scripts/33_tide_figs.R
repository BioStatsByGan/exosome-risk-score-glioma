## 步骤 33 图：TIDE 各指标与风险分层的比较
## 输入：04_results/33.TIDE/TIDE.csv、04_results/20.SHAP/risk.TCGA.txt
## 输出：04_results/33.TIDE/tide_diff.txt、03_figures/Fig7A_tide.pdf、Fig7B_tide_response.pdf
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(ggplot2); library(ggpubr); library(reshape2); library(data.table)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT); source("scripts/plot_style.R")

tide <- read.csv("04_results/33.TIDE/TIDE.csv", row.names = 1, check.names = FALSE)
risk <- read.table("04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t",
                   check.names = FALSE, row.names = 1)
same <- intersect(rownames(tide), rownames(risk))
tide <- tide[same, ]; risk <- risk[same, ]
cat("配对样本:", length(same), "\n")

vars <- c("TIDE", "Dysfunction", "Exclusion", "MSI Score", "IFNG", "CD274", "CD8", "CTL",
          "MDSC", "CAF", "TAM M2")
res <- do.call(rbind, lapply(vars, function(v) {
  hi <- as.numeric(tide[risk$Risk == "high", v]); lo <- as.numeric(tide[risk$Risk == "low", v])
  data.frame(Indicator = v, mean_high = mean(hi, na.rm = TRUE), mean_low = mean(lo, na.rm = TRUE),
             p = wilcox.test(hi, lo)$p.value)
}))
res$FDR <- p.adjust(res$p, "BH")
res <- res[order(res$p), ]
fwrite(res, "04_results/33.TIDE/tide_diff.txt", sep = "\t", quote = FALSE)
cat("各指标组间比较（前 6）:\n"); print(head(res, 6), row.names = FALSE, digits = 3)

d <- reshape2::melt(data.frame(sample = rownames(tide), tide[, vars], Risk = risk$Risk, check.names = FALSE), id = c("sample", "Risk"))
colnames(d)[3] <- "Indicator"
d$Risk <- factor(d$Risk, levels = c("low", "high"), labels = c("Low risk", "High risk"))
d$Indicator <- factor(d$Indicator, levels = res$Indicator)
lab <- setNames(sprintf("%s\nP = %s", res$Indicator,
                        ifelse(res$p < 0.001, "<0.001", sprintf("%.3f", res$p))), res$Indicator)
p1 <- ggplot(d, aes(Indicator, value, fill = Risk)) +
  geom_boxplot(width = 0.66, outlier.size = 0.3, linewidth = 0.3, colour = "grey25") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  scale_fill_manual(values = c("Low risk" = "#0072B2", "High risk" = "#D55E00"), name = "") +
  scale_x_discrete(labels = lab) +
  labs(x = NULL, y = "TIDE-based indicator score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5), legend.position = "top")
save_pdf(p1, "03_figures/Fig7A_tide.pdf", 9.5, 5.2)

## 免疫治疗获益预测比例
d2 <- data.frame(Risk = risk$Risk, Benefit = ifelse(tide[, "No benefits"], "No benefit", "Benefit"))
tb <- table(d2)
pv <- chisq.test(tb)$p.value
cat("获益预测列联表:\n"); print(tb); cat("卡方 p =", signif(pv, 3), "\n")
pd <- as.data.frame(prop.table(tb, 1) * 100)
pd$Risk <- factor(pd$Risk, levels = c("low", "high"), labels = c("Low risk", "High risk"))
p2 <- ggplot(pd, aes(Risk, Freq, fill = Benefit)) +
  geom_col(width = 0.6, colour = "grey25", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", Freq)), position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = c(Benefit = "#0072B2", `No benefit` = "#D55E00"), name = "") +
  labs(x = NULL, y = "比例 (%)") + theme(legend.position = "top")
save_pdf(p2, "03_figures/Fig7B_tide_response.pdf", 4.6, 4.2)
cat("STEP33 FIG DONE\n")
