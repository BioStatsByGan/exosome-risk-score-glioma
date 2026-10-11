## 步骤 31 + 32：TMB 与风险分层（31）及 TMB×风险联合生存（32）
## 输入：04_results/31.riskTMB/TMB.txt、04_results/20.SHAP/risk.TCGA.txt
## 输出：04_results/31.riskTMB/tmb_risk.txt、tmb_surv.txt、03_figures/Fig8B_tmb.pdf、Fig8C_tmb_surv.pdf
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(ggplot2); library(ggpubr); library(survival); library(survminer)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT); source("scripts/plot_style.R")
## 控制台输出同时落盘（风险评分与 log2(TMB) 的 Spearman 相关、四组 log-rank 等关键数字）
sink("04_results/31.riskTMB/run.log", split = TRUE)

tmb <- fread("04_results/31.riskTMB/TMB.txt")
risk <- read.table("04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t",
                   check.names = FALSE, row.names = 1)
m <- merge(tmb, data.frame(patient = rownames(risk), riskScore = risk$riskScore, Risk = risk$Risk,
                           futime = risk$futime, fustat = risk$fustat), by = "patient")
cat("TMB 与风险评分匹配样本:", nrow(m), "\n")
m$logTMB <- log2(m$TMB + 1)
m$Risk <- factor(m$Risk, levels = c("low", "high"), labels = c("Low risk", "High risk"))

## 组间比较与相关性
pv <- wilcox.test(logTMB ~ Risk, data = m)$p.value
ct <- cor.test(m$riskScore, m$logTMB, method = "spearman")
## 注意：生存/检验在 log2(TMB+1) 尺度上做，但报告的组间中位 TMB 用原始 TMB（突变/Mb）尺度
cat(sprintf("TMB：高风险组中位 %.2f vs 低风险组 %.2f 突变/Mb，Wilcoxon P(基于 log2(TMB+1)) = %.3g\n",
            median(m$TMB[m$Risk == "High risk"]), median(m$TMB[m$Risk == "Low risk"]), pv))
cat(sprintf("      log2(TMB+1) 中位：高风险 %.2f vs 低风险 %.2f\n",
            median(m$logTMB[m$Risk == "High risk"]), median(m$logTMB[m$Risk == "Low risk"])))
cat(sprintf("风险评分与 log2(TMB) Spearman rho = %.3f, P = %.3g\n", ct$estimate, ct$p.value))
fwrite(m, "04_results/31.riskTMB/tmb_risk.txt", sep = "\t", quote = FALSE)

p1 <- ggplot(m, aes(Risk, logTMB, fill = Risk)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.35, colour = "grey25") +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.4, colour = "grey35") +
  scale_fill_manual(values = c("Low risk" = "#0072B2", "High risk" = "#D55E00"), guide = "none") +
  annotate("text", x = 1.5, y = max(m$logTMB) * 1.02, size = 3.2,
           label = sprintf("Wilcoxon P = %s", ifelse(pv < 0.001, "<0.001", sprintf("%.3f", pv)))) +
  labs(x = NULL, y = expression(log[2]~"(TMB + 1), mutations/Mb"))
save_pdf(p1, "03_figures/Fig8B_tmb.pdf", 4.6, 4.4)

## 步骤 32：TMB 高低 × 风险高低 四组生存
m$TMBgroup <- ifelse(m$logTMB > median(m$logTMB), "High TMB", "Low TMB")
m$Group <- paste(m$Risk, m$TMBgroup, sep = " + ")
s <- summary(coxph(Surv(futime, fustat) ~ riskScore + logTMB, data = m))
cat("\n多因素 Cox（风险评分 + log2TMB）:\n")
print(data.frame(Variable = rownames(s$coefficients), HR = sprintf("%.2f (%.2f-%.2f)", s$conf.int[, 1], s$conf.int[, 3], s$conf.int[, 4]),
                 P = signif(s$coefficients[, 5], 3), row.names = NULL), row.names = FALSE)
fwrite(data.frame(Variable = rownames(s$coefficients), HR = s$conf.int[, 1], L = s$conf.int[, 3],
                  H = s$conf.int[, 4], p = s$coefficients[, 5]),
       "04_results/31.riskTMB/tmb_surv.txt", sep = "\t", quote = FALSE)

m$Group <- factor(m$Group, levels = c("Low risk + Low TMB", "Low risk + High TMB",
                                      "High risk + Low TMB", "High risk + High TMB"))
kv <- survfit(Surv(futime, fustat) ~ Group, data = m)
sd <- survdiff(Surv(futime, fustat) ~ Group, data = m)
pv2 <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
cat(sprintf("\n四组 log-rank P = %.3g\n", pv2))
g <- ggsurvplot(kv, data = m, pval = fmt_p_lab(pv2, "Log-rank P"), pval.size = 3.4,
                risk.table = TRUE, risk.table.height = 0.28, censor.size = 1.6, size = 0.7,
                palette = c("#0072B2", "#56B4E9", "#E69F00", "#D55E00"),
                legend.title = "", legend.labs = levels(m$Group),
                xlab = "Time (years)", ylab = "Overall survival probability",
                ggtheme = theme_pub(), break.time.by = 2, xlim = c(0, 10))
g$plot <- g$plot + labs(title = NULL) + theme(legend.text = element_text(size = 7.5))
save_pdf(g$plot, "03_figures/Fig8C_tmb_surv.pdf", 6.4, 5.2)
sink()
cat("written: 04_results/31.riskTMB/run.log\n")
cat("STEP31-32 DONE\n")
