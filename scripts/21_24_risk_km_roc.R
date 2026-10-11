## 步骤 21-24：风险评分三联图 + KM + ROC（TCGA 训练队列）+ LASSO CV + 多因素森林图
## 输出：03_figures/Fig2A_lasso_cv.pdf, Fig2B_multicox_forest.pdf,
##       03_figures/Fig4A_riskscore.pdf, Fig4B_survstatus.pdf, Fig4C_km_TCGA.pdf, Fig4D_roc_TCGA.pdf
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({
  library(data.table); library(ggplot2); library(survival); library(survminer)
  library(glmnet); library(patchwork); library(grid)
  library(forestploter); library(timeROC)
})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
setwd(ROOT); source("scripts/plot_style.R")
dir.create("04_results/21_24", recursive = TRUE, showWarnings = FALSE)

d <- read.table("04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t",
                check.names = FALSE, row.names = 1)
d$Risk <- factor(d$Risk, levels = c("low", "high"))
cat("训练队列:", nrow(d), "例 | high", sum(d$Risk == "high"), " low", sum(d$Risk == "low"), "\n")

## ---- Fig2A LASSO 交叉验证曲线 ----
uni <- read.table("04_results/19.uniCox/TCGA.uniSigExp.txt", header = TRUE, sep = "\t",
                  check.names = FALSE, row.names = 1)
x <- as.matrix(uni[, -(1:2)]); y <- Surv(uni$futime, uni$fustat)
set.seed(20260919)
cvfit <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10)
cvdf <- data.frame(loglambda = log(cvfit$lambda), cvm = cvfit$cvm, cvsd = cvfit$cvsd,
                   nzero = cvfit$nzero)
pA <- ggplot(cvdf, aes(loglambda, cvm)) +
  geom_errorbar(aes(ymin = cvm - cvsd, ymax = cvm + cvsd), width = 0.06, colour = "grey55", linewidth = 0.3) +
  geom_point(size = 1.1, colour = "#D55E00") +
  geom_vline(xintercept = log(cvfit$lambda.min), linetype = "dashed", colour = "#0072B2", linewidth = 0.35) +
  geom_vline(xintercept = log(cvfit$lambda.1se), linetype = "dotted", colour = "grey40", linewidth = 0.35) +
  labs(x = expression(log(lambda)), y = "Partial likelihood deviance") +
  annotate("text", x = log(cvfit$lambda.min), y = max(cvdf$cvm + cvdf$cvsd),
           label = sprintf("lambda.min (%d genes)", cvfit$nzero[which.min(abs(cvdf$loglambda - log(cvfit$lambda.min)))]),
           hjust = -0.05, vjust = 1, size = 2.7)
save_pdf(pA, "03_figures/Fig2A_lasso_cv.pdf", 5.2, 4.2)

## ---- Fig2B 多因素 Cox 森林图（forestploter）----
mc <- read.table("04_results/20.SHAP/multiCox.txt", header = TRUE, sep = "\t", check.names = FALSE)
mc <- mc[order(mc$HR), ]
mc$`Hazard ratio (95% CI)` <- sprintf("%.2f (%.2f-%.2f)", mc$HR, mc$HR.95L, mc$HR.95H)
mc$`P value` <- ifelse(mc$pvalue < 0.001, "<0.001", sprintf("%.3f", mc$pvalue))
mc$` ` <- paste(rep(" ", 22), collapse = " ")
tab <- mc[, c("id", "Hazard ratio (95% CI)", "P value", " ")]
names(tab)[1] <- "Gene"
tm <- forest_theme(base_size = 10, refline_col = "grey45", ci_col = "#0072B2",
                   ci_lty = 1, ci_lwd = 1.5, ci_Theight = 0.22,
                   row_fill = c("white", "grey96"), core = list(bg_params = list(fill = c("white", "grey97"), alpha = 1)))
pB <- forest(tab, est = mc$HR, lower = mc$HR.95L, upper = mc$HR.95H,
             sizes = 0.9, ci_column = 4, ref_line = 1, xlim = c(0, 2.2),
             ticks_at = c(0.5, 1, 1.5, 2), theme = tm)
save_pdf(pB, "03_figures/Fig2B_multicox_forest.pdf", 6.4, 4.6)

## ---- Fig4A/B/C/D 风险评分 + 生存状态 + KM + ROC ----
ds <- d[order(d$riskScore), ]
ds$rank <- seq_len(nrow(ds))
p1 <- ggplot(ds, aes(rank, riskScore)) +
  geom_point(aes(colour = Risk), size = 0.6) +
  geom_vline(xintercept = sum(ds$Risk == "low"), linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  scale_colour_manual(values = c("low" = "#0072B2", "high" = "#D55E00"), name = "Risk") +
  labs(x = "Patients (ranked by risk score)", y = "Risk score") +
  theme(legend.position = c(0.25, 0.85))
ds$status <- ifelse(ds$fustat == 1, "Dead", "Alive")
p2 <- ggplot(ds, aes(rank, futime)) +
  geom_point(aes(colour = status), size = 0.7) +
  geom_vline(xintercept = sum(ds$Risk == "low"), linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  scale_colour_manual(values = c("Alive" = "#0072B2", "Dead" = "#D55E00"), name = "Status") +
  labs(x = "Patients (ranked by risk score)", y = "Survival time (years)") +
  theme(legend.position = c(0.25, 0.85))
pkm <- km_pub(d)$plot + labs(title = NULL)
p4 <- roc_step(d, times = c(1, 3, 5), marker_col = "riskScore")
save_pdf(p1, "03_figures/Fig4A_riskscore.pdf", 5.6, 4.0)
save_pdf(p2, "03_figures/Fig4B_survstatus.pdf", 5.6, 4.0)
save_pdf(pkm, "03_figures/Fig4C_km_TCGA.pdf", 5.6, 4.6)
save_pdf(p4, "03_figures/Fig4D_roc_TCGA.pdf", 5.6, 4.6)

## 训练集统计
sd <- survdiff(Surv(futime, fustat) ~ Risk, data = d)
cox <- coxph(Surv(futime, fustat) ~ riskScore, data = d)
r <- timeROC(T = d$futime, delta = d$fustat, marker = d$riskScore, cause = 1, times = c(1, 3, 5), iid = TRUE)
cis <- confint(r, parm = NULL, n.sim = 500); cis$CI <- cis$CI / 100
stat <- data.frame(
  cohort = "TCGA (train)", n = nrow(d),
  HR = sprintf("%.2f (%.2f-%.2f)", exp(coef(cox)), exp(confint(cox))[1], exp(confint(cox))[2]),
  p = summary(cox)$coefficients[, "Pr(>|z|)"],
  Cindex = summary(cox)$concordance[1],
  logrank_p = 1 - pchisq(sd$chisq, 1),
  AUC1 = r$AUC[1], AUC3 = r$AUC[2], AUC5 = r$AUC[3])
fwrite(stat, "04_results/21_24/train_stats.txt", sep = "\t", quote = FALSE)
print(stat, row.names = FALSE, digits = 4)
cat("STEP21-24 DONE\n")
