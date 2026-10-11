## 审计补救 B 组：决策曲线（DCA）、时间依赖 Brier、训练集切点套外部队列
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(survival); library(rms); library(pec); library(reshape2)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT); source("scripts/plot_style.R")
dir.create("04_results/audit", recursive = TRUE, showWarnings = FALSE)
## 控制台输出同时落盘（含 DCA 纳入例数、Brier/IPA、训练切点分组计数）
sink("04_results/audit/audit_b_group_log.txt", split = TRUE)

risk <- read.table("04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t",
                   check.names = FALSE, row.names = 1)
cli <- read.table("02_inputs/clinical_glioma.txt", header = TRUE, sep = "\t", check.names = FALSE)
d <- merge(data.frame(patient = rownames(risk), risk, check.names = FALSE),
           cli[, c("patient", "age", "sex", "grade", "idh")], by = "patient")
d$Grade <- factor(d$grade); d$IDH <- factor(d$idh); d$Sex <- factor(d$sex)
units(d$futime) <- "Year"

## ---------- B1 决策曲线（3 年死亡的二分类版本）----------
cat("========== B1 决策曲线（3 年死亡，限制在随访信息完整的样本）==========\n")
d3 <- d[(d$futime >= 3) | (d$fustat == 1 & d$futime < 3), ]
d3$death3 <- ifelse(d3$fustat == 1 & d3$futime <= 3, 1, 0)
d3 <- d3[complete.cases(d3[, c("death3", "riskScore", "age", "Grade", "IDH")]), ]
cat("纳入样本:", nrow(d3), "| 3 年内死亡:", sum(d3$death3),
    "（定义：随访≥3 年或 3 年内死亡，且 death3/riskScore/age/Grade/IDH 完整）\n")
set.seed(1)
m <- glm(death3 ~ riskScore + age + Grade + IDH, data = d3, family = binomial)
d3$pred <- predict(m, newdata = d3, type = "response")
## 仅临床变量模型（年龄 + 分级 + IDH）——判断净获益是否来自新增评分
m_cli <- glm(death3 ~ age + Grade + IDH, data = d3, family = binomial)
d3$pred_cli <- predict(m_cli, newdata = d3, type = "response")
## 二分类 DCA 标准公式：NB = TP/n - FP/n * pt/(1-pt)
n <- nrow(d3); P <- sum(d3$death3)
nb_fun <- function(pt) {
  tr <- d3$pred >= pt
  TP <- sum(tr & d3$death3 == 1); FP <- sum(tr & d3$death3 == 0)
  nb_model <- TP / n - FP / n * (pt / (1 - pt))
  trc <- d3$pred_cli >= pt
  TPc <- sum(trc & d3$death3 == 1); FPc <- sum(trc & d3$death3 == 0)
  nb_cli <- TPc / n - FPc / n * (pt / (1 - pt))
  nb_all <- P / n - (n - P) / n * (pt / (1 - pt))
  c(model = nb_model, clinical = nb_cli, all = nb_all, none = 0)
}
pts <- seq(0.05, 0.90, by = 0.01)
nb <- as.data.frame(t(sapply(pts, nb_fun)))
nb$thresholds <- pts
fwrite(nb, "04_results/audit/dca_3year.txt", sep = "\t", quote = FALSE)
cat("净获益（阈值 0.3/0.5/0.7）:\n")
print(nb[nb$thresholds %in% c(0.3, 0.5, 0.7), c("model", "clinical", "all", "none", "thresholds")], row.names = FALSE, digits = 3)
nb2 <- reshape2::melt(nb, id.vars = "thresholds", variable.name = "model", value.name = "nb")
p <- ggplot(nb2, aes(thresholds, nb, colour = model)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c(model = "#0072B2", clinical = "#009E73", all = "#D55E00", none = "grey55"),
                      labels = c("Full model", "Clinical only", "Treat all", "Treat none"), name = "") +
  labs(x = "Threshold probability", y = "Net benefit") +
  theme(legend.position = c(0.75, 0.85))
save_pdf(p, "03_figures/Fig12_dca.pdf", 4.8, 4.4)

## ---------- B2 时间依赖 Brier 分数与 IPA ----------
cat("\n========== B2 时间依赖 Brier 分数 ==========\n")
d <- d[complete.cases(d[, c("futime", "fustat", "riskScore", "age", "Grade", "IDH")]), ]
cat("Brier 分析纳入完整病例:", nrow(d), "\n")
dd <- datadist(d); options(datadist = "dd")
f <- cph(Surv(futime, fustat) ~ riskScore + age + Grade + IDH, data = d, x = TRUE, y = TRUE, surv = TRUE)
## Graf 型 IPW Brier 分数（手算，公式透明）
kmc <- survfit(Surv(futime, 1 - fustat) ~ 1, data = d)   # 删失分布
Gt <- function(t) {
  s <- summary(kmc, times = t, extend = TRUE)$surv
  ifelse(is.na(s) | s <= 0, min(s[!is.na(s) & s > 0], 1), s)
}
## 预测生存率用基线累积风险 × exp(lp) 手算（survest 对缺失协变量会返回 NA）
bh <- basehaz(f, centered = FALSE)
H0 <- approx(bh$time, bh$hazard, xout = c(1, 3, 5), rule = 2)$y
lp <- predict(f, type = "lp")
surv_pred <- exp(-outer(exp(lp), H0))
resB <- do.call(rbind, lapply(seq_along(c(1, 3, 5)), function(k) {
  t <- c(1, 3, 5)[k]
  Shat <- surv_pred[, k]
  ev <- d$fustat == 1 & d$futime <= t
  at_risk <- d$futime > t
  w_ev <- ifelse(ev, 1 / Gt(d$futime), 0)
  bs_full <- mean(w_ev * (0 - Shat)^2 + ifelse(at_risk, 1 / Gt(t), 0) * (1 - Shat)^2)
  km <- summary(survfit(Surv(futime, fustat) ~ 1, data = d), times = t, extend = TRUE)$surv
  bs_null <- mean(w_ev * (0 - km)^2 + ifelse(at_risk, 1 / Gt(t), 0) * (1 - km)^2)
  data.frame(Time = t, Brier_full = round(bs_full, 3), Brier_null = round(bs_null, 3),
             IPA = round(1 - bs_full / bs_null, 3))
}))
print(resB, row.names = FALSE)
fwrite(resB, "04_results/audit/brier_scores.txt", sep = "\t", quote = FALSE)

## ---------- B3 训练集切点套外部队列 ----------
cat("\n========== B3 训练集 cutpoint 直接套 CGGA 两队列 ==========\n")
cutoff <- as.numeric(readLines("04_results/20.SHAP/cutoff.txt"))
cat("训练集 cutoff =", round(cutoff, 3), "（CGGA 验证使用的是各队列中位数）\n")
for (nm in c("325", "693")) {
  cg <- read.table(sprintf("04_results/18.CGGA/risk.CGGA%s.txt", nm), header = TRUE, sep = "\t",
                   check.names = FALSE)
  cg$RiskT <- ifelse(cg$riskScore > cutoff, "high", "low")
  tb <- table(cg$RiskT)
  cox <- coxph(Surv(futime, fustat) ~ riskScore, data = cg)
  sdT <- survdiff(Surv(futime, fustat) ~ RiskT, data = cg)
  cat(sprintf("CGGA-%s：训练切点分组 high=%d low=%d | HR=%.2f (%.2f-%.2f) | log-rank P=%.3g | C-index=%.3f\n",
              nm, tb["high"], tb["low"], exp(coef(cox)), exp(confint(cox))[1], exp(confint(cox))[2],
              1 - pchisq(sdT$chisq, 1), summary(cox)$concordance[1]))
  med <- median(cg$riskScore)
  cat(sprintf("        该队列中位数为 %.3f，训练切点与之差 %.3f（差距越大说明跨平台分值不可直接搬）\n",
              med, cutoff - med))
}
sink()
cat("written: 04_results/audit/audit_b_group_log.txt\n")
cat("B-GROUP DONE\n")
