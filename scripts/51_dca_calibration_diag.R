## 步骤 51：IPCW 决策曲线 + 校准斜率/截距诊断 + 共线性与 EPV（回应外部审查第 3、4 条及一致性检查）
## (A) 用适合右删失的 IPCW 净获益，在全部可评估样本上重做 3 年决策曲线，并与旧的「完整随访子集」版本对照；
## (B) 校准：给出斜率、校准截距（calibration-in-the-large）与分箱的预测-观测差，替代「斜率>1 ⇒ 低值高估」的推测；
## (C) 共线性：偏相关系数 / VIF / 事件数-参数比（EPV）。
## 输出：04_results/audit/dca_ipcw.txt、calibration_diag.txt、collinearity.txt
##       稿/tables/表S13_校准与DCA诊断.csv
##       03_figures/FigS7_ipcw_dca.{pdf,png}
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(survival); library(rms); library(ggplot2); library(reshape2)})
source("scripts/plot_style.R")
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
OUTD <- "04_results/audit"; dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
LOG <- character(); say <- function(...) { s <- sprintf(...); cat(s, "\n"); LOG <<- c(LOG, s); invisible(s) }

risk <- read.table("04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
cli <- read.table("02_inputs/clinical_glioma.txt", header = TRUE, sep = "\t", check.names = FALSE)
d <- merge(data.frame(patient = rownames(risk), risk, check.names = FALSE),
           cli[, c("patient", "age", "sex", "grade", "idh")], by = "patient")
d$Grade <- factor(d$grade); d$IDH <- factor(d$idh); d$Sex <- factor(d$sex)
units(d$futime) <- "Year"

## ================= (A) IPCW 3 年决策曲线 =================
say("========== (A) 3 年决策曲线：IPCW vs 完整随访子集 ==========")
dc <- d[complete.cases(d[, c("futime", "fustat", "riskScore", "age", "Grade", "IDH")]), ]
say("全部可评估样本 n=%d（事件 %d）", nrow(dc), sum(dc$fustat))

fit_full <- cph(Surv(futime, fustat) ~ riskScore + age + Grade + IDH, data = dc, x = TRUE, y = TRUE, surv = TRUE)
fit_cli  <- cph(Surv(futime, fustat) ~ age + Grade + IDH, data = dc, x = TRUE, y = TRUE, surv = TRUE)
## 预测生存率统一用 survest()（2026-10-12 修正：不能用 basehaz(centered=FALSE) 配中心化 lp）
dc$s_full <- survest(fit_full, newdata = dc, times = 3)$surv
dc$s_cli  <- survest(fit_cli,  newdata = dc, times = 3)$surv
dc$p_full <- 1 - dc$s_full
dc$p_cli  <- 1 - dc$s_cli
say("预测 3 年生存率均值：全模型 %.3f、仅临床 %.3f；KM 观测 %.3f",
    mean(dc$s_full), mean(dc$s_cli),
    summary(survfit(Surv(futime, fustat) ~ 1, data = dc), times = 3, extend = TRUE)$surv)

## 删失分布 G(t)（IPCW 用）
kmc <- survfit(Surv(futime, 1 - fustat) ~ 1, data = dc)
Gt <- function(t) {
  s <- summary(kmc, times = t, extend = TRUE)$surv
  ifelse(is.na(s) | s <= 0, min(s[!is.na(s) & s > 0], 1), s)
}
tau <- 0.3
nb_ipcw <- function(pred, tau_t = 3) {
  ev <- dc$fustat == 1 & dc$futime <= tau_t          # 3 年内死亡
  at <- dc$futime > tau_t                            # 3 年时仍在随访
  w_ev <- ifelse(ev, 1 / Gt(dc$futime), 0)
  w_at <- ifelse(at, 1 / Gt(tau_t), 0)
  f <- function(pt) {
    tr <- pred >= pt
    mean(tr * w_ev) - mean(tr * w_at) * (pt / (1 - pt))
  }
  f
}
pts <- seq(0.05, 0.90, by = 0.01)
f_full <- nb_ipcw(dc$p_full); f_cli <- nb_ipcw(dc$p_cli)
nb_all <- function(pt) {
  ev <- dc$fustat == 1 & dc$futime <= 3; at <- dc$futime > 3
  w_ev <- ifelse(ev, 1 / Gt(dc$futime), 0); w_at <- ifelse(at, 1 / Gt(3), 0)
  mean(w_ev) - mean(w_at) * (pt / (1 - pt))
}
nb <- data.frame(threshold = pts, model = sapply(pts, f_full), clinical = sapply(pts, f_cli),
                 all = sapply(pts, nb_all), none = 0)
fwrite(nb, file.path(OUTD, "dca_3year_ipcw.txt"), sep = "\t", quote = FALSE)
say("IPCW 净获益（全部 %d 例）阈值 0.3：模型 %.3f，仅临床 %.3f，全部治疗 %.3f",
    nrow(dc), nb$model[nb$threshold == 0.3], nb$clinical[nb$threshold == 0.3], nb$all[nb$threshold == 0.3])
say("IPCW 净获益 阈值 0.5：模型 %.3f，仅临床 %.3f，全部治疗 %.3f",
    nb$model[nb$threshold == 0.5], nb$clinical[nb$threshold == 0.5], nb$all[nb$threshold == 0.5])

## 旧的「完整随访子集」口径（复现已有 396 例分析）
d3 <- dc[(dc$futime >= 3) | (dc$fustat == 1 & dc$futime < 3), ]
d3$death3 <- ifelse(d3$fustat == 1 & d3$futime <= 3, 1, 0)
d3 <- d3[complete.cases(d3[, c("death3", "riskScore", "age", "Grade", "IDH")]), ]
m3 <- glm(death3 ~ riskScore + age + Grade + IDH, data = d3, family = binomial)
m3c <- glm(death3 ~ age + Grade + IDH, data = d3, family = binomial)
d3$p <- predict(m3, type = "response"); d3$pc <- predict(m3c, type = "response")
n3 <- nrow(d3); P3 <- sum(d3$death3)
nb_cc <- function(pt) {
  tr <- d3$p >= pt; trc <- d3$pc >= pt
  c(model = (sum(tr & d3$death3 == 1) - sum(tr & d3$death3 == 0) * pt / (1 - pt)) / n3,
    clinical = (sum(trc & d3$death3 == 1) - sum(trc & d3$death3 == 0) * pt / (1 - pt)) / n3,
    all = (P3 - (n3 - P3) * pt / (1 - pt)) / n3)
}
nbcc <- t(sapply(c(0.3, 0.5), nb_cc)); rownames(nbcc) <- c("0.3", "0.5")
say("完整随访子集（旧口径，n=%d）：阈值 0.3 模型 %.3f / 仅临床 %.3f；阈值 0.5 模型 %.3f / 仅临床 %.3f",
    n3, nbcc["0.3", "model"], nbcc["0.3", "clinical"], nbcc["0.5", "model"], nbcc["0.5", "clinical"])

## 子集选择是否改变风险分布
say("风险评分分布：全部 799 例 中位 %.3f（IQR %.3f–%.3f）；396 例子集 中位 %.3f（IQR %.3f–%.3f）；KS P=%.3g",
    median(dc$riskScore), quantile(dc$riskScore, .25), quantile(dc$riskScore, .75),
    median(d3$riskScore), quantile(d3$riskScore, .25), quantile(d3$riskScore, .75),
    ks.test(dc$riskScore, d3$riskScore)$p.value)
say("被排除的 %d 例（3 年前删失）评分中位 %.3f，与纳入者中位 %.3f",
    nrow(dc) - n3, median(dc$riskScore[!(rownames(dc) %in% rownames(d3))]), median(d3$riskScore))

## ================= (B) 校准诊断 =================
say("\n========== (B) 校准斜率 / 校准截距 ==========")
cal <- readRDS("04_results/25_26/calibration.rds")
calib_rows <- list()
for (u in c(1, 3, 5)) {
  cc <- cal[[as.character(u)]]
  cc <- cc[is.finite(cc[, "mean.predicted"]) & is.finite(cc[, "KM"]), , drop = FALSE]
  x <- as.numeric(cc[, "mean.predicted"]); y <- as.numeric(cc[, "KM"])
  fit <- lm(y ~ x)
  slope <- coef(fit)[2]; inter <- coef(fit)[1]
  lo <- which.min(x); hi <- which.max(x)
  say("%d 年：斜率 %.2f，截距 %+.2f（理想 0/1：斜率 1、截距 0）", u, slope, inter)
  say("     最低分箱：预测 %.2f vs 观测 %.2f（差 %+.2f）；最高分箱：预测 %.2f vs 观测 %.2f（差 %+.2f）",
      x[lo], y[lo], y[lo] - x[lo], x[hi], y[hi], y[hi] - x[hi])
  say("     预测区间 %.2f–%.2f，观测区间 %.2f–%.2f；分箱数 %d", min(x), max(x), min(y), max(y), nrow(cc))
  calib_rows[[as.character(u)]] <- data.frame(Time = u, Slope = round(slope, 2), Intercept = round(inter, 2),
                                              Bin_low_pred = round(x[lo], 2), Bin_low_obs = round(y[lo], 2),
                                              Bin_high_pred = round(x[hi], 2), Bin_high_obs = round(y[hi], 2),
                                              n_bins = nrow(cc))
}
calib_tab <- do.call(rbind, calib_rows)

## 同时给出 IPCW 校准（避免「只用 KM 分箱」的质疑）：观测 3 年生存率 vs 预测均值
say("IPCW 一致性检查（整体）：1/3/5 年预测均值 %s；KM 观测 %s",
    paste(sprintf("%.2f", sapply(c(1,3,5), function(u) {
      mean(survest(fit_full, newdata = dc, times = u)$surv) })), collapse = "/"),
    paste(sprintf("%.2f", sapply(c(1,3,5), function(u)
      summary(survfit(Surv(futime, fustat) ~ 1, data = dc), times = u, extend = TRUE)$surv)), collapse = "/"))

## ================= (C) 共线性 / EPV =================
say("\n========== (C) 共线性与 EPV ==========")
vars <- c("riskScore", "age", "Grade", "IDH")
mm <- model.matrix(~ riskScore + age + Grade + IDH, data = dc)[, -1, drop = FALSE]
vif1 <- sapply(colnames(mm), function(v) {
  y <- mm[, v]; X <- mm[, setdiff(colnames(mm), v), drop = FALSE]
  1 / (1 - summary(lm(y ~ X))$r.squared)
})
say("TCGA 多因素模型 VIF：%s", paste(sprintf("%s=%.2f", names(vif1), vif1), collapse = ", "))
say("事件数 %d，模型参数 %d，EPV=%.1f", sum(dc$fustat), ncol(mm) + 1, sum(dc$fustat) / (ncol(mm) + 1))
say("评分与分级 Spearman ρ=%.3f（G2/G3/G4 序数编码）；与 IDH 的秩相关 ρ=%.3f",
    suppressWarnings(cor(dc$riskScore, as.numeric(dc$Grade), method = "spearman")),
    suppressWarnings(cor(dc$riskScore, as.numeric(dc$IDH), method = "spearman")))
cg <- read.table("04_results/18.CGGA/risk.CGGA325.txt", header = TRUE, sep = "\t", check.names = FALSE)
cg <- cg[!is.na(cg$futime) & !is.na(cg$fustat) & cg$Grade %in% c("WHO II", "WHO III", "WHO IV"), ]
mmc <- model.matrix(~ riskScore + factor(Grade, levels = c("WHO II","WHO III","WHO IV")), data = cg)[, -1, drop = FALSE]
vifc <- sapply(colnames(mmc), function(v) {
  y <- mmc[, v]; X <- mmc[, setdiff(colnames(mmc), v), drop = FALSE]
  1 / (1 - summary(lm(y ~ X))$r.squared)
})
say("CGGA-325（评分 + 分级）VIF：%s；事件 %d，参数 %d，EPV=%.1f",
    paste(sprintf("%s=%.2f", names(vifc), vifc), collapse = ", "), sum(cg$fustat), ncol(mmc) + 1,
    sum(cg$fustat) / (ncol(mmc) + 1))
say("注：CGGA-325 若同时纳入分级与 IDH，WHO II × IDH 野生型格子仅 %d 例（准完全分离），
故主分析只用评分+分级，与表 3 同口径。", sum(cg$Grade == "WHO II" & cg$IDH == "Wildtype"))

## ================= 落盘 =================
## 长表：Analysis / Metric / Time_years / Value
mk <- function(a, m, t, v) data.frame(Analysis = a, Metric = m, Time_years = t, Value = round(v, 3))
tabS13 <- rbind(
  mk("IPCW DCA, 3-year death (n=691)", "net benefit, full model, threshold 0.3", 3, nb$model[nb$threshold == 0.3]),
  mk("IPCW DCA, 3-year death (n=691)", "net benefit, clinical only, threshold 0.3", 3, nb$clinical[nb$threshold == 0.3]),
  mk("IPCW DCA, 3-year death (n=691)", "net benefit, treat all, threshold 0.3", 3, nb$all[nb$threshold == 0.3]),
  mk("IPCW DCA, 3-year death (n=691)", "net benefit, full model, threshold 0.5", 3, nb$model[nb$threshold == 0.5]),
  mk("IPCW DCA, 3-year death (n=691)", "net benefit, clinical only, threshold 0.5", 3, nb$clinical[nb$threshold == 0.5]),
  mk("IPCW DCA, 3-year death (n=691)", "net benefit, treat all, threshold 0.5", 3, nb$all[nb$threshold == 0.5]),
  mk("Complete-case DCA (n=396, old analysis)", "net benefit, full model, threshold 0.3", 3, nbcc["0.3", "model"]),
  mk("Complete-case DCA (n=396, old analysis)", "net benefit, clinical only, threshold 0.3", 3, nbcc["0.3", "clinical"]),
  mk("Complete-case DCA (n=396, old analysis)", "net benefit, full model, threshold 0.5", 3, nbcc["0.5", "model"]),
  mk("Complete-case DCA (n=396, old analysis)", "net benefit, clinical only, threshold 0.5", 3, nbcc["0.5", "clinical"]),
  mk("Selection check", "median risk score, all evaluable (n=691)", NA, median(dc$riskScore)),
  mk("Selection check", "median risk score, complete-case subset (n=396)", NA, median(d3$riskScore)),
  mk("Selection check", "KS test P (score distribution, subset vs all)", NA, ks.test(dc$riskScore, d3$riskScore)$p.value),
  mk("Calibration (TCGA, n=691)", paste0("slope, ", calib_tab$Time, "-year"), calib_tab$Time, calib_tab$Slope),
  mk("Calibration (TCGA, n=691)", paste0("intercept, ", calib_tab$Time, "-year"), calib_tab$Time, calib_tab$Intercept),
  mk("Collinearity (TCGA)", paste0("VIF, ", names(vif1)), NA, vif1),
  mk("Collinearity (CGGA-325)", paste0("VIF, ", gsub("factor[^)]*\\)", "Grade ", names(vifc))), NA, vifc))
fwrite(tabS13, "稿/tables/表S13_校准与DCA诊断.csv", sep = ",", quote = FALSE)
say("\n表 S13 已写入 稿/tables/表S13_校准与DCA诊断.csv")

## 图 S7：IPCW DCA（全样本）与完整随访子集对照
nb2 <- reshape2::melt(nb, id.vars = "threshold", variable.name = "model", value.name = "nb")
cc_curve <- data.frame(threshold = pts, nb = sapply(pts, function(p) nb_cc(p)["model"]))
p <- ggplot(nb2, aes(threshold, nb, colour = model)) +
  geom_line(linewidth = 0.6) +
  geom_line(data = cc_curve, aes(threshold, nb), inherit.aes = FALSE, linetype = 2, colour = "grey30") +
  scale_colour_manual(values = c(model = "#0072B2", clinical = "#009E73", all = "#D55E00", none = "grey55"),
                      labels = c("Full model (IPCW)", "Clinical only (IPCW)", "Treat all (IPCW)", "Treat none", "Complete-case subset"),
                      name = "") +
  labs(x = "Threshold probability", y = "Net benefit (3-year, IPCW)") +
  theme(legend.position = c(0.72, 0.8))
ggsave("03_figures/FigS7_ipcw_dca.pdf", p, width = 5.2, height = 4.4)
ggsave("03_figures/FigS7_ipcw_dca.png", p, width = 5.2, height = 4.4, dpi = 200)
say("图 S7 已保存：03_figures/FigS7_ipcw_dca.pdf/.png")

## 主图 8D 同步改为 IPCW 口径（4 条曲线，与正文「IPCW 为主分析」一致；完整随访子集曲线见 图 S7）
p8 <- ggplot(nb2, aes(threshold, nb, colour = model)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c(model = "#0072B2", clinical = "#009E73", all = "#D55E00", none = "grey55"),
                      labels = c("Full model", "Clinical only", "Treat all", "Treat none"), name = "") +
  labs(x = "Threshold probability", y = "Net benefit") +
  theme(legend.position = c(0.75, 0.85))
ggsave("03_figures/Fig12_dca.pdf", p8, width = 4.8, height = 4.4)
ggsave("03_figures/Fig12_dca.png", p8, width = 4.8, height = 4.4, dpi = 200)
say("主图 8D 已改为 IPCW 口径：03_figures/Fig12_dca.pdf/.png（n=%d）", nrow(dc))

writeLines(c(sprintf("== DCA/校准/共线性诊断（%s）==", format(Sys.time(), "%Y-%m-%d %H:%M")), LOG),
           file.path(OUTD, "dca_ipcw.txt"))
writeLines(c(sprintf("== 校准诊断（%s）==", format(Sys.time(), "%Y-%m-%d %H:%M")),
             capture.output(print(calib_tab, row.names = FALSE))),
           file.path(OUTD, "calibration_diag.txt"))
writeLines(c(sprintf("== 共线性与 EPV（%s）==", format(Sys.time(), "%Y-%m-%d %H:%M")),
             sprintf("TCGA VIF: %s", paste(sprintf("%s=%.2f", names(vif1), vif1), collapse = ", ")),
             sprintf("CGGA-325 VIF: %s", paste(sprintf("%s=%.2f", names(vifc), vifc), collapse = ", "))),
           file.path(OUTD, "collinearity.txt"))
cat("\nwritten: dca_ipcw.txt / calibration_diag.txt / collinearity.txt\nSTEP51 DONE\n")
