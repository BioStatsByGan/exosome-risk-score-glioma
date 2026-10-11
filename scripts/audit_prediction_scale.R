## 审计脚本：预测生存率的两种算法对照，证明 Brier/IPA 之前为何被高估（2026-10-12）
## 结论：rms::cph 的 predict(type="lp") 是「按 datadist 均值中心化」的线性预测子，
##       必须配 basehaz(centered = TRUE) 或直接用 survest()；
##       旧代码用 basehaz(centered = FALSE)（对应协变量全取 0 的基线风险）配中心化 lp，
##       使预测生存率整体抬高（3 年 0.76 vs KM 0.52），导致 Brier 被高估、IPA 被低估。
## 输出：04_results/audit/prediction_scale_check.txt
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(survival); library(rms)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
OUTD <- "04_results/audit"; dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
LOG <- character(); say <- function(...) { s <- sprintf(...); cat(s, "\n"); LOG <<- c(LOG, s); invisible(s) }

risk <- read.table("04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
cli <- read.table("02_inputs/clinical_glioma.txt", header = TRUE, sep = "\t", check.names = FALSE)
d <- merge(data.frame(patient = rownames(risk), risk, check.names = FALSE),
           cli[, c("patient", "age", "sex", "grade", "idh")], by = "patient")
d$Grade <- factor(d$grade); d$IDH <- factor(d$idh); units(d$futime) <- "Year"
d <- d[complete.cases(d[, c("futime", "fustat", "riskScore", "age", "Grade", "IDH")]), ]
dd <- datadist(d); options(datadist = "dd")
f <- cph(Surv(futime, fustat) ~ riskScore + age + Grade + IDH, data = d, x = TRUE, y = TRUE, surv = TRUE)
say("n=%d，事件 %d；mean(lp)=%.3g（说明 lp 已中心化）", nrow(d), sum(d$fustat), mean(predict(f, type = "lp")))

TIMES <- c(1, 3, 5)
sv <- survest(f, newdata = d, times = TIMES)$surv
KM <- sapply(TIMES, function(u) summary(survfit(Surv(futime, fustat) ~ 1, data = d), times = u, extend = TRUE)$surv)
say("survest()（正确）             预测均值：%s", paste(sprintf("%.3f", colMeans(sv)), collapse = "/"))
for (cen in c(TRUE, FALSE)) {
  bh <- basehaz(f, centered = cen)
  H0 <- approx(bh$time, bh$hazard, xout = TIMES, rule = 2)$y
  S <- sapply(seq_along(TIMES), function(k) mean(exp(-exp(predict(f, type = "lp")) * H0[k])))
  say("basehaz(centered=%-5s)  预测均值：%s  与 survest 最大差 %.3f",
      cen, paste(sprintf("%.3f", S), collapse = "/"), max(abs(S - colMeans(sv))))
}
say("Kaplan-Meier 观测             ：%s", paste(sprintf("%.3f", KM), collapse = "/"))

kmc <- survfit(Surv(futime, 1 - fustat) ~ 1, data = d)
Gt <- function(t) { s <- summary(kmc, times = t, extend = TRUE)$surv
  ifelse(is.na(s) | s <= 0, min(s[!is.na(s) & s > 0], 1), s) }
brier <- function(S, t) {
  ev <- d$fustat == 1 & d$futime <= t; at <- d$futime > t
  w <- ifelse(ev, 1 / Gt(d$futime), 0)
  full <- mean(w * (0 - S)^2 + ifelse(at, 1 / Gt(t), 0) * (1 - S)^2)
  km <- summary(survfit(Surv(futime, fustat) ~ 1, data = d), times = t, extend = TRUE)$surv
  null <- mean(w * (0 - km)^2 + ifelse(at, 1 / Gt(t), 0) * (1 - km)^2)
  c(Brier = full, KMref = null, IPA = 1 - full / null)
}
say("\n%-6s %-28s %-28s", "Time", "survest（修正后）", "basehaz(centered=FALSE)（旧）")
bh_old <- basehaz(f, centered = FALSE); H0o <- approx(bh_old$time, bh_old$hazard, xout = TIMES, rule = 2)$y
S_old <- sapply(seq_along(TIMES), function(k) exp(-exp(predict(f, type = "lp")) * H0o[k]))
for (k in seq_along(TIMES)) {
  a <- brier(sv[, k], TIMES[k]); b <- brier(S_old[, k], TIMES[k])
  say("%-6d Brier %.3f / IPA %.3f        Brier %.3f / IPA %.3f", TIMES[k], a[1], a[3], b[1], b[3])
}
writeLines(c(sprintf("== 预测生存率量纲核对（%s）==", format(Sys.time(), "%Y-%m-%d %H:%M")), LOG),
           file.path(OUTD, "prediction_scale_check.txt"))
cat("\nwritten: 04_results/audit/prediction_scale_check.txt\nSTEP audit_prediction_scale DONE\n")
