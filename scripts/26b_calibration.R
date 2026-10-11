## 步骤 26 补充：校准曲线（时间单位修正为年）+ 完整打印独立预后模型
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(survival); library(rms); library(data.table)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
## 控制台输出同时写入日志文件，供审计回源核对（多因素模型、C 指数、校准斜率）
sink("04_results/25_26/model_and_calibration_stats.txt", split = TRUE)

risk <- read.table("04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t",
                   check.names = FALSE, row.names = 1)
cli <- read.table("02_inputs/clinical_glioma.txt", header = TRUE, sep = "\t", check.names = FALSE)
d <- merge(data.frame(patient = rownames(risk), risk, check.names = FALSE),
           cli[, c("patient", "age", "sex", "grade", "idh")], by = "patient")
d$Grade <- factor(d$grade, levels = c("G2", "G3", "G4")); d$IDH <- factor(d$idh); d$Sex <- factor(d$sex)
units(d$futime) <- "Year"          # 关键：把生存时间单位标成年，calibrate 的 u 才是年

cat("===== TCGA 多因素 Cox (riskScore + age + Sex + Grade + IDH) =====\n")
m <- summary(coxph(Surv(futime, fustat) ~ riskScore + age + Sex + Grade + IDH, data = d))
print(data.frame(Variable = rownames(m$coefficients),
                 HR = sprintf("%.2f (%.2f-%.2f)", m$conf.int[, 1], m$conf.int[, 3], m$conf.int[, 4]),
                 p = signif(m$coefficients[, 5], 3), row.names = NULL), row.names = FALSE)
cat(sprintf("C-index = %.3f\n", m$concordance[1]))

cg <- read.table("04_results/18.CGGA/risk.CGGA325.txt", header = TRUE, sep = "\t", check.names = FALSE)
## 生存分析集：仅剔除缺随访时间/生存状态者；本模型只含 riskScore + Grade，不额外要求 IDH 标签完整
## （与 04_results/25_26/independence_multivariable_CGGA325.txt / 表 3 同口径，保证两处数值一致）
cg <- cg[!is.na(cg$futime) & !is.na(cg$fustat), ]
## 分级缺失写作空串；限定三个 WHO 级别并把参照组固定为 WHO II（与 25_26 / 表 3 完全同口径）
cg <- cg[cg$Grade %in% c("WHO II", "WHO III", "WHO IV"), ]
cg$Grade <- factor(cg$Grade, levels = c("WHO II", "WHO III", "WHO IV"))
cg$IDH <- factor(cg$IDH, levels = c("Mutant", "Wildtype"))
cat(sprintf("CGGA-325 有 WHO 分级且可分析 n = %d（参照组 WHO II）\n", nrow(cg)))
cat("\n===== CGGA-325 多因素 Cox (riskScore + WHO grade) =====\n")
m2 <- summary(coxph(Surv(futime, fustat) ~ riskScore + Grade, data = cg))
print(data.frame(Variable = rownames(m2$coefficients),
                 HR = sprintf("%.2f (%.2f-%.2f)", m2$conf.int[, 1], m2$conf.int[, 3], m2$conf.int[, 4]),
                 p = signif(m2$coefficients[, 5], 3), row.names = NULL), row.names = FALSE)
cat(sprintf("C-index = %.3f\n", m2$concordance[1]))
cat("\n===== CGGA-325 单因素 IDH 与 分级 的分离性检查 =====\n")
cgL <- cg[cg$Grade %in% c("WHO II", "WHO III", "WHO IV") & cg$IDH %in% c("Mutant", "Wildtype"), ]
cat(sprintf("（分级与 IDH 标签齐全的子集 n = %d）\n", nrow(cgL)))
print(table(cgL$Grade, cgL$IDH))
s_idh <- summary(coxph(Surv(futime, fustat) ~ IDH, data = cg))
cat(sprintf("仅 IDH 的单因素模型 HR=%.3f (%.3f-%.3f) p=%.3g\n",
            s_idh$conf.int[1, 1], s_idh$conf.int[1, 3], s_idh$conf.int[1, 4], s_idh$coefficients[1, 5]))

## ---------- 校准曲线（cph + bootstrap KM）----------
dd <- datadist(d); options(datadist = "dd")
f <- cph(Surv(futime, fustat) ~ riskScore + age + Grade + IDH, data = d, x = TRUE, y = TRUE, surv = TRUE)
set.seed(7)
res <- list()
for (u in c(1, 3, 5)) {
  cal <- calibrate(f, cmethod = "KM", method = "boot", u = u, m = 100, B = 200)
  res[[as.character(u)]] <- cal
  cc <- cal[is.finite(cal[, "mean.predicted"]) & is.finite(cal[, "KM"]), , drop = FALSE]
  x <- as.numeric(cc[, "mean.predicted"]); y <- as.numeric(cc[, "KM"])
  cat(sprintf("%d 年校准: 预测区间 %.2f-%.2f | 观测区间 %.2f-%.2f | 校准斜率(线性拟合)=%.2f\n",
              u, min(x), max(x), min(y), max(y), coef(lm(y ~ x))[2]))
}
saveRDS(res, "04_results/25_26/calibration.rds")

draw_cal <- function() {
  par(mfrow = c(1, 3), mar = c(4, 4, 2.5, 1))
  for (u in c(1, 3, 5)) {
    cc <- res[[as.character(u)]]; cc <- cc[is.finite(cc[, "mean.predicted"]) & is.finite(cc[, "KM"]), , drop = FALSE]
    x <- as.numeric(cc[, "mean.predicted"]); y <- as.numeric(cc[, "KM"])
    plot(x, y, xlim = c(0.2, 1), ylim = c(0.2, 1), pch = 16, cex = 0.7, col = "#0072B2",
         xlab = sprintf("Predicted %d-year survival", u), ylab = "Observed survival",
         axes = FALSE, main = NULL, type = "b", lwd = 0.8)
    abline(0, 1, lty = 2, col = "grey50"); box(); axis(1); axis(2)
    mtext(sprintf("%d-year", u), side = 3, line = 0.4, cex = 0.85)
  }
}
cairo_pdf("03_figures/Fig5E_calibration.pdf", width = 8.4, height = 3.0); draw_cal(); dev.off()
png("03_figures/Fig5E_calibration.png", width = 8.4, height = 3.0, units = "in", res = 200, type = "cairo"); draw_cal(); dev.off()
cat("saved: Fig5E_calibration.pdf\n")
sink()
cat("written: 04_results/25_26/model_and_calibration_stats.txt\n")
cat("STEP26b DONE\n")
