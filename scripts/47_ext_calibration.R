## 步骤 47：外部队列（CGGA-325 / CGGA-693）的校准与预测误差
## 回应外部意见第 4 条子项：校准（bootstrap 斜率）与 Brier/IPA 原先只有 TCGA，外部队列需要同样的输出。
## 模型与 TCGA 完全同口径（risk score + age + WHO grade + IDH，即图 5E/5F 的列线图模型），
## 指标实现也沿用 TCGA 的同一套代码（rms::calibrate 的 KM 分箱 + Graf 型 IPW Brier）。
## 输出：04_results/audit/external_calibration.txt（日志）
##       稿/tables/表S9_外部队列校准与预测误差.csv
##       03_figures/FigS6_ext_calibration.{pdf,png}
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(survival); library(rms)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
OUTD <- "04_results/audit"; dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
LOG <- character(); say <- function(...) { s <- sprintf(...); cat(s, "\n"); LOG <<- c(LOG, s); invisible(s) }

## ---------- 数据准备 ----------
## TCGA（与 26b / audit_b_group.R 同口径：grade、idh 直接取因子）
rk <- read.table("04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
cl <- read.table("02_inputs/clinical_glioma.txt", header = TRUE, sep = "\t", check.names = FALSE)
dt <- merge(data.frame(id = rownames(rk), futime = rk$futime, fustat = rk$fustat, riskScore = rk$riskScore),
            cl[, c("patient", "age", "grade", "idh")], by.x = "id", by.y = "patient")
tcga <- data.frame(Cohort = "TCGA", id = dt$id, futime = dt$futime, fustat = dt$fustat,
                   riskScore = dt$riskScore, Age = as.numeric(dt$age),
                   Grade = factor(dt$grade), IDH = factor(dt$idh), stringsAsFactors = FALSE)

## CGGA-325 / 693（生存分析集 = 有随访时间与生存状态者，与 18_cgga_validate.R 一致）
load_cgga <- function(nm) {
  cg <- read.table(sprintf("04_results/18.CGGA/risk.CGGA%s.txt", nm), header = TRUE, sep = "\t", check.names = FALSE)
  cgc <- fread(sprintf("00_data/cgga/CGGA.mRNAseq_%s_clinical.20200506.txt", nm), data.table = FALSE)
  names(cgc) <- gsub(" .*", "", names(cgc))
  d <- merge(cg[, c("id", "futime", "fustat", "riskScore")],
             cgc[, c("CGGA_ID", "Age", "Grade", "IDH_mutation_status")], by.x = "id", by.y = "CGGA_ID")
  d <- d[!is.na(d$futime) & !is.na(d$fustat), ]
  d$Age <- as.numeric(d$Age)
  d$Grade <- factor(d$Grade, levels = c("WHO II", "WHO III", "WHO IV"))
  d$IDH <- factor(ifelse(grepl("^mut", d$IDH_mutation_status, ignore.case = TRUE), "Mutant", "Wildtype"),
                  levels = c("Mutant", "Wildtype"))
  data.frame(Cohort = sprintf("CGGA-%s", nm), id = d$id, futime = d$futime, fustat = d$fustat,
             riskScore = d$riskScore, Age = d$Age, Grade = d$Grade, IDH = d$IDH, stringsAsFactors = FALSE)
}

VARS <- c("futime", "fustat", "riskScore", "Age", "Grade", "IDH")
TIMES <- c(1, 3, 5)
CAL <- list()          # 保存各队列的 calibrate 输出，供绘图复用（避免重复 bootstrap）

## ---------- 指标（与 TCGA 同一实现） ----------
metrics <- function(d, label) {
  d <- d[complete.cases(d[, VARS]), ]
  units(d$futime) <- "Year"
  dd <<- datadist(d); options(datadist = "dd")   # <<- 必须写到全局环境，rms 才能按名找到 datadist
  f <- cph(Surv(futime, fustat) ~ riskScore + Age + Grade + IDH, data = d, x = TRUE, y = TRUE, surv = TRUE)
  ## C 指数一律走 coxph 的 concordance（survival::concordance(Surv~x) 会把 x 当生存时间，方向反）
  C <- as.numeric(summary(coxph(Surv(futime, fustat) ~ riskScore + Age + Grade + IDH, data = d))$concordance[1])
  say("== %s：完整病例 n=%d，事件 %d，列线图模型 C 指数=%.3f ==", label, nrow(d), sum(d$fustat), C)

  ## 1) 校准斜率（rms::calibrate，KM 分箱 + bootstrap，与 TCGA 同设置 m=100, B=200）
  set.seed(7); cal_list <- list(); slopes <- setNames(numeric(3), as.character(TIMES))
  for (u in TIMES) {
    cal <- calibrate(f, cmethod = "KM", method = "boot", u = u, m = 100, B = 200)
    cal_list[[as.character(u)]] <- cal
    cc <- cal[is.finite(cal[, "mean.predicted"]) & is.finite(cal[, "KM"]), , drop = FALSE]
    x <- as.numeric(cc[, "mean.predicted"]); y <- as.numeric(cc[, "KM"])
    slopes[as.character(u)] <- coef(lm(y ~ x))[2]
    say("   %d 年校准：斜率 %.2f（预测 %.2f–%.2f，观测 %.2f–%.2f，分箱 %d）",
        u, slopes[as.character(u)], min(x), max(x), min(y), max(y), nrow(cc))
  }
  CAL[[label]] <<- cal_list

  ## 2) 时间依赖 Brier 分数与 IPA（Graf 型 IPW，手算，公式与 audit_b_group.R 相同）
  kmc <- survfit(Surv(futime, 1 - fustat) ~ 1, data = d)   # 删失分布
  Gt <- function(t) {
    s <- summary(kmc, times = t, extend = TRUE)$surv
    ifelse(is.na(s) | s <= 0, min(s[!is.na(s) & s > 0], 1), s)
  }
  bh <- basehaz(f, centered = FALSE)
  H0 <- approx(bh$time, bh$hazard, xout = TIMES, rule = 2)$y
  sp <- exp(-outer(exp(predict(f, type = "lp")), H0))
  bs <- do.call(rbind, lapply(seq_along(TIMES), function(k) {
    t <- TIMES[k]; Sh <- sp[, k]
    ev <- d$fustat == 1 & d$futime <= t
    at <- d$futime > t
    w <- ifelse(ev, 1 / Gt(d$futime), 0)
    full <- mean(w * (0 - Sh)^2 + ifelse(at, 1 / Gt(t), 0) * (1 - Sh)^2)
    km <- summary(survfit(Surv(futime, fustat) ~ 1, data = d), times = t, extend = TRUE)$surv
    null <- mean(w * (0 - km)^2 + ifelse(at, 1 / Gt(t), 0) * (1 - km)^2)
    data.frame(Time = t, Brier_model = full, Brier_KM = null, IPA = 1 - full / null)
  }))
  ## 注意：R 4.4 的 sprintf 在「格式串含非 ASCII 字符 + 传入向量参数」时会误报 too few arguments，
  ## 所以下面这条格式串一律用 ASCII，中文标签放在外层 say() 里（外层只传标量，安全）。
  say("   Brier model / KM / IPA (1,3,5 y): %s",
      paste(sprintf("%.3f / %.3f / %.3f", bs$Brier_model, bs$Brier_KM, bs$IPA), collapse = " ;  "))

  data.frame(Cohort = label, n = nrow(d), Events = sum(d$fustat), C_index = C,
             Time = bs$Time, Calibration_slope = as.numeric(slopes[as.character(bs$Time)]),
             Brier_model = bs$Brier_model, Brier_KM = bs$Brier_KM, IPA = bs$IPA)
}

DATA <- list("TCGA" = tcga, "CGGA-325" = load_cgga("325"), "CGGA-693" = load_cgga("693"))
res <- do.call(rbind, lapply(names(DATA), function(k) metrics(DATA[[k]], k)))
res[, 5:9] <- lapply(res[, 5:9], function(x) round(as.numeric(x), 3))
fwrite(res, "稿/tables/表S9_外部队列校准与预测误差.csv", sep = ",", quote = FALSE)
say("\n汇总（写入 稿/tables/表S9_外部队列校准与预测误差.csv）：")
print(res, row.names = FALSE, digits = 3)

## ---------- 图 S6：两队列的校准曲线（复用上面的 calibrate 输出） ----------
draw_ext_cal <- function(nm) {
  for (u in TIMES) {
    cal <- CAL[[nm]][[as.character(u)]]
    cc <- cal[is.finite(cal[, "mean.predicted"]) & is.finite(cal[, "KM"]), , drop = FALSE]
    x <- as.numeric(cc[, "mean.predicted"]); y <- as.numeric(cc[, "KM"])
    plot(x, y, xlim = c(0, 1), ylim = c(0, 1), pch = 16, cex = 0.7, col = "#0072B2",
         xlab = sprintf("Predicted %d-year survival", u), ylab = "Observed survival",
         axes = FALSE, type = "b", lwd = 0.8, main = NULL)
    abline(0, 1, lty = 2, col = "grey50"); box(); axis(1); axis(2)
    mtext(sprintf("%s (%d-year, slope %.2f)", nm, u, coef(lm(y ~ x))[2]), side = 3, line = 0.4, cex = 0.8)
  }
}
cairo_pdf("03_figures/FigS6_ext_calibration.pdf", width = 8.4, height = 6.0)
par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1))
draw_ext_cal("CGGA-325"); draw_ext_cal("CGGA-693"); dev.off()
png("03_figures/FigS6_ext_calibration.png", width = 8.4, height = 6.0, units = "in", res = 200, type = "cairo")
par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1)); draw_ext_cal("CGGA-325"); draw_ext_cal("CGGA-693"); dev.off()
say("\n图 S6 已保存：03_figures/FigS6_ext_calibration.pdf/.png")

writeLines(c(sprintf("== 外部队列校准与预测误差（%s）==", format(Sys.time(), "%Y-%m-%d %H:%M")), LOG),
           file.path(OUTD, "external_calibration.txt"))
cat("\nwritten: 04_results/audit/external_calibration.txt\n")
cat("STEP47 DONE\n")
