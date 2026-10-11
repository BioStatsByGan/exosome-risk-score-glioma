## 步骤 48：统一「固定系数 13 基因评分」口径（回应外部审查 P0 第 1 条）
## 目的：TCGA 中直接算出与 CGGA 完全一致的固定系数 13 基因评分，并与 15 基因模型、
##       TCGA 内重拟合的 13 基因模型三方对照；同时量化 TPM vs FPKM 的评分尺度差异。
## 输出：04_results/audit/gene13_fixed.txt
##       稿/tables/表S10_13基因固定系数对照.csv
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(survival)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
OUTD <- "04_results/audit"; dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
LOG <- character(); say <- function(...) { s <- sprintf(...); cat(s, "\n"); LOG <<- c(LOG, s); invisible(s) }

ci <- function(f) as.numeric(summary(f)$concordance[1])

mc <- fread("04_results/20.SHAP/multiCox.txt")
coefs <- setNames(mc$coef, mc$id)
say("模型基因 %d 个；CGGA 缺失的两基因：%s", length(coefs),
    paste(c("ADGRE5", "AQP1")[c("ADGRE5","AQP1") %in% names(coefs)], collapse = ", "))
KEEP13 <- setdiff(names(coefs), c("ADGRE5", "AQP1"))

## ---------- 1) TCGA ----------
rk <- read.table("04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t",
                 check.names = FALSE, row.names = 1)
E <- as.matrix(rk[, names(coefs)])
stopifnot(!any(is.na(E)))
S15 <- as.numeric(E %*% coefs)                      # 15 基因（原发布评分）
S13 <- as.numeric(E[, KEEP13] %*% coefs[KEEP13])    # 固定系数 13 基因（与 CGGA 完全同口径）
say("TCGA：重算 15 基因评分与文件内 riskScore 的最大绝对差 = %.3g",
    max(abs(S15 - rk$riskScore)))
d <- data.frame(futime = rk$futime, fustat = rk$fustat, S15 = S15, S13 = S13)

fit15 <- coxph(Surv(futime, fustat) ~ S15, data = d)
fit13 <- coxph(Surv(futime, fustat) ~ S13, data = d)
d13r <- data.frame(d, riskRefit = predict(coxph(Surv(futime, fustat) ~ ., data = data.frame(
  futime = d$futime, fustat = d$fustat, E[, KEEP13]))))

r15 <- c(C = ci(fit15), HR = unname(exp(coef(fit15))), L = exp(confint(fit15))[1], H = exp(confint(fit15))[2],
         P = summary(fit15)$coefficients[1, 5])
r13 <- c(C = ci(fit13), HR = unname(exp(coef(fit13))), L = exp(confint(fit13))[1], H = exp(confint(fit13))[2],
         P = summary(fit13)$coefficients[1, 5])
r13r <- c(C = ci(coxph(Surv(futime, fustat) ~ riskRefit, data = d13r)),
          HR = unname(exp(coef(coxph(Surv(futime, fustat) ~ riskRefit, data = d13r)))),
          L = exp(confint(coxph(Surv(futime, fustat) ~ riskRefit, data = d13r)))[1],
          H = exp(confint(coxph(Surv(futime, fustat) ~ riskRefit, data = d13r)))[2],
          P = summary(coxph(Surv(futime, fustat) ~ riskRefit, data = d13r))$coefficients[1, 5])
rho <- cor(d$S15, d$S13, method = "spearman")
say("TCGA 15 基因评分      ：C=%.3f，每单位 HR=%.2f (%.2f–%.2f)，P=%.2g", r15["C"], r15["HR"], r15["L"], r15["H"], r15["P"])
say("TCGA 固定系数 13 基因  ：C=%.3f，每单位 HR=%.2f (%.2f–%.2f)，P=%.2g（ΔC=%.3f）",
    r13["C"], r13["HR"], r13["L"], r13["H"], r13["P"], r13["C"] - r15["C"])
say("TCGA 重拟合 13 基因    ：C=%.3f（ΔC=%.3f）", r13r["C"], r13r["C"] - r15["C"])
say("15 基因与固定系数 13 基因评分的 Spearman ρ = %.3f", rho)

## ---------- 2) CGGA：确认外部队列用的就是固定系数 13 基因评分 ----------
apply_score <- function(nm) {
  ## 注意：不能写成 fread(sprintf(...))，data.table 在函数内重新求值 input 表达式会丢实参
  fp <- sprintf("00_data/cgga/CGGA.mRNAseq_%s.RSEM-genes.20200506.txt", nm)
  ex <- fread(fp, showProgress = FALSE)
  gn <- ex[[1]]; m <- as.matrix(ex[, -1]); rownames(m) <- gn
  use <- intersect(KEEP13, rownames(m)); m <- log2(m[use, , drop = FALSE] + 1)
  sc <- colSums(coefs[use] * m)
  list(score = sc, n_gene = length(use))
}
cg <- list()
for (nm in c("325", "693")) {
  a <- apply_score(nm)
  rr <- read.table(sprintf("04_results/18.CGGA/risk.CGGA%s.txt", nm), header = TRUE, sep = "\t", check.names = FALSE)
  same <- intersect(names(a$score), rr$id)
  rho_cg <- cor(as.numeric(a$score[same]), rr$riskScore[match(same, rr$id)], method = "spearman")
  dd <- data.frame(futime = rr$futime[match(same, rr$id)], fustat = rr$fustat[match(same, rr$id)],
                   S13 = as.numeric(a$score[same]))
  dd <- dd[!is.na(dd$futime) & !is.na(dd$fustat), ]
  f <- coxph(Surv(futime, fustat) ~ S13, data = dd)
  cg[[nm]] <- dd
  say("CGGA-%s：用 %d 个基因重算固定系数评分，n=%d，C=%.3f，每单位 HR=%.2f (%.2f–%.2f)；与文件内 riskScore 的 ρ=%.3f",
      nm, a$n_gene, nrow(dd), ci(f), exp(coef(f)), exp(confint(f))[1], exp(confint(f))[2],
      rho_cg)
}

## ---------- 3) 尺度：TCGA(log2 TPM+1) vs CGGA(log2 FPKM+1) ----------
gmean <- function(E, genes) colMeans(E[, genes])
say("\n各队列评分分布（固定系数 13 基因口径）：")
dist_line <- function(x, label) say("  %-10s 中位 %.3f，IQR %.3f–%.3f，SD %.3f，范围 %.3f–%.3f",
                                    label, median(x), quantile(x, .25), quantile(x, .75), sd(x), min(x), max(x))
dist_line(S13, "TCGA")
dist_line(cg[["325"]]$S13, "CGGA-325")
dist_line(cg[["693"]]$S13, "CGGA-693")
cuts <- quantile(S13, c(.25, .5, .75))
say("TCGA 评分四分位（%s）", paste(sprintf("%.3f", cuts), collapse = ", "))
for (nm in c("325", "693")) {
  x <- cg[[nm]]$S13
  say("  CGGA-%s 高于 TCGA 中位切点的比例 = %.1f%%；高于 TCGA 上四分位 = %.1f%%",
      nm, 100 * mean(x > cuts[2]), 100 * mean(x > cuts[3]))
}
## 基因层面尺度（TCGA 用 risk.TCGA.txt 内的表达值，即 log2(TPM+1)）
gm_t <- colMeans(E[, KEEP13])
ex325 <- fread("00_data/cgga/CGGA.mRNAseq_325.RSEM-genes.20200506.txt", showProgress = FALSE)
g325 <- ex325[[1]]; m325 <- as.matrix(ex325[, -1]); rownames(m325) <- g325
gm_c <- rowMeans(log2(m325[KEEP13, , drop = FALSE] + 1))
say("\n基因层面平均表达（log2(TPM+1) vs log2(FPKM+1)）：")
say("  TCGA 均值 %.2f（范围 %.2f–%.2f）；CGGA-325 均值 %.2f（范围 %.2f–%.2f）；两者相关 ρ=%.3f",
    mean(gm_t), min(gm_t), max(gm_t), mean(gm_c), min(gm_c), max(gm_c),
    cor(gm_t, gm_c[match(names(gm_t), names(gm_c))], method = "spearman", use = "complete.obs"))

## ---------- 4) 落盘 ----------
tb <- data.frame(
  Cohort = c("TCGA (15-gene, published)", "TCGA (fixed-coefficient 13-gene)",
             "TCGA (13-gene refitted)", "CGGA-325 (fixed-coefficient 13-gene)", "CGGA-693 (fixed-coefficient 13-gene)"),
  n = c(nrow(d), nrow(d), nrow(d), nrow(cg[["325"]]), nrow(cg[["693"]])),
  C_index = round(c(r15["C"], r13["C"], r13r["C"], ci(coxph(Surv(futime, fustat) ~ S13, data = cg[["325"]])),
                    ci(coxph(Surv(futime, fustat) ~ S13, data = cg[["693"]]))), 3),
  HR_per_unit = round(c(r15["HR"], r13["HR"], r13r["HR"],
                        exp(coef(coxph(Surv(futime, fustat) ~ S13, data = cg[["325"]]))),
                        exp(coef(coxph(Surv(futime, fustat) ~ S13, data = cg[["693"]])))), 3),
  Score_median = round(c(median(S15), median(S13), NA, median(cg[["325"]]$S13), median(cg[["693"]]$S13)), 3),
  Score_SD = round(c(sd(S15), sd(S13), NA, sd(cg[["325"]]$S13), sd(cg[["693"]]$S13)), 3))
fwrite(tb, "稿/tables/表S10_13基因固定系数对照.csv", sep = ",", quote = FALSE)
say("\n表 S10 已写入 稿/tables/表S10_13基因固定系数对照.csv")
print(tb, row.names = FALSE)

writeLines(c(sprintf("== 固定系数 13 基因评分口径核对（%s）==", format(Sys.time(), "%Y-%m-%d %H:%M")), LOG),
           file.path(OUTD, "gene13_fixed.txt"))
cat("\nwritten: 04_results/audit/gene13_fixed.txt\nSTEP48 DONE\n")
