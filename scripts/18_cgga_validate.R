## 步骤 18 + 外部验证：15 基因模型在 CGGA 325 / 693 队列上验证
## 输入：04_results/20.SHAP/multiCox.txt（系数）、00_data/cgga/*（CGGA 表达与临床）
## 输出：04_results/18.CGGA/{risk.CGGA325.txt, risk.CGGA693.txt, validate_stats.txt}
##       03_figures/Fig3{A,B,C,D}_*.pdf（KM + 时间 ROC）
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({
  library(data.table); library(survival); library(survminer)
  library(timeROC); library(ggplot2); library(patchwork)
})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
setwd(ROOT); dir.create("04_results/18.CGGA", recursive = TRUE, showWarnings = FALSE)
source("scripts/plot_style.R")   # 提供 roc_step / auc_time / OI / theme_pub

## 关键输出同时落盘，供审计回源核对（原本只打印到控制台）
LOG <- character()
say <- function(...) { s <- sprintf(...); cat(s, "\n"); LOG <<- c(LOG, s); invisible(s) }

mc <- fread("04_results/20.SHAP/multiCox.txt")
genes <- mc$id; coefs <- mc$coef; names(coefs) <- genes
cat("模型基因:", length(genes), "\n")

## ---- 通用：套用系数算分 ----
apply_model <- function(expFile, cliFile, cohortLabel) {
  ex <- fread(expFile, showProgress = FALSE)
  gn <- ex[[1]]; m <- as.matrix(ex[, -1]); rownames(m) <- gn; rm(ex)
  miss <- setdiff(genes, rownames(m))
  if (length(miss)) cat("  缺失基因:", paste(miss, collapse = ", "), "\n")
  use <- intersect(genes, rownames(m))
  m <- m[use, , drop = FALSE]
  # 同一基因多探针取平均，log2(FPKM+1)
  m <- log2(m + 1)
  score <- colSums(coefs[use] * m)
  cli <- as.data.frame(fread(cliFile))
  cn <- colnames(cli)
  osCol <- cn[grep("^OS$", cn)]; cenCol <- cn[grep("^Censor", cn)]
  cli$futime <- cli[[osCol]] / 365; cli$fustat <- cli[[cenCol]]
  cli$sample <- cli$CGGA_ID
  same <- intersect(colnames(m), cli$sample)
  # 表达列名去下划线匹配
  if (!length(same)) {
    colnames(m) <- gsub("_", "-", colnames(m)); cli$sample <- gsub("_", "-", cli$sample)
    same <- intersect(colnames(m), cli$sample)
  }
  cli <- cli[match(same, cli$sample), ]
  s <- score[same]
  data.frame(id = same, futime = cli$futime, fustat = cli$fustat,
             riskScore = as.numeric(s), Risk = ifelse(s > median(s), "high", "low"),
             Grade = cli[[grep("^Grade", colnames(cli))[1]]],
             IDH = cli[[grep("IDH", colnames(cli))[1]]], stringsAsFactors = FALSE)
}

report <- function(d, label) {
  fit <- survfit(Surv(futime, fustat) ~ Risk, data = d)
  sd <- survdiff(Surv(futime, fustat) ~ Risk, data = d)
  cox <- coxph(Surv(futime, fustat) ~ riskScore, data = d)
  ci <- summary(cox)$concordance
  r <- timeROC(T = d$futime, delta = d$fustat, marker = d$riskScore, cause = 1,
               times = c(1, 3, 5), iid = TRUE)
  cin <- confint(r, parm = NULL, n.sim = 500); cin$CI <- cin$CI / 100
  say("%s: 入组 n=%d / 可分析 n=%d | HR/单位log2=%.2f (%.2f-%.2f) p=%.2g | C-index=%.3f | logrank p=%.2g",
      label, nrow(d), sum(!is.na(d$futime) & !is.na(d$fustat)),
      exp(coef(cox)), exp(confint(cox))[1], exp(confint(cox))[2],
      summary(cox)$coefficients[, "Pr(>|z|)"], ci[1], 1 - pchisq(sd$chisq, 1))
  say("   1/3/5 年 AUC: %s",
      paste(sprintf("%.3f(%.2f-%.2f)", r$AUC, cin$CI[, 1], cin$CI[, 2]), collapse = " "))
  list(d = d, fit = fit, p = 1 - pchisq(sd$chisq, 1), cox = cox, auc = r$AUC)
}

cp <- "00_data/cgga/"
d325 <- apply_model(paste0(cp, "CGGA.mRNAseq_325.RSEM-genes.20200506.txt"),
                    paste0(cp, "CGGA.mRNAseq_325_clinical.20200506.txt"), "CGGA325")
d693 <- apply_model(paste0(cp, "CGGA.mRNAseq_693.RSEM-genes.20200506.txt"),
                    paste0(cp, "CGGA.mRNAseq_693_clinical.20200506.txt"), "CGGA693")
fwrite(cbind(id = d325$id, d325[, -1]), "04_results/18.CGGA/risk.CGGA325.txt", sep = "\t", quote = FALSE)
fwrite(cbind(id = d693$id, d693[, -1]), "04_results/18.CGGA/risk.CGGA693.txt", sep = "\t", quote = FALSE)

cat("\n===== 外部验证结果 =====\n")
r325 <- report(d325, "CGGA-325 (主验证)")
r693 <- report(d693, "CGGA-693 (敏感性验证)")

## 分层：按 WHO 分级校正的多因素 Cox
for (nm in c("325", "693")) {
  d <- if (nm == "325") d325 else d693
  d$Grade <- as.factor(d$Grade); d$IDH <- as.factor(d$IDH)
  f <- tryCatch(coxph(Surv(futime, fustat) ~ riskScore + Grade + IDH, data = d), error = function(e) NULL)
  if (!is.null(f)) {
    s <- summary(f)$coefficients
    say("CGGA-%s 校正 Grade+IDH 后 riskScore: HR=%.2f (%.2f-%.2f) p=%.3g | C-index=%.3f", nm,
        exp(s["riskScore", "coef"]), exp(confint(f))["riskScore", 1], exp(confint(f))["riskScore", 2],
        s["riskScore", "Pr(>|z|)"], summary(f)$concordance[1])
  }
}

## ---- 图：KM + 时间 ROC（4 panel）----
km_plot <- function(d, label) {
  kv <- survfit(Surv(futime, fustat) ~ Risk, data = d)
  g <- ggsurvplot(kv, data = d, pval = TRUE, risk.table = TRUE, conf.int = FALSE,
                  palette = c("#D55E00", "#0072B2"), xlab = "Time (years)",
                  ylab = "Overall survival probability", legend.title = "",
                  legend.labs = c("High risk", "Low risk"), tables.height = 0.22,
                  ggtheme = theme_pub(), risk.table.y.text = FALSE, censor.size = 1.6,
                  break.time.by = 2, xlim = c(0, 10))$plot + labs(title = NULL)
  g <- g + theme(legend.position = c(0.78, 0.85))
  g
}
pA <- km_plot(d325, "CGGA-325"); pB <- km_plot(d693, "CGGA-693")
pC <- roc_step(d325, times = c(1, 3, 5), marker_col = "riskScore")
pD <- roc_step(d693, times = c(1, 3, 5), marker_col = "riskScore")
fig <- (pA | pC) / (pB | pD) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "bold"))
ggsave("03_figures/Fig3_external_validation.pdf", fig, width = 11, height = 9, device = cairo_pdf)
writeLines(c(sprintf("== CGGA 外部验证关键输出（%s，R %s，timeROC CI n.sim=500）==",
                     format(Sys.time(), "%Y-%m-%d %H:%M"), getRversion()), LOG),
           "04_results/18.CGGA/validate_stats.txt")
cat("saved: Fig3_external_validation.pdf\n")
cat("written: 04_results/18.CGGA/validate_stats.txt\n")
cat("STEP18 DONE\n")
