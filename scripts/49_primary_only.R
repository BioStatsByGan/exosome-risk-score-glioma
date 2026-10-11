## 步骤 49：CGGA 原发肿瘤限定分析（回应外部审查 P0 第 2 条）
## 目的：CGGA 两队列混合原发/复发/继发肿瘤且无取材-诊断时间信息，故限定原发肿瘤重复主要验证分析，
##       检验主结论是否依赖复发/继发样本；并核对临床表中是否存在可用的时间起点字段。
## 输出：04_results/audit/primary_only.txt
##       稿/tables/表S11_原发肿瘤限定分析.csv
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(survival); library(survminer)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
OUTD <- "04_results/audit"; dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
LOG <- character(); say <- function(...) { s <- sprintf(...); cat(s, "\n"); LOG <<- c(LOG, s); invisible(s) }
ci <- function(f) as.numeric(summary(f)$concordance[1])

## ---------- 0) 先看临床表里到底有没有「时间起点」字段 ----------
for (nm in c("325", "693")) {
  fp <- sprintf("00_data/cgga/CGGA.mRNAseq_%s_clinical.20200506.txt", nm)
  h <- names(fread(fp, nrows = 0))
  say("CGGA-%s 临床表字段：%s", nm, paste(h, collapse = " | "))
  say("   含日期/时间起点类字段：%s", paste(grep("date|Date|time|Time|year|Year|diagnos|surg", h, value = TRUE), collapse = ", "))
}

## ---------- 1) 逐队列：全部 vs 仅原发 ----------
res <- list()
for (nm in c("325", "693")) {
  rk <- read.table(sprintf("04_results/18.CGGA/risk.CGGA%s.txt", nm), header = TRUE, sep = "\t", check.names = FALSE)
  cl <- fread(sprintf("00_data/cgga/CGGA.mRNAseq_%s_clinical.20200506.txt", nm), data.table = FALSE)
  names(cl) <- gsub(" .*", "", names(cl))
  d <- merge(rk, cl[, c("CGGA_ID", "PRS_type")], by.x = "id", by.y = "CGGA_ID")
  d <- d[!is.na(d$futime) & !is.na(d$fustat), ]   # 与主分析同口径：313 / 657（含 4 例类型未记录）
  d$Grade <- factor(d$Grade, levels = c("WHO II", "WHO III", "WHO IV"))
  say("\n== CGGA-%s：可分析 %d 例（%s）", nm, nrow(d),
      paste(sprintf("%s=%d", names(table(d$PRS_type)), as.integer(table(d$PRS_type))), collapse = " "))

  one <- function(dd, label) {
    if (nrow(dd) < 30 || sum(dd$fustat) < 15) { say("   %s：样本过少，跳过", label); return(NULL) }
    f <- coxph(Surv(futime, fustat) ~ riskScore, data = dd)
    sg <- survdiff(Surv(futime, fustat) ~ Risk, data = dd)
    ddg <- dd[!is.na(dd$Grade), ]
    fg <- if (nrow(ddg) > 30) coxph(Surv(futime, fustat) ~ riskScore + Grade, data = ddg) else NULL
    fmt_p <- function(p) if (is.na(p)) "NA" else if (p < 2.2e-16) "<2.2e-16" else format(signif(p, 3))
    say("   %-12s n=%3d，事件 %3d | C=%.3f | 每单位 HR=%.2f (%.2f–%.2f) P=%.2g | log-rank P=%.2g%s",
        label, nrow(dd), sum(dd$fustat), ci(f), exp(coef(f))[1], exp(confint(f))[1], exp(confint(f))[2],
        summary(f)$coefficients[1, 5], 1 - pchisq(sg$chisq, 1),
        if (!is.null(fg)) sprintf(" | 校正分级后 HR=%.2f (%.2f–%.2f) P=%s",
                                  exp(coef(fg))[1], exp(confint(fg))[1, 1], exp(confint(fg))[1, 2],
                                  fmt_p(summary(fg)$coefficients[1, 5])) else "")
    data.frame(Cohort = sprintf("CGGA-%s %s", nm, label), n = nrow(dd), Events = sum(dd$fustat),
               C_index = round(ci(f), 3), HR = round(exp(coef(f))[1], 2),
               HR_low = round(exp(confint(f))[1], 2), HR_high = round(exp(confint(f))[2], 2),
               P = signif(summary(f)$coefficients[1, 5], 3),
               logrank_P = format(1 - pchisq(sg$chisq, 1), scientific = TRUE, digits = 3),
               HR_grade_adj = if (!is.null(fg)) round(exp(coef(fg))[1], 2) else NA)
  }
  res[[paste0(nm, "_all")]] <- one(d, "全部")
  res[[paste0(nm, "_pri")]] <- one(d[which(d$PRS_type == "Primary"), ], "仅原发")
  if (nm == "325") res[[paste0(nm, "_rec")]] <- one(d[which(d$PRS_type == "Recurrent"), ], "仅复发")
  if (nm == "325") res[[paste0(nm, "_sec")]] <- one(d[which(d$PRS_type == "Secondary"), ], "仅继发")
}

tb <- do.call(rbind, res)
tb <- tb[!duplicated(tb$Cohort), ]
fwrite(tb, "稿/tables/表S11_原发肿瘤限定分析.csv", sep = ",", quote = FALSE)
say("\n表 S11 已写入 稿/tables/表S11_原发肿瘤限定分析.csv")
print(tb, row.names = FALSE)

writeLines(c(sprintf("== CGGA 原发肿瘤限定分析（%s）==", format(Sys.time(), "%Y-%m-%d %H:%M")), LOG),
           file.path(OUTD, "primary_only.txt"))
cat("\nwritten: 04_results/audit/primary_only.txt\nSTEP49 DONE\n")
