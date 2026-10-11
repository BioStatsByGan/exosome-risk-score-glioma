## 表 1 基线临床特征（TCGA 训练 / CGGA-325 验证，按风险组）
## 表 2 多因素 Cox 模型（15 基因）
## 表 3 独立预后分析
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(survival)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
dir.create("稿/tables", recursive = TRUE, showWarnings = FALSE)

fmt_p <- function(p) ifelse(is.na(p), NA, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

## ---------- 表1 基线 ----------
risk <- read.table("04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t",
                   check.names = FALSE, row.names = 1)
cli <- read.table("02_inputs/clinical_glioma.txt", header = TRUE, sep = "\t", check.names = FALSE)
tcga <- merge(data.frame(id = rownames(risk), risk, check.names = FALSE),
              cli[, c("patient", "age", "sex", "grade", "idh", "mgmt")], by.x = "id", by.y = "patient")

cg <- read.table("04_results/18.CGGA/risk.CGGA325.txt", header = TRUE, sep = "\t", check.names = FALSE)
cgc <- fread("00_data/cgga/CGGA.mRNAseq_325_clinical.20200506.txt")
cgc <- as.data.frame(cgc); names(cgc) <- gsub(" ", ".", names(cgc))
cgc2 <- cgc[, c("CGGA_ID", "Age", "Gender", "Grade", "IDH_mutation_status", "MGMTp_methylation_status")]
names(cgc2) <- c("id", "age2", "sex2", "grade2", "idh2", "mgmt2")
cg <- merge(cg, cgc2, by = "id", all.x = TRUE)

baseline <- function(d, vars, label) {
  ## 生存分析集：缺随访时间或生存状态者不纳入（与 Cox/KM 的完整病例口径一致）
  d <- d[!is.na(d$futime) & !is.na(d$fustat), ]
  out <- data.frame(Characteristic = character(), High = character(), Low = character(), P = character())
  add <- function(name, hi, lo, p) {
    out <<- rbind(out, data.frame(Characteristic = name, High = hi, Low = lo, P = fmt_p(p)))
  }
  hi <- d[d$Risk == "high", ]; lo <- d[d$Risk == "low", ]
  add("No. of patients", as.character(nrow(hi)), as.character(nrow(lo)), NA)
  a_hi <- hi$age; a_lo <- lo$age
  wt <- suppressWarnings(wilcox.test(a_hi, a_lo))
  add("Age, median (IQR)", sprintf("%.0f (%.0f-%.0f)", median(a_hi, na.rm = TRUE),
                                   quantile(a_hi, .25, na.rm = TRUE), quantile(a_hi, .75, na.rm = TRUE)),
      sprintf("%.0f (%.0f-%.0f)", median(a_lo, na.rm = TRUE),
              quantile(a_lo, .25, na.rm = TRUE), quantile(a_lo, .75, na.rm = TRUE)), wt$p.value)
  cat_ <- function(name, v_hi, v_lo) {
    tb <- table(c(v_hi, v_lo), c(rep("high", length(v_hi)), rep("low", length(v_lo))))
    tb <- tb[rownames(tb) %in% unique(c(v_hi, v_lo)), , drop = FALSE]
    p <- tryCatch(suppressWarnings(chisq.test(tb)$p.value), error = function(e) NA)
    if (is.na(p)) p <- tryCatch(fisher.test(tb)$p.value, error = function(e) NA)
    for (lv in rownames(tb)) {
      if (lv %in% c("NA", "")) next
      add(sprintf("%s, n (%%)", lv),
          sprintf("%d (%.1f)", tb[lv, "high"], 100 * tb[lv, "high"] / length(v_hi)),
          sprintf("%d (%.1f)", tb[lv, "low"], 100 * tb[lv, "low"] / length(v_lo)), p)
    }
  }
  cat_(paste0(label, ": sex"), as.character(hi$sex), as.character(lo$sex))
  cat_(paste0(label, ": GRADE"), as.character(hi[[vars[1]]]), as.character(lo[[vars[1]]]))
  cat_(paste0(label, ": IDH"), as.character(hi[[vars[2]]]), as.character(lo[[vars[2]]]))
  if (length(vars) > 2 && vars[3] %in% names(d)) cat_(paste0(label, ": MGMT"), as.character(hi[[vars[3]]]), as.character(lo[[vars[3]]]))
  add("Dead, n (%)", sprintf("%d (%.1f)", sum(hi$fustat == 1, na.rm = TRUE), 100 * mean(hi$fustat == 1, na.rm = TRUE)),
      sprintf("%d (%.1f)", sum(lo$fustat == 1, na.rm = TRUE), 100 * mean(lo$fustat == 1, na.rm = TRUE)),
      suppressWarnings(chisq.test(table(c(hi$fustat, lo$fustat), c(rep(1, nrow(hi)), rep(0, nrow(lo)))))$p.value))
  add("Follow-up, median (IQR), years",
      sprintf("%.2f (%.2f-%.2f)", median(hi$futime, na.rm = TRUE), quantile(hi$futime, .25, na.rm = TRUE), quantile(hi$futime, .75, na.rm = TRUE)),
      sprintf("%.2f (%.2f-%.2f)", median(lo$futime, na.rm = TRUE), quantile(lo$futime, .25, na.rm = TRUE), quantile(lo$futime, .75, na.rm = TRUE)),
      suppressWarnings(wilcox.test(hi$futime, lo$futime)$p.value))
  out
}
t1a <- baseline(tcga, c("grade", "idh", "mgmt"), "TCGA")
t1b <- baseline(data.frame(cg, grade = cg$grade2, idh = cg$idh2, sex = cg$sex2, age = cg$age2,
                           mgmt = cg$mgmt2),
                c("grade", "idh", "mgmt"), "CGGA")
write.csv(t1a, "稿/tables/表1a_基线_TCGA.csv", row.names = FALSE)
write.csv(t1b, "稿/tables/表1b_基线_CGGA325.csv", row.names = FALSE)

## ---------- 表2 多因素 Cox（15 基因）----------
mc <- fread("04_results/20.SHAP/multiCox.txt")
t2 <- data.frame(Gene = mc$id,
                 Coefficient = round(mc$coef, 3),
                 HR = sprintf("%.3f", mc$HR),
                 `95% CI` = sprintf("%.3f-%.3f", mc$HR.95L, mc$HR.95H),
                 P = fmt_p(mc$pvalue), check.names = FALSE)
write.csv(t2, "稿/tables/表2_多因素Cox_15基因.csv", row.names = FALSE)

## ---------- 表3 独立预后分析 ----------
u <- fread("04_results/25_26/independence_univariable_TCGA.txt")
m <- fread("04_results/25_26/independence_multivariable_TCGA.txt")
g <- fread("04_results/25_26/independence_multivariable_CGGA325.txt")
mk <- function(t, model) data.frame(Model = model, Variable = t$Variable,
                                    HR = sprintf("%.2f", t$HR), `95% CI` = sprintf("%.2f-%.2f", t$L, t$H),
                                    P = fmt_p(t$p), check.names = FALSE)
t3 <- rbind(mk(u, "TCGA univariable"), mk(m, "TCGA multivariable"),
            mk(g[Model == "Grade"], "CGGA-325 multivariable (adjusted for WHO grade)"))
write.csv(t3, "稿/tables/表3_独立预后分析.csv", row.names = FALSE)

cat("表1a 行数:", nrow(t1a), "| 表1b 行数:", nrow(t1b), "| 表2 行数:", nrow(t2), "| 表3 行数:", nrow(t3), "\n")
cat("\n=== 表1a TCGA ===\n"); print(t1a, row.names = FALSE)
cat("\n=== 表1b CGGA ===\n"); print(t1b, row.names = FALSE)
cat("\n=== 表2 ===\n"); print(t2, row.names = FALSE)
cat("\n=== 表3 ===\n"); print(t3, row.names = FALSE)
cat("TABLES DONE\n")
