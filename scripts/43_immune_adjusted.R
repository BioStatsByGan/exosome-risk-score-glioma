## 步骤 43：22 类 CIBERSORT 细胞统一按「未校正 / 校正分级 / 校正分级+IDH」三种口径检验
##         并加一个「免疫含量校正」的敏感性分析（用 29 项 ssGSEA 免疫特征均值作免疫含量代理）
## 输出：04_results/audit/immune_all_cells_adjusted.txt / _log.txt
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
OUTD <- "04_results/audit"; dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
LOG <- character(); say <- function(...) { s <- sprintf(...); cat(s, "\n"); LOG <<- c(LOG, s); invisible(s) }

imm <- read.table("04_results/27.CIBERSORT/CIBERSORT-Results.txt", header = TRUE, sep = "\t",
                  check.names = FALSE, row.names = 1)
rownames(imm) <- sub("^(.*?)-(.*?)-(.*?)-.*", "\\1-\\2-\\3", rownames(imm))
cells <- setdiff(colnames(imm), c("P-value", "Correlation", "RMSE"))
rk <- read.table("04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
cli <- read.table("02_inputs/clinical_glioma.txt", header = TRUE, sep = "\t", check.names = FALSE)
ss <- fread("04_results/29.immFunction/ssgsea_scores.txt"); setnames(ss, 1, "Feature")
immune_content <- data.table(patient = names(ss)[-1], imm_content = colMeans(as.matrix(ss[, -1])))

d <- data.frame(patient = rownames(imm), imm[, cells], check.names = FALSE)
d <- merge(d, data.frame(patient = rownames(rk), Risk = rk$Risk), by = "patient")
d <- merge(d, cli[, c("patient", "grade", "idh")], by = "patient")
d <- merge(d, immune_content, by = "patient")
d$Grade <- factor(d$grade); d$IDH <- factor(d$idh)
say("纳入样本: %d（高危 %d / 低危 %d）", nrow(d), sum(d$Risk == "high"), sum(d$Risk == "low"))

one <- function(cc) {
  f <- function(rhs) summary(lm(as.formula(sprintf("`%s` ~ Risk + %s", cc, rhs)), data = d))$coefficients
  s1 <- f("1"); s2 <- f("Grade"); s3 <- f("Grade + IDH"); s4 <- f("Grade + IDH + imm_content")
  data.frame(Cell = cc, beta_unadj = s1["Risklow", 1], p_unadj = s1["Risklow", 4],
             beta_grade = s2["Risklow", 1], p_grade = s2["Risklow", 4],
             beta_grade_idh = s3["Risklow", 1], p_grade_idh = s3["Risklow", 4],
             beta_plus_immune = s4["Risklow", 1], p_plus_immune = s4["Risklow", 4])
}
res <- do.call(rbind, lapply(cells, one))
res$FDR_unadj <- p.adjust(res$p_unadj, "BH"); res$FDR_grade_idh <- p.adjust(res$p_grade_idh, "BH")
res <- res[order(res$p_unadj), ]
fwrite(res, file.path(OUTD, "immune_all_cells_adjusted.txt"), sep = "\t", quote = FALSE)
say("\n== 22 类 CIBERSORT 细胞三种口径（按未校正 P 排序，前 12）==")
print(head(res[, c("Cell","beta_unadj","p_unadj","p_grade","p_grade_idh","p_plus_immune","FDR_grade_idh")], 12), digits = 3, row.names = FALSE)
say("\n三种口径下 P<0.05 的细胞数：未校正 %d | 校正分级 %d | 校正分级+IDH %d | 再加免疫含量 %d",
    sum(res$p_unadj < 0.05), sum(res$p_grade < 0.05), sum(res$p_grade_idh < 0.05), sum(res$p_plus_immune < 0.05))
writeLines(c(sprintf("== 免疫细胞三口径校正（%s）==", format(Sys.time(), "%Y-%m-%d %H:%M")), LOG),
           file.path(OUTD, "immune_all_cells_adjusted_log.txt"))
cat("written:", file.path(OUTD, "immune_all_cells_adjusted.txt"), "\n")