## 步骤 25-26：独立预后分析 + 临床相关性 + 列线图与校准（TCGA 训练 + CGGA 验证）
## 输出：03_figures/Fig5A_independence_TCGA.pdf, Fig5B_risk_by_clinical.pdf,
##       03_figures/Fig5C_clinical_heatmap.pdf, Fig5D_nomogram.pdf, Fig5E_calibration.pdf
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({
  library(data.table); library(ggplot2); library(survival); library(survminer)
  library(forestploter); library(patchwork); library(grid); library(ComplexHeatmap)
  library(circlize); library(rms)
})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
setwd(ROOT); source("scripts/plot_style.R")
dir.create("04_results/25_26", recursive = TRUE, showWarnings = FALSE)

risk <- read.table("04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t",
                   check.names = FALSE, row.names = 1)
cli <- read.table("02_inputs/clinical_glioma.txt", header = TRUE, sep = "\t", check.names = FALSE)
cli$grade_num <- as.numeric(gsub("G", "", cli$grade))
d <- merge(data.frame(patient = rownames(risk), risk, check.names = FALSE),
           cli[, c("patient", "age", "sex", "grade", "grade_num", "idh", "mgmt", "cohort")],
           by = "patient")
d$Grade <- factor(d$grade, levels = c("G2", "G3", "G4"))
d$IDH <- factor(d$idh); d$Sex <- factor(d$sex); d$MGMT <- factor(d$mgmt)
cat("纳入临床信息样本:", nrow(d), "\n")

## ---------- Fig5A 独立预后：单因素 + 多因素 ----------
vars <- c("riskScore", "age", "Sex", "Grade", "IDH")
uni <- do.call(rbind, lapply(vars, function(v) {
  f <- as.formula(sprintf("Surv(futime, fustat) ~ %s", v))
  s <- summary(coxph(f, data = d))
  data.frame(Variable = v, HR = s$conf.int[1, 1], L = s$conf.int[1, 3],
             H = s$conf.int[1, 4], p = s$coefficients[1, 5], row.names = NULL)
}))
mul <- summary(coxph(Surv(futime, fustat) ~ riskScore + age + Sex + Grade + IDH, data = d))
mul <- data.frame(Variable = rownames(mul$coefficients), HR = mul$conf.int[, 1],
                  L = mul$conf.int[, 3], H = mul$conf.int[, 4], p = mul$coefficients[, 5], row.names = NULL)
fwrite(uni, "04_results/25_26/independence_univariable_TCGA.txt", sep = "\t", quote = FALSE)
fwrite(mul, "04_results/25_26/independence_multivariable_TCGA.txt", sep = "\t", quote = FALSE)
cat("\n多因素 Cox (TCGA):\n"); print(mul, row.names = FALSE, digits = 3)

fp <- function(tab, title) {
  t2 <- data.frame(Variable = tab$Variable,
                   `Hazard ratio (95% CI)` = sprintf("%.2f (%.2f-%.2f)", tab$HR, tab$L, tab$H),
                   `P value` = ifelse(tab$p < 0.001, "<0.001", sprintf("%.3f", tab$p)),
                   ` ` = paste(rep(" ", 20), collapse = ""), check.names = FALSE)
  forest(t2, est = tab$HR, lower = tab$L, upper = tab$H, ci_column = 4, ref_line = 1,
         sizes = 0.9, xlim = c(0, 4), ticks_at = c(0.5, 1, 2, 3, 4),
         theme = forest_theme(base_size = 10, ci_col = "#0072B2", ci_lwd = 1.5,
                              ci_Theight = 0.22, core = list(bg_params = list(fill = c("white", "grey97")))))
}
save_pdf(fp(uni, "Univariable"), "03_figures/Fig5A_independence_TCGA.pdf", 6.2, 3.4)
save_pdf(fp(mul, "Multivariable"), "03_figures/Fig5A2_independence_multi_TCGA.pdf", 6.2, 3.4)

## CGGA 验证队列的独立预后（Grade 与 IDH 在低级别胶质瘤中高度共线，分开建模再合并）
cg <- read.table("04_results/18.CGGA/risk.CGGA325.txt", header = TRUE, sep = "\t", check.names = FALSE)
## 注意：缺失值在 risk 文件里写作空串 ""（不是 NA），若不显式剔除会变成一个假参照类别；
## 每个模型各自用其完整病例，并把参照组固定为 WHO II / IDH 突变型
cg$Grade <- factor(cg$Grade, levels = c("WHO II", "WHO III", "WHO IV"))
cg$IDH <- factor(cg$IDH, levels = c("Mutant", "Wildtype"))
mk_multi <- function(form, label) {
  vars <- all.vars(as.formula(form)); d <- cg
  if ("Grade" %in% vars) d <- d[d$Grade %in% c("WHO II", "WHO III", "WHO IV"), ]
  if ("IDH" %in% vars) d <- d[d$IDH %in% c("Mutant", "Wildtype"), ]
  d$Grade <- droplevels(d$Grade); d$IDH <- droplevels(d$IDH)
  cat(sprintf("   [%s] n=%d（参照组 %s）\n", label, nrow(d),
              paste(c(if ("Grade" %in% vars) levels(d$Grade)[1], if ("IDH" %in% vars) levels(d$IDH)[1]), collapse = " / ")))
  m <- summary(coxph(as.formula(form), data = d))
  t <- data.frame(Variable = rownames(m$coefficients), HR = m$conf.int[, 1],
                  L = m$conf.int[, 3], H = m$conf.int[, 4], p = m$coefficients[, 5], row.names = NULL)
  cat("\n", label, "\n"); print(t, row.names = FALSE, digits = 3)
  t
}
t_g <- mk_multi("Surv(futime, fustat) ~ riskScore + Grade", "CGGA-325 校正 WHO 分级")
t_i <- mk_multi("Surv(futime, fustat) ~ riskScore + IDH", "CGGA-325 校正 IDH 状态")
t_gi <- mk_multi("Surv(futime, fustat) ~ riskScore + Grade + IDH", "CGGA-325 校正 分级+IDH（两者共线，仅供参考）")
fwrite(rbind(cbind(Model = "Grade", t_g), cbind(Model = "IDH", t_i), cbind(Model = "Grade+IDH", t_gi)),
       "04_results/25_26/independence_multivariable_CGGA325.txt", sep = "\t", quote = FALSE)

## ---------- Fig5B 风险评分按临床分层 ----------
mk_box <- function(gene_col, lab, pal) {
  ggplot(d[!is.na(d[[gene_col]]), ], aes(.data[[gene_col]], riskScore, fill = .data[[gene_col]])) +
    geom_boxplot(width = 0.55, outlier.size = 0.5, linewidth = 0.35, colour = "grey25") +
    geom_jitter(width = 0.16, size = 0.35, alpha = 0.35, colour = "grey35") +
    scale_fill_manual(values = pal, guide = "none") +
    labs(x = lab, y = "Risk score") +
    stat_compare_means(method = "kruskal.test", label = "p.format", size = 2.8, vjust = -0.5)
}
suppressMessages(library(ggpubr))
pb1 <- mk_box("Grade", "WHO grade", c("#0072B2", "#E69F00", "#D55E00"))
pb2 <- mk_box("IDH", "IDH status", c("#009E73", "#CC79A7"))
pb3 <- mk_box("Sex", "Sex", c("#56B4E9", "#E69F00"))
pb4 <- mk_box("MGMT", "MGMT promoter", c("#0072B2", "#D55E00"))
save_pdf((pb1 | pb2) / (pb3 | pb4), "03_figures/Fig5B_risk_by_clinical.pdf", 8.4, 7.2)

## ---------- Fig5C 临床特征热图 ----------
hd <- d[order(d$riskScore), ]
mat <- t(scale(as.matrix(data.frame(
  `Risk score` = hd$riskScore, Age = hd$age, `WHO grade` = hd$grade_num,
  `IDH mutation` = as.numeric(hd$IDH == "Mutant"), `MGMT methylated` = as.numeric(hd$MGMT == "Methylated"),
  check.names = FALSE))))
colnames(mat) <- NULL
ha <- HeatmapAnnotation(`Risk group` = hd$Risk, col = list(`Risk group` = c(low = "#0072B2", high = "#D55E00")),
                        annotation_name_gp = gpar(fontsize = 8))
ht <- Heatmap(mat, name = "Z score", top_annotation = ha, cluster_columns = FALSE, show_column_names = FALSE,
              row_names_side = "left", row_names_gp = gpar(fontsize = 8.5),
              col = colorRamp2(c(-2, 0, 2), c("#0072B2", "white", "#D55E00")),
              column_title = NULL, heatmap_legend_param = list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)))
cairo_pdf("03_figures/Fig5C_clinical_heatmap.pdf", width = 6.4, height = 2.6); draw(ht); dev.off()
png("03_figures/Fig5C_clinical_heatmap.png", width = 6.4, height = 2.6, units = "in", res = 200, type = "cairo"); draw(ht); dev.off()
cat("saved: Fig5C_clinical_heatmap.pdf\n")

## ---------- Fig5D/E 列线图 + 校准 ----------
dd <- datadist(d); options(datadist = "dd")
f <- cph(Surv(futime, fustat) ~ riskScore + age + Grade + IDH, data = d, x = TRUE, y = TRUE, surv = TRUE)
surv <- Survival(f)
nom <- nomogram(f, fun = list(function(x) surv(1, x), function(x) surv(3, x), function(x) surv(5, x)),
                funlabel = c("1-year survival", "3-year survival", "5-year survival"),
                lp = FALSE, fun.at = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1))
cairo_pdf("03_figures/Fig5D_nomogram.pdf", width = 7.6, height = 6.2); plot(nom, xfrac = 0.32, cex.axis = 0.7, cex.var = 0.85); dev.off()
png("03_figures/Fig5D_nomogram.png", width = 7.6, height = 6.2, units = "in", res = 200, type = "cairo"); plot(nom, xfrac = 0.32, cex.axis = 0.7, cex.var = 0.85); dev.off()
cat("saved: Fig5D_nomogram.pdf\n")

## ---------- Fig5E 校准曲线 ----------
## 注：本脚本原先在这里重复做了一次 calibrate(u = c(1,3,5))，但 d$futime 的时间单位未标成年，
## 该调用必然报错（"the condition has length > 1"）且其输出被 26b_calibration.R 完全覆盖，
## 故此处不再重复，校准曲线与校准斜率统一由 scripts/26b_calibration.R 产出
## （该脚本设 units(d$futime) <- "Year"，并给出 1/3/5 年斜率 1.09 / 1.80 / 1.53）。
cat("Fig5E_calibration 由 scripts/26b_calibration.R 产出，此处跳过\n")
cat("STEP25-26 DONE\n")
