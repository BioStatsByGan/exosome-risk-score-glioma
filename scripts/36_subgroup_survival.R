## 步骤 36：分层（亚组）验证 —— 连续评分法
## 设计说明：
##   本评分与 WHO 分级 / IDH 状态高度共线（R²=0.822），因此按这些变量分层后，
##   中位数切分的两组会严重失衡（如 TCGA IDH 野生型分层内 347 高 vs 23 低），
##   二分类 KM 在分层内不可靠。故分层内主分析采用连续评分（每单位 HR + C 指数），
##   分层内 KM 采用三分位（每组例数均衡）。
## 输出：03_figures/Fig6A_subgroup_forest.pdf, 03_figures/Fig6B_km_strata.pdf
##       04_results/36.subgroup/subgroup_continuous_TCGA.txt / _CGGA325.txt
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({
  library(data.table); library(ggplot2); library(survival); library(survminer)
  library(patchwork); library(forestploter); library(grid)
})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
setwd(ROOT); source("scripts/plot_style.R")
OUTD <- "04_results/36.subgroup"; dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)

## ================= 数据 =================
cl <- fread("02_inputs/clinical_glioma.txt")
rk <- fread("04_results/20.SHAP/risk.TCGA.txt")
tc <- merge(rk[, .(id, futime, fustat, riskScore, Risk)],
            cl[, .(patient, age, sex, grade, idh, mgmt)],
            by.x = "id", by.y = "patient")
tc <- tc[!is.na(futime) & !is.na(fustat)]
tc[, Risk := factor(Risk, levels = c("low", "high"))]
tc[, `:=`(Age = ifelse(age > 50, "Age > 50", "Age <= 50"),
          # 注意：idh/mgmt 存在 NA，必须先判空，否则 NA 会被归入野生型（曾经的 bug）
          IDH = ifelse(is.na(idh) | idh == "", NA_character_,
                       ifelse(grepl("mut", idh, ignore.case = TRUE), "IDH-mutant", "IDH-wildtype")),
          MGMT = ifelse(is.na(mgmt) | mgmt == "", NA_character_,
                        ifelse(grepl("unm", mgmt, ignore.case = TRUE), "MGMT unmethylated", "MGMT methylated")),
          Grade = paste0("WHO grade ", sub("^G", "", grade)), Sex = sex)]

cgcl <- fread("00_data/cgga/CGGA.mRNAseq_325_clinical.20200506.txt")
setnames(cgcl, gsub(" .*", "", names(cgcl)))
cgcl <- cgcl[, .(CGGA_ID, Gender, Age, Radio_status, Chemo_status, IDH_mutation_status,
                 `1p19q_codeletion_status`, MGMTp_methylation_status)]
cg <- fread("04_results/18.CGGA/risk.CGGA325.txt")
cg <- merge(cg[, .(id, futime, fustat, riskScore, Risk, Grade, IDH)], cgcl,
            by.x = "id", by.y = "CGGA_ID")
cg <- cg[!is.na(futime) & !is.na(fustat)]
cg[, Risk := factor(Risk, levels = c("low", "high"))]
cg[, `:=`(Age = ifelse(Age > 50, "Age > 50", "Age <= 50"),
          IDH = ifelse(is.na(IDH_mutation_status) | IDH_mutation_status == "", NA_character_,
                       ifelse(grepl("^mut", IDH_mutation_status, ignore.case = TRUE), "IDH-mutant", "IDH-wildtype")),
          MGMTp = ifelse(is.na(MGMTp_methylation_status) | MGMTp_methylation_status == "", NA_character_,
                         ifelse(grepl("^un", MGMTp_methylation_status, ignore.case = TRUE),
                                "MGMT promoter unmethylated", "MGMT promoter methylated")),
          Codel = ifelse(is.na(`1p19q_codeletion_status`) | `1p19q_codeletion_status` == "",
                         NA_character_,
                         ifelse(grepl("^codel", `1p19q_codeletion_status`, ignore.case = TRUE),
                                "1p/19q codeleted", "1p/19q intact")),
          Radio = ifelse(Radio_status == 1, "Radiotherapy treated", "Radiotherapy untreated"),
          Chemo = ifelse(Chemo_status == 1, "Chemotherapy treated", "Chemotherapy untreated"),
          Grade = as.character(Grade))]
cat("TCGA:", nrow(tc), "例 events", sum(tc$fustat), "| CGGA-325:", nrow(cg), "例 events", sum(cg$fustat), "\n")

## ================= 分层内统计 =================
strata_stats <- function(d, vars, cohort) {
  rbindlist(lapply(vars, function(v) {
    lv <- sort(unique(d[[v]])); lv <- lv[!is.na(lv) & !grepl("NA", as.character(lv))]
    rbindlist(lapply(lv, function(g) {
      s <- d[get(v) == g]
      if (nrow(s) < 40 || sum(s$fustat) < 15) return(NULL)
      fit <- coxph(Surv(futime, fustat) ~ riskScore, data = s)
      ci <- summary(fit)$conf.int
      cc <- summary(coxph(Surv(futime, fustat) ~ riskScore, data = s))$concordance
      n_hi <- sum(s$Risk == "high"); n_lo <- sum(s$Risk == "low")
      b_hr <- b_l <- b_h <- b_p <- NA_real_
      if (min(n_hi, n_lo) >= 20) {
        bf <- coxph(Surv(futime, fustat) ~ Risk, data = s); bc <- summary(bf)$conf.int
        b_hr <- bc[1, 1]; b_l <- bc[1, 3]; b_h <- bc[1, 4]; b_p <- summary(bf)$coefficients[1, 5]
      }
      data.table(Cohort = cohort, Stratum = v, Level = as.character(g), n = nrow(s),
                 Events = sum(s$fustat), HR = ci[1, 1], L = ci[1, 3], H = ci[1, 4],
                 P = summary(fit)$coefficients[1, 5], Cindex = as.numeric(cc["C"]),
                 n_high = n_hi, n_low = n_lo, HR_binary = b_hr, Lb = b_l, Hb = b_h, Pb = b_p)
    }))
  }))
}
overall_row <- function(d, cohort) {
  fit <- coxph(Surv(futime, fustat) ~ riskScore, data = d); ci <- summary(fit)$conf.int
  cc <- summary(fit)$concordance
  data.table(Cohort = cohort, Stratum = "Overall", Level = "All patients", n = nrow(d),
             Events = sum(d$fustat), HR = ci[1, 1], L = ci[1, 3], H = ci[1, 4],
             P = summary(fit)$coefficients[1, 5], Cindex = as.numeric(cc["C"]),
             n_high = sum(d$Risk == "high"), n_low = sum(d$Risk == "low"),
             HR_binary = NA_real_, Lb = NA_real_, Hb = NA_real_, Pb = NA_real_)
}
st_tc <- rbind(overall_row(tc, "TCGA"), strata_stats(tc, c("IDH", "MGMT", "Grade", "Age", "Sex"), "TCGA"))
st_cg <- rbind(overall_row(cg, "CGGA-325"),
               strata_stats(cg, c("IDH", "MGMTp", "Grade", "Age", "Codel", "Radio", "Chemo"), "CGGA-325"))
fmt <- function(x) { x <- copy(x); x[, `:=`(HR = round(HR, 2), L = round(L, 2), H = round(H, 2),
                                            Cindex = round(Cindex, 3), P = signif(P, 3))]; x }
st_tc <- fmt(st_tc); st_cg <- fmt(st_cg)
fwrite(st_tc, file.path(OUTD, "subgroup_continuous_TCGA.txt"), sep = "\t")
fwrite(st_cg, file.path(OUTD, "subgroup_continuous_CGGA325.txt"), sep = "\t")
cat("\n== TCGA ==\n"); print(st_tc[, .(Stratum, Level, n, Events, HR, L, H, P, Cindex, n_high, n_low)])
cat("\n== CGGA-325 ==\n"); print(st_cg[, .(Stratum, Level, n, Events, HR, L, H, P, Cindex, n_high, n_low)])

## ================= Fig6A 森林图（连续评分）=================
fp <- function(st) {
  d <- as.data.frame(st)
  d$` ` <- paste(rep(" ", 22), collapse = " ")
  d$`n (events)` <- sprintf("%d (%d)", d$n, d$Events)
  d$`HR per unit (95% CI)` <- sprintf("%.2f (%.2f-%.2f)", d$HR, d$L, d$H)
  d$`C-index` <- sprintf("%.3f", d$Cindex)
  d$`P value` <- ifelse(d$P < 0.001, "<0.001", sprintf("%.3f", d$P))
  d$Label <- paste0(d$Level)
  tm <- forest_theme(base_size = 8.5, ci_pch = 16, ci_col = OI[1], ci_lty = 1, ci_lwd = 1.4,
                     refline_gp = gpar(lty = "dashed", col = "grey40"),
                     core = list(bg_params = list(fill = c("white", "grey95"))))
  forest(d[, c("Label", "n (events)", " ", "HR per unit (95% CI)", "C-index", "P value")],
         est = d$HR, lower = d$L, upper = d$H, ci_column = 3, ref_line = 1,
         xlim = c(0.5, max(d$H) * 1.08), ticks_at = c(1, 2, 3),
         theme = tm,
         footnote = "Hazard ratio per unit of the risk score within each stratum (continuous score).")
}
g6a <- wrap_elements(fp(st_tc)) / wrap_elements(fp(st_cg))
save_pdf(g6a, "03_figures/Fig6A_subgroup_forest.pdf", 9.4, 7.8)   # 加高：亚层增到 16 行后坐标轴才不会被挤出页面

## ================= Fig6B 分层内三分位 KM =================
km_tertile <- function(d, title) {
  cat("  panel:", title, "n =", nrow(d), "\n")
  if (nrow(d) < 20) return(ggplot() + theme_void())
  d <- copy(d)
  q <- quantile(d$riskScore, c(1/3, 2/3))
  d[, Grp := factor(ifelse(riskScore <= q[1], "Low", ifelse(riskScore <= q[2], "Mid", "High")),
                    levels = c("Low", "Mid", "High"))]
  f <- Surv(futime, fustat) ~ Grp
  s <- survfit(f, data = d); sd <- survdiff(f, data = d)
  pv <- 1 - pchisq(sd$chisq, length(sd$n) - 1L)
  df <- data.frame(time = s$time, surv = s$surv,
                   grp = factor(rep(names(s$strata), times = as.integer(s$strata)),
                                levels = c("Grp=Low", "Grp=Mid", "Grp=High"),
                                labels = c("Low", "Mid", "High")))
  ggplot(df, aes(time, surv, colour = grp)) +
    geom_step(linewidth = 0.6) +
    labs(x = "Time (years)", y = "Overall survival probability",
         subtitle = sprintf("%s  (n=%d, deaths=%d)\nTertiles within stratum, %s",
                            title, nrow(d), sum(d$fustat), fmt_p_lab(pv, "log-rank p"))) +
    scale_colour_manual(values = c(Low = OI[1], Mid = OI[7], High = OI[2]), name = "") +
    theme_pub() +
    theme(legend.position = c(0.78, 0.9), legend.key.size = unit(0.32, "cm"),
          legend.text = element_text(size = 6.8), plot.subtitle = element_text(size = 7, lineheight = 1.05))
}
panels <- list(
  km_tertile(tc[IDH == "IDH-wildtype"], "TCGA IDH-wildtype"),
  km_tertile(tc[IDH == "IDH-mutant"], "TCGA IDH-mutant"),
  km_tertile(tc[Grade == "WHO grade 3"], "TCGA WHO grade 3"),
  km_tertile(cg[IDH == "IDH-wildtype"], "CGGA-325 IDH-wildtype"),
  km_tertile(cg[IDH == "IDH-mutant"], "CGGA-325 IDH-mutant"),
  km_tertile(cg[Grade == "WHO IV"], "CGGA-325 WHO grade IV"),
  km_tertile(cg[Codel == "1p/19q intact"], "CGGA-325 1p/19q intact"),
  km_tertile(cg[Radio == "Radiotherapy treated"], "CGGA-325 radiotherapy treated"))
save_pdf(wrap_plots(panels, ncol = 4), "03_figures/Fig6B_km_strata.pdf", 13.0, 6.4)
cat("\n完成：Fig6A_subgroup_forest.pdf / Fig6B_km_strata.pdf\n")
