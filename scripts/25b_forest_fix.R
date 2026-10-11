## 图5A / 图5A2 森林图重画（修复空白页）
## 病因一：IDH 的 95%CI 上限 12.12 超出 xlim=4，forestploter 直接画空
## 病因二：forestploter 1.1.4 在已打开的绘图设备内部构建对象会输出空页，必须先构建再打印
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(forestploter); library(grid)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT); source("scripts/plot_style.R")

LAB <- c(riskScore = "Risk score", age = "Age", Sex = "Sex (male vs female)",
         Grade = "WHO grade", IDH = "IDH status", SexMale = "Male (vs female)",
         GradeG3 = "WHO grade G3 (vs G2)", GradeG4 = "WHO grade G4 (vs G2)",
         IDHWT = "IDH wildtype (vs mutant)")

mk_tab <- function(tab) {
  lab <- ifelse(tab$Variable %in% names(LAB), LAB[tab$Variable], tab$Variable)
  data.frame(Variable = lab,
             `Hazard ratio (95% CI)` = sprintf("%.2f (%.2f-%.2f)", tab$HR, tab$L, tab$H),
             `P value` = ifelse(tab$p < 0.001, "<0.001", sprintf("%.3f", tab$p)),
             ` ` = paste(rep(" ", 22), collapse = ""), check.names = FALSE)
}
th <- forest_theme(base_size = 12, refline_gp = gpar(col = "grey45", lty = 2),
                   ci_col = "#0072B2", ci_lty = 1, ci_lwd = 1.5, ci_Theight = 0.22,
                   core = list(bg_params = list(fill = c("white", "grey97"))))

## 关键：先在无设备状态下构建 grob
u <- fread("04_results/25_26/independence_univariable_TCGA.txt")
m <- fread("04_results/25_26/independence_multivariable_TCGA.txt")
cat("单因素 HR 范围:", sprintf("%.2f-%.2f", min(u$L), max(u$H)),
    "| 多因素 HR 范围:", sprintf("%.2f-%.2f", min(m$L), max(m$H)), "\n")
pA <- forest(mk_tab(u), est = u$HR, lower = u$L, upper = u$H, sizes = 0.75, ci_column = 4,
             ref_line = 1, xlim = c(0, 13), ticks_at = c(0, 4, 8, 12), theme = th)
pB <- forest(mk_tab(m), est = m$HR, lower = m$L, upper = m$H, sizes = 0.75, ci_column = 4,
             ref_line = 1, xlim = c(0, 3), ticks_at = c(0, 1, 2, 3), theme = th)

draw <- function(p, file, w, h) {
  grDevices::cairo_pdf(file, width = w, height = h, family = "sans"); print(p); dev.off()
  grDevices::png(sub("\\.pdf$", ".png", file), width = w, height = h, units = "in", res = 200, type = "cairo")
  print(p); dev.off()
  cat("saved:", file, "\n")
}
draw(pA, "03_figures/Fig5A_independence_TCGA.pdf", 6.6, 3.2)
draw(pB, "03_figures/Fig5A2_independence_multi_TCGA.pdf", 6.6, 3.6)
cat("FOREST FIX DONE\n")
