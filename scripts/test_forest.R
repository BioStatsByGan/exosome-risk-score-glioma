options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(forestploter); library(grid); library(ggplot2)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
u <- fread(file.path(ROOT, "04_results/25_26/independence_univariable_TCGA.txt"))
print(u)
t2 <- data.frame(Variable = u$Variable,
                 `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", u$HR, u$L, u$H),
                 `P` = ifelse(u$p < 0.001, "<0.001", sprintf("%.3f", u$p)),
                 ` ` = paste(rep(" ", 20), collapse = ""), check.names = FALSE)
cat("列数:", ncol(t2), "行数:", nrow(t2), "\n")

try_pdf <- function(expr, file, w = 6.4, h = 3.4) {
  pdf(file, width = w, height = h)
  r <- try(expr, silent = TRUE)
  dev.off()
  cat(file, "->", if (inherits(r, "try-error")) paste("ERROR:", attr(r, "condition")$message) else "OK", "\n")
}
try_pdf(print(forest(t2, est = u$HR, lower = u$L, upper = u$H, ci_column = 4, ref_line = 1)), file.path(ROOT, "04_results/t1.pdf"))
try_pdf(print(forest(t2, est = u$HR, lower = u$L, upper = u$H, ci_column = 4, ref_line = 1, xlim = c(0, 13), ticks_at = c(0, 3, 6, 9, 12))), file.path(ROOT, "04_results/t2.pdf"))
tm <- forest_theme(base_size = 10, ci_col = "#0072B2", ci_lwd = 1.5, ci_Theight = 0.22,
                   core = list(bg_params = list(fill = c("white", "grey97"))))
try_pdf(print(forest(t2, est = u$HR, lower = u$L, upper = u$H, ci_column = 4, ref_line = 1, xlim = c(0, 13), ticks_at = c(0, 3, 6, 9, 12), theme = tm)), file.path(ROOT, "04_results/t3.pdf"))
