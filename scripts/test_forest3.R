options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(forestploter); library(grid)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
u <- fread(file.path(ROOT, "04_results/25_26/independence_univariable_TCGA.txt"))
t2 <- data.frame(Variable = u$Variable,
                 `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", u$HR, u$L, u$H),
                 `P` = ifelse(u$p < 0.001, "<0.001", sprintf("%.3f", u$p)),
                 ` ` = paste(rep(" ", 20), collapse = ""), check.names = FALSE)

themeA <- forest_theme(base_size = 10, ci_col = "#0072B2", ci_lwd = 1.5, ci_Theight = 0.22)
themeB <- forest_theme(base_size = 10, refline_col = "grey45", ci_col = "#0072B2",
                       ci_lty = 1, ci_lwd = 1.5, ci_Theight = 0.22,
                       row_fill = c("white", "grey96"),
                       core = list(bg_params = list(fill = c("white", "grey97"), alpha = 1)))

mk <- function(th) forest(t2, est = u$HR, lower = u$L, upper = u$H, sizes = 0.9, ci_column = 4,
                          ref_line = 1, xlim = c(0, 13), ticks_at = c(0, 3, 6, 9, 12), theme = th)
dev_list <- list(cairo_pdf = grDevices::cairo_pdf, pdf = grDevices::pdf)
for (dn in names(dev_list)) {
  for (tn in c("A", "B")) {
    f <- file.path(ROOT, sprintf("04_results/dev_%s_theme%s.pdf", dn, tn))
    d <- dev_list[[dn]]
    d(f, width = 6.4, height = 3.4)
    r <- try(print(mk(if (tn == "A") themeA else themeB)), silent = TRUE)
    dev.off()
    cat(sprintf("%-10s theme%s -> %s\n", dn, tn, if (inherits(r, "try-error")) attr(r, "condition")$message else "OK"))
  }
}
