options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(forestploter); library(grid); library(ggplot2)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
u <- fread(file.path(ROOT, "04_results/25_26/independence_univariable_TCGA.txt"))
t2 <- data.frame(Variable = u$Variable,
                 `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", u$HR, u$L, u$H),
                 `P` = ifelse(u$p < 0.001, "<0.001", sprintf("%.3f", u$p)),
                 ` ` = paste(rep(" ", 20), collapse = ""), check.names = FALSE)
cat("forestploter 版本:", as.character(packageVersion("forestploter")), "\n")
p <- forest(t2, est = u$HR, lower = u$L, upper = u$H, ci_column = 4, ref_line = 1, xlim = c(0, 13))
cat("class(p):", paste(class(p), collapse = ", "), "\n")
for (m in c("print", "plot", "grid.draw", "draw")) {
  f <- file.path(ROOT, sprintf("04_results/m_%s.pdf", m))
  pdf(f, width = 6.4, height = 3.4)
  r <- try(switch(m, print = print(p), plot = plot(p), grid.draw = grid.draw(p), draw = grid.draw(ggplotify::as.grob(p))), silent = TRUE)
  dev.off()
  cat(sprintf("%-10s -> %s\n", m, if (inherits(r, "try-error")) paste("ERROR:", attr(r, "condition")$message) else "OK"))
}
