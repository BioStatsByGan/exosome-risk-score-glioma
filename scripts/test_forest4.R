options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(forestploter); library(grid)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
u <- fread(file.path(ROOT, "04_results/25_26/independence_univariable_TCGA.txt"))
t2 <- data.frame(Variable = u$Variable,
                 `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", u$HR, u$L, u$H),
                 `P` = ifelse(u$p < 0.001, "<0.001", sprintf("%.3f", u$p)),
                 ` ` = paste(rep(" ", 20), collapse = ""), check.names = FALSE)

combos <- expand.grid(xlim = c(FALSE, TRUE), ticks = c(FALSE, TRUE), sizes = c(FALSE, TRUE))
for (i in seq_len(nrow(combos))) {
  cc <- combos[i, ]
  args <- list(data = t2, est = u$HR, lower = u$L, upper = u$H, ci_column = 4, ref_line = 1)
  if (cc$xlim) args$xlim <- c(0, 13)
  if (cc$ticks) args$ticks_at <- c(0, 3, 6, 9, 12)
  if (cc$sizes) args$sizes <- 0.9
  f <- file.path(ROOT, sprintf("04_results/grid_%s%s%s.pdf", cc$xlim, cc$ticks, cc$sizes))
  pdf(f, width = 6.4, height = 3.4)
  r <- try(print(do.call(forest, args)), silent = TRUE)
  dev.off()
  cat(sprintf("xlim=%s ticks=%s sizes=%s -> %s\n", cc$xlim, cc$ticks, cc$sizes,
              if (inherits(r, "try-error")) attr(r, "condition")$message else "OK"))
}
