options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(forestploter); library(grid)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
u <- fread(file.path(ROOT, "04_results/25_26/independence_univariable_TCGA.txt"))
t2 <- data.frame(Variable = u$Variable,
                 `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", u$HR, u$L, u$H),
                 `P` = ifelse(u$p < 0.001, "<0.001", sprintf("%.3f", u$p)),
                 ` ` = paste(rep(" ", 20), collapse = ""), check.names = FALSE)

## 配方 A：先在无设备状态下构建 grob，再开设备打印
pA <- forest(t2, est = u$HR, lower = u$L, upper = u$H, ci_column = 4, ref_line = 1, xlim = c(0, 13))
pdf(file.path(ROOT, "04_results/recipeA.pdf"), width = 6.4, height = 3.4); print(pA); dev.off()

## 配方 B：设备内构建并打印
pdf(file.path(ROOT, "04_results/recipeB.pdf"), width = 6.4, height = 3.4)
print(forest(t2, est = u$HR, lower = u$L, upper = u$H, ci_column = 4, ref_line = 1, xlim = c(0, 13)))
dev.off()

## 配方 C：cairo_pdf 设备 + 预先构建
pC <- forest(t2, est = u$HR, lower = u$L, upper = u$H, ci_column = 4, ref_line = 1, xlim = c(0, 13),
             theme = forest_theme(base_size = 10, ci_col = "#0072B2", ci_lwd = 1.5, ci_Theight = 0.22))
grDevices::cairo_pdf(file.path(ROOT, "04_results/recipeC.pdf"), width = 6.4, height = 3.4, family = "sans")
print(pC); dev.off()

## 配方 D：ggsave（forestploter 官方推荐之一，需要 ggplot 化）
suppressMessages(library(ggplot2))
ok <- requireNamespace("ggplotify", quietly = TRUE)
cat("ggplotify 可用:", ok, "\n")
if (ok) {
  g <- ggplotify::as.ggplot(function() grid.draw(pC))
  ggsave(file.path(ROOT, "04_results/recipeD.pdf"), g, width = 6.4, height = 3.4, device = cairo_pdf)
}
cat("DONE\n")
