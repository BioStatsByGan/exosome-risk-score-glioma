## 检查 TCGA 免疫分型注释文件内容（aedin/TCGAimmune）
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(SummarizedExperiment)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
e <- new.env()
load(file.path(ROOT, "00_data/TCGA_immune_MAE.rda"), envir = e)
cat("对象:", ls(e), "\n")
obj <- get(ls(e)[1], envir = e)
cat("类:", paste(class(obj), collapse = ", "), "\n")
cd <- tryCatch(as.data.frame(SummarizedExperiment::colData(obj)), error = function(x) NULL)
if (!is.null(cd)) {
  cat("colData 列名:\n"); print(colnames(cd))
  sub <- cd[, grep("Subtype|Immune", colnames(cd), value = TRUE), drop = FALSE]
  print(head(sub, 5))
  cat("样本数:", nrow(cd), "\n")
  for (cn in colnames(sub)) { cat("\n", cn, "取值分布:\n"); print(table(sub[[cn]], useNA = "ifany")) }
} else {
  cat("colData 读取失败，改为查看顶层结构\n"); print(str(obj, max.level = 2))
}
