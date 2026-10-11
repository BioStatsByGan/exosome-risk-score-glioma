## 步骤 13-15（胶质瘤版）：外泌体相关基因差异分析 + 富集 + 发表级出图
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({
  library(data.table); library(limma); library(sva); library(ggplot2); library(pheatmap)
  library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot); library(patchwork)
  library(ggrepel); library(ComplexHeatmap); library(circlize)
})
source("F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤/scripts/plot_style.R")
base <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; inp <- file.path(base, "02_inputs")
fig  <- file.path(base, "03_figures"); dir.create(fig, recursive = TRUE, showWarnings = FALSE)
out  <- file.path(base, "04_results/13.diff"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
rd <- function(f) fread(f, data.table = FALSE, check.names = FALSE)

## ---------- 数据 ----------
tum <- rd(file.path(inp, "TCGA_glioma_tumor_log2tpm.csv"))
gex <- rd(file.path(inp, "GTEx_brain_log2tpm.csv"))
eg  <- readLines(file.path(inp, "exosome_genes.txt"))
eg  <- unique(eg[nzchar(eg)])
cat("肿瘤矩阵:", dim(tum), " GTEx 脑:", dim(gex), " 外泌体基因:", length(eg), "\n")

rownames(tum) <- tum[[1]]; tum <- as.matrix(tum[, -1])
rownames(gex) <- gex[[1]]; gex <- as.matrix(gex[, -1])
gene <- intersect(intersect(rownames(tum), rownames(gex)), eg)
cat("三者交集基因:", length(gene), "\n")
X <- cbind(tum[gene, , drop = FALSE], gex[gene, , drop = FALSE])
grp <- factor(c(rep("Tumor", ncol(tum)), rep("Normal", ncol(gex))), levels = c("Normal", "Tumor"))
cat("合并矩阵:", dim(X), " 分组:", paste(levels(grp), as.vector(table(grp)), sep = "=", collapse = ", "), "\n")

## ---------- 批次校正（TCGA vs GTEx，保留分组效应）与差异分析（带缓存） ----------
cache <- file.path(out, "deg_cache.rds")
if (file.exists(cache)) {
  z <- readRDS(cache); Xc <- z$Xc; deg <- z$deg
} else {
  X <- X[rowMeans(X) > 0.5, , drop = FALSE]
  design <- model.matrix(~ grp)
  Xc <- tryCatch(ComBat(X, batch = ifelse(grp == "Tumor", "TCGA", "GTEx"), mod = design, par.prior = TRUE),
                 error = function(e) { cat("ComBat 失败，改用未校正矩阵\n"); X })
  fit <- eBayes(lmFit(Xc, design))
  deg <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  deg$gene <- rownames(deg)
  saveRDS(list(Xc = Xc, deg = deg), cache)
}
cat("差异基因 FDR<0.05:", sum(deg$adj.P.Val < 0.05), " |logFC|>1 且 FDR<0.05:", sum(deg$adj.P.Val < 0.05 & abs(deg$logFC) > 1), "\n")
write.table(deg, file.path(out, "glioma_exosome_all.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
dif <- deg[deg$adj.P.Val < 0.05 & abs(deg$logFC) > 1, ]
write.table(dif, file.path(out, "glioma_exosome_diff.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

## ---------- 图 1 火山图（渐变配色） ----------
d <- data.frame(gene = deg$gene, logFC = deg$logFC, fdr = deg$adj.P.Val)
p1 <- volcano_pub(d, fc_cut = 1, fdr_cut = 0.05, label_n = 16)
save_pdf(p1, file.path(fig, "Fig1A_volcano.pdf"), 5.4, 4.6)

## ---------- 图 2 热图（前 50 差异基因，ComplexHeatmap） ----------
hm <- head(dif$gene[order(-abs(dif$logFC))], 50)
z <- t(scale(t(Xc[hm, , drop = FALSE])))
z[z > 2] <- 2; z[z < -2] <- -2
gfac <- factor(grp, levels = c("Normal", "Tumor"))
ha <- HeatmapAnnotation(Group = gfac,
                        col = list(Group = c(Normal = "#0072B2", Tumor = "#D55E00")),
                        annotation_name_gp = gpar(fontsize = 8))
ht <- Heatmap(z, name = "z-score", top_annotation = ha, show_column_names = FALSE,
              cluster_columns = TRUE, column_split = gfac,
              row_names_gp = gpar(fontsize = 6.5), column_title_gp = gpar(fontsize = 9),
              col = colorRamp2(c(-2, 0, 2), c("#2c7bb6", "white", "#d7191c")),
              border = TRUE, border_gp = gpar(col = "black", lwd = 0.4))
pdf(file.path(fig, "Fig1B_heatmap.pdf"), width = 6.4, height = 6.8)
draw(ht); dev.off()
png(sub("\\.pdf$", ".png", file.path(fig, "Fig1B_heatmap.pdf")), width = 6.4, height = 6.8, units = "in", res = 200)
draw(ht); dev.off()
cat("saved: Fig1B_heatmap.pdf\n")

## ---------- 富集 ----------
eg_sym <- bitr(dif$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go <- enrichGO(eg_sym$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
               pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
kk <- enrichKEGG(eg_sym$ENTREZID, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
write.table(as.data.frame(go), file.path(out, "GO.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(as.data.frame(kk), file.path(out, "KEGG.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

## ---------- 图 3 GSEA 曲线（替代气泡图） ----------
ranked <- deg$gene[order(deg$logFC, decreasing = TRUE)]
curves <- list()
if (nrow(as.data.frame(kk)) > 0) {
  kkdf <- as.data.frame(kk)
  for (i in 1:min(3, nrow(kkdf))) {
    gs <- eg_sym$SYMBOL[match(strsplit(kkdf$geneID[i], "/")[[1]], eg_sym$ENTREZID)]
    gs <- gs[!is.na(gs)]
    curves[[i]] <- gsea_curve(ranked, gs, col = OI[i]) + ggtitle(paste0("KEGG: ", kkdf$Description[i]))
  }
}
if (length(curves) > 0) {
  p3 <- wrap_plots(curves, ncol = 1)
  save_pdf(p3, file.path(fig, "Fig1C_gsea_curves.pdf"), 5.6, 7.2)
}

## ---------- 图 4 富集网络图（cnetplot，替代气泡） ----------
if (nrow(as.data.frame(go)) > 0) {
  pdf(file.path(fig, "Fig1D_cnetplot.pdf"), width = 7.5, height = 6.2)
  print(cnetplot(go, showCategory = 5) + theme_pub())
  dev.off()
  png(sub("\\.pdf$", ".png", file.path(fig, "Fig1D_cnetplot.pdf")), width = 7.5, height = 6.2, units = "in", res = 200)
  print(cnetplot(go, showCategory = 5) + theme_pub()); dev.off()
  cat("saved: Fig1D_cnetplot.pdf\n")
}
cat("STEP13 DONE\n")
