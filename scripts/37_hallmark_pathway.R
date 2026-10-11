## 步骤 37：高/低风险组 Hallmark 通路活性（GSVA）+ 分级校正
## 输出：03_figures/Fig6C_hallmark_heatmap.pdf, 03_figures/Fig6D_hallmark_violin.pdf
##       04_results/37.hallmark/hallmark_diff.txt, hallmark_scores.txt
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({
  library(data.table); library(ggplot2); library(limma); library(msigdbr)
  library(GSVA); library(patchwork); library(tidyr)
})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
setwd(ROOT); source("scripts/plot_style.R")
OUTD <- "04_results/37.hallmark"; dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)

## ---- 表达矩阵（基因 × 样本）----
ex <- fread("02_inputs/TCGA_glioma_tumor_log2tpm.csv")
gn <- ex[[1]]; m <- as.matrix(ex[, -1]); rownames(m) <- gn
m <- m[!duplicated(rownames(m)), ]
cat("表达矩阵:", nrow(m), "基因 ×", ncol(m), "样本\n")

## ---- 风险分组 ----
rk <- fread("04_results/20.SHAP/risk.TCGA.txt")
rk <- rk[!is.na(futime) & !is.na(fustat)]
cl <- fread("02_inputs/clinical_glioma.txt")
rk <- merge(rk[, .(id, riskScore, Risk)], cl[, .(patient, grade, idh)],
            by.x = "id", by.y = "patient", all.x = TRUE)
rk[, Risk := factor(ifelse(Risk == "high", "High", "Low"), levels = c("Low", "High"))]
sam <- intersect(rk$id, colnames(m))
cat("交集样本:", length(sam), "\n")
m <- m[, sam]; rk <- rk[match(sam, rk$id)]

## ---- Hallmark 基因集 ----
hm <- as.data.frame(msigdbr(species = "Homo sapiens", category = "H"))
sets <- split(hm$gene_symbol, hm$gs_name)
sets <- lapply(sets, function(s) intersect(unique(s), rownames(m)))
sets <- sets[lengths(sets) >= 10]
cat("Hallmark 通路数:", length(sets), "\n")

## ---- GSVA ----
## 仅保留 Hallmark 基因，降低 GSVA 内存/时间开销（约 4k 基因）
hg <- unique(unlist(sets))
m2 <- m[intersect(rownames(m), hg), ]
cat("Hallmark 基因矩阵:", nrow(m2), "基因 ×", ncol(m2), "样本\n")

cache <- file.path(OUTD, "hallmark_scores.txt")
res <- if (file.exists(cache)) {
  cat("复用已算好的 GSVA 分数缓存\n")
  c2 <- fread(cache); mm <- as.matrix(c2[, -1]); rownames(mm) <- c2[[1]]; mm
} else {
  tryCatch({
    gp <- gsvaParam(m2, sets, kcdf = "Gaussian", minSize = 10, maxSize = 500)
    gsva(gp, verbose = FALSE)
  }, error = function(e) {
    gsva(m2, sets, method = "gsva", kcdf = "Gaussian", min.sz = 10, max.sz = 500, verbose = FALSE)
  })
}
cat("GSVA 输出:", nrow(res), "通路 ×", ncol(res), "样本\n")
write.table(data.frame(Pathway = rownames(res), res, check.names = FALSE),
            file.path(OUTD, "hallmark_scores.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

## ---- 组间差异 + 分级校正 ----
gd <- data.frame(Pathway = rownames(res), res, check.names = FALSE)
long <- melt(as.data.table(gd), id.vars = "Pathway", variable.name = "id", value.name = "score")
long <- merge(long, rk[, .(id, Risk, riskScore, grade, idh)], by = "id")
stat <- long[, {
  w <- wilcox.test(score ~ Risk)
  d <- median(score[Risk == "High"]) - median(score[Risk == "Low"])
  pv <- tryCatch({ summary(lm(score ~ Risk + factor(grade) + factor(idh)))$coefficients["RiskHigh", 4] },
                 error = function(e) NA_real_)
  .(delta = d, P = w$p.value, P_grade_adj = pv,
    rho = cor(score, riskScore, method = "spearman"))
}, by = Pathway]
stat[, FDR := p.adjust(P, "BH")]
setorder(stat, -delta)
fwrite(stat, file.path(OUTD, "hallmark_diff.txt"), sep = "\t")
cat("\n== 高风险的 Top10 通路 ==\n"); print(head(stat, 10))
cat("\n== 低风险的 Top5 通路 ==\n"); print(tail(stat, 5))
cat("\n分级校正后仍 P<0.05 的通路:", sum(stat$P_grade_adj < 0.05, na.rm = TRUE), "/", nrow(stat), "\n")

## ================= Fig6C 热图：Top24 通路 × 样本 =================
EXCL <- c("HALLMARK_PANCREAS_BETA_CELLS", "HALLMARK_HEME_METABOLISM")
s_plot <- stat[!Pathway %in% EXCL]
top <- head(s_plot[order(-abs(delta))], 24)$Pathway
sub <- as.data.table(t(scale(t(res[top, ]))))
sub[, Pathway := rownames(res[top, ])]
lg <- melt(sub, id.vars = "Pathway", variable.name = "id", value.name = "z")
lg <- merge(lg, rk[, .(id, Risk, riskScore)], by = "id")
lg[, id := factor(id, levels = rk[order(riskScore), id])]
lg[, Pathway := factor(Pathway, levels = rev(stat$Pathway[stat$Pathway %in% top]))]
lab <- gsub("^HALLMARK_", "", levels(lg$Pathway))
lab <- gsub("_", " ", lab); lab <- paste0(toupper(substr(lab, 1, 1)), tolower(substr(lab, 2, nchar(lab))))
pC <- ggplot(lg, aes(id, Pathway, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0,
                       name = "GSVA z-score") +
  scale_y_discrete(labels = rev(lab)) +
  labs(x = "Patients ordered by risk score", y = NULL) +
  theme_pub() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank(), axis.text.y = element_text(size = 6.2),
        legend.key.height = unit(0.7, "cm"), legend.key.width = unit(0.28, "cm"))
# 顶部风险条
bar <- ggplot(unique(lg[, .(id, riskScore, Risk)]), aes(id, 1, fill = Risk)) +
  geom_tile() + scale_fill_manual(values = c(Low = OI[1], High = OI[2]), name = "") +
  theme_void() + theme(legend.position = "bottom", legend.key.size = unit(0.3, "cm"),
                       legend.text = element_text(size = 6.5))
save_pdf(bar / pC + plot_layout(heights = c(1, 22)), "03_figures/Fig6C_hallmark_heatmap.pdf", 7.4, 5.6)

## ================= Fig6D 小提琴：Top8 通路 =================
top8 <- head(s_plot[order(-abs(delta))], 8)$Pathway
dd <- long[Pathway %in% top8]
dd[, Plab := gsub("_", " ", sub("^HALLMARK_", "", Pathway))]
dd[, Plab := factor(Plab, levels = gsub("_", " ", sub("^HALLMARK_", "", s_plot[order(-abs(delta))]$Pathway[s_plot$Pathway %in% top8])))]
ann <- stat[Pathway %in% top8, .(Pathway, P, P_grade_adj)]
ann[, Plab := factor(gsub("_", " ", sub("^HALLMARK_", "", Pathway)),
                     levels = levels(dd$Plab))]
pD <- ggplot(dd, aes(Plab, score, fill = Risk)) +
  geom_violin(scale = "width", width = 0.85, colour = NA, alpha = 0.85) +
  geom_boxplot(width = 0.14, outlier.size = 0.25, linewidth = 0.28, fill = "white", colour = "grey25") +
  scale_fill_manual(values = c(Low = OI[1], High = OI[2]), name = "") +
  labs(x = NULL, y = "GSVA score") +
  coord_flip() +
  geom_text(data = ann, aes(x = Plab, y = Inf, label = sprintf("p=%.1g | adj p=%.2g", P, P_grade_adj)),
            inherit.aes = FALSE, hjust = 1.02, vjust = -0.5, size = 2.3, colour = "grey25") +
  theme_pub() + theme(legend.position = c(0.85, 0.12))
save_pdf(pD, "03_figures/Fig6D_hallmark_violin.pdf", 6.6, 4.4)
cat("\nFig6C / Fig6D 完成\n")
