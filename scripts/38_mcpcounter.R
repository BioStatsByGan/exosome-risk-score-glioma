## 步骤 38：MCP-counter 第三方去卷积 —— 解决中性粒细胞的 CIBERSORT/ssGSEA 结论冲突
## 输出：03_figures/Fig6E_mcpcounter.pdf
##       04_results/38.MCPcounter/mcp_diff.txt, mcp_scores.txt, neutrophil_method_compare.txt
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({
  library(data.table); library(ggplot2); library(MCPcounter); library(patchwork)
})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
setwd(ROOT); source("scripts/plot_style.R")
OUTD <- "04_results/38.MCPcounter"; dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)

## ---- 表达矩阵 ----
ex <- fread("02_inputs/TCGA_glioma_tumor_log2tpm.csv")
m <- as.matrix(ex[, -1]); rownames(m) <- ex[[1]]
m <- m[!duplicated(rownames(m)), ]

## ---- 风险与临床 ----
rk <- fread("04_results/20.SHAP/risk.TCGA.txt"); rk <- rk[!is.na(futime) & !is.na(fustat)]
cl <- fread("02_inputs/clinical_glioma.txt")
rk <- merge(rk[, .(id, riskScore, Risk)], cl[, .(patient, grade, idh)],
            by.x = "id", by.y = "patient", all.x = TRUE)
rk[, Risk := factor(ifelse(Risk == "high", "High", "Low"), levels = c("Low", "High"))]
sam <- intersect(rk$id, colnames(m)); m <- m[, sam]; rk <- rk[match(sam, rk$id)]
cat("样本:", length(sam), "\n")

## ---- MCP-counter ----
mcp <- MCPcounter.estimate(m, featuresType = "HUGO_symbols")
cat("MCP-counter 细胞群体:", paste(rownames(mcp), collapse = " | "), "\n")
write.table(data.frame(Population = rownames(mcp), mcp, check.names = FALSE),
            file.path(OUTD, "mcp_scores.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

lg <- as.data.table(as.table(t(mcp)))
setnames(lg, c("id", "Population", "score"))
lg[, score := as.numeric(score)]
lg <- merge(lg, rk[, .(id, Risk, riskScore, grade, idh)], by = "id")
stat <- lg[, {
  w <- wilcox.test(score ~ Risk)
  d <- median(score[Risk == "High"]) - median(score[Risk == "Low"])
  pv <- tryCatch(summary(lm(score ~ Risk + factor(grade) + factor(idh)))$coefficients["RiskHigh", 4],
                 error = function(e) NA_real_)
  .(delta = d, P = w$p.value, P_grade_adj = pv, rho = cor(score, riskScore, method = "spearman"))
}, by = Population]
stat[, FDR := p.adjust(P, "BH")]
setorder(stat, -delta)
fwrite(stat, file.path(OUTD, "mcp_diff.txt"), sep = "\t")
cat("\n== MCP-counter 组间差异 ==\n"); print(stat)

## ---- 与 CIBERSORT / ssGSEA 的中性粒一致性 ----
cb <- fread("04_results/27.CIBERSORT/CIBERSORT-Results.txt")
setnames(cb, 1, "id")
neu_cb <- cb[, .(id, CIBERSORT_neutrophil = `Neutrophils`)]
neu_mcp <- lg[Population == "Neutrophils", .(id, MCPcounter_neutrophil = score)]
ss <- fread("04_results/29.immFunction/ssgsea_scores.txt")
neu_col <- grep("eutrophil", names(ss), value = TRUE)[1]
neu_ss <- if (!is.na(neu_col)) {
  tmp <- as.data.table(ss); setnames(tmp, 1, "id"); tmp[, .(id, ssGSEA_neutrophil = get(neu_col))]
} else NULL
cmp <- merge(neu_mcp, neu_cb, by = "id")
if (!is.null(neu_ss)) cmp <- merge(cmp, neu_ss, by = "id")
cmp[, MCPcounter_neutrophil := as.numeric(MCPcounter_neutrophil)]
r1 <- cor(cmp$MCPcounter_neutrophil, cmp$CIBERSORT_neutrophil, method = "spearman")
cat(sprintf("\nMCP-counter 中性粒 vs CIBERSORT 中性粒: rho = %.3f (n=%d)\n", r1, nrow(cmp)))
if (!is.null(neu_ss)) {
  r2 <- cor(cmp$MCPcounter_neutrophil, cmp$ssGSEA_neutrophil, method = "spearman")
  r3 <- cor(cmp$CIBERSORT_neutrophil, cmp$ssGSEA_neutrophil, method = "spearman")
  cat(sprintf("MCP-counter vs ssGSEA: rho = %.3f | CIBERSORT vs ssGSEA: rho = %.3f\n", r2, r3))
}
fwrite(cmp, file.path(OUTD, "neutrophil_method_compare.txt"), sep = "\t")

## ================= Fig6E =================
d <- copy(stat)
d[, Plab := gsub("_", " ", Population)]
d[, Sig := ifelse(P_grade_adj < 0.05, "Grade-adjusted P<0.05", "Not significant")]
d[, Plab := factor(Plab, levels = rev(Plab[order(delta)]))]
pE1 <- ggplot(d, aes(x = delta, y = Plab)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = delta, yend = Plab), colour = "grey70", linewidth = 0.4) +
  geom_point(aes(colour = Sig), size = 2.1) +
  scale_colour_manual(values = c("Grade-adjusted P<0.05" = OI[2], "Not significant" = OI[7]), name = "") +
  labs(x = "Difference in abundance score (high - low risk)", y = NULL) +
  theme_pub() + theme(legend.position = "bottom", legend.key.size = unit(0.3, "cm"),
                      legend.text = element_text(size = 6.6))

pE2 <- ggplot(cmp, aes(MCPcounter_neutrophil, CIBERSORT_neutrophil)) +
  geom_point(size = 0.35, colour = OI[1], alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = OI[2], linewidth = 0.5) +
  annotate("text", x = min(cmp$MCPcounter_neutrophil), y = max(cmp$CIBERSORT_neutrophil),
           label = sprintf("Spearman rho = %.2f", r1), hjust = 0, vjust = 1, size = 2.6) +
  labs(x = "MCP-counter neutrophil (log scale)", y = "CIBERSORT neutrophil fraction") +
  theme_pub()
save_pdf(pE1 + pE2 + plot_layout(widths = c(1.15, 1)), "03_figures/Fig6E_mcpcounter.pdf", 9.4, 4.2)
cat("\nFig6E 完成\n")
