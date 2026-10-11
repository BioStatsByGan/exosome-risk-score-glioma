## MR 层（最终版）：脑组织 lead cis-eQTL → FinnGen R12 胶质瘤，两样本 MR + 共定位
## 设计：每个基因×脑组织取最显著的 lead cis-eQTL 作为工具变量（单 SNP，Wald ratio），
##       规避 LD 冗余；工具强度用 F=(beta/se)^2 表示；共定位用该基因全部可得 SNP。
## 输出：04_results/MR/{mr_results.txt, coloc_results.txt}、03_figures/Fig10A_mr_forest.pdf、Fig10B_mr_coloc.pdf
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({
  library(data.table); library(TwoSampleMR); library(coloc)
  library(ggplot2); library(forestploter); library(grid)
})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT); source("scripts/plot_style.R")

N_TISSUE <- c(Brain_Amygdala = 129, Brain_Anterior_cingulate_cortex_BA24 = 147,
              Brain_Caudate_basal_ganglia = 194, Brain_Cerebellar_Hemisphere = 175,
              Brain_Cerebellum = 209, Brain_Cortex = 205, Brain_Frontal_Cortex_BA9 = 175,
              Brain_Hippocampus = 165, Brain_Hypothalamus = 170,
              Brain_Nucleus_accumbens_basal_ganglia = 202, Brain_Putamen_basal_ganglia = 170,
              `Brain_Spinal_cord_cervical_c-1` = 126, Brain_Substantia_nigra = 114)
OUTCOME_META <- list(GBM_ASTROCYTOMA = list(ncase = 544, nctrl = 378749, label = "Glioma (GBM and astrocytoma)"),
                     BRAIN_EXALLC = list(ncase = 1894, nctrl = 378749, label = "Malignant neoplasm of brain"))

inst <- fread("04_results/MR/instruments_fixed.txt")   # 已用 Ensembl 校正染色体编号
inst <- inst[!is.na(chr) & !is.na(pos)]
inst[, chr_i := as.integer(chr)]        # 统一成整数，避免 chr 字符型与结局整数型比较失败
inst[, pos_i := as.integer(pos)]
cat("inst$chr class:", class(inst$chr), "-> chr_i class:", class(inst$chr_i), "\n")
read_out <- function(tag) {
  f <- paste0("04_results/MR/outcome_", tag, ".txt")
  if (!file.exists(f)) return(NULL)
  o <- fread(f); setnames(o, "#chrom", "chrom", skip_absent = TRUE)
  o[, chrom := as.integer(gsub("chr", "", chrom))]
  o[, pos := as.integer(pos)]
  o
}
out1 <- read_out("GBM_ASTROCYTOMA"); out2 <- read_out("BRAIN_EXALLC")
cat("工具变量行:", nrow(inst), "| 结局位点(主):", nrow(out1), "| 敏感性:", if (is.null(out2)) 0 else nrow(out2), "\n")

res <- list(); coloc <- list(); nmatch <- 0
for (g in unique(inst$external_gene_name)) {
  for (t in unique(inst$tissue)) {
    d <- inst[external_gene_name == g & tissue == t]
    if (!nrow(d)) next
    lead <- d[which.min(p)]
    n1 <- if (t %in% names(N_TISSUE)) N_TISSUE[[t]] else NA
    for (tag in c("GBM_ASTROCYTOMA", "BRAIN_EXALLC")) {
      O <- if (tag == "GBM_ASTROCYTOMA") out1 else out2
      if (is.null(O)) next
      om <- O[chrom == lead$chr_i & pos == lead$pos_i]
      if (!nrow(om)) next
      nmatch <- nmatch + 1
      om <- om[1]
      ma <- OUTCOME_META[[tag]]
      exp <- data.frame(SNP = paste0(g, "_", t), effect_allele.exposure = lead$A1, other_allele.exposure = lead$A2,
                        eaf.exposure = lead$eaf, beta.exposure = lead$beta, se.exposure = lead$se,
                        pval.exposure = lead$p, exposure = paste0(g, " (", t, ")"),
                        id.exposure = paste0(g, "_", t), stringsAsFactors = FALSE)
      oc <- data.frame(SNP = exp$SNP, effect_allele.outcome = toupper(om$alt), other_allele.outcome = toupper(om$ref),
                       eaf.outcome = om$af_alt, beta.outcome = om$beta, se.outcome = om$sebeta,
                       pval.outcome = om$pval, outcome = ma$label, id.outcome = tag,
                       stringsAsFactors = FALSE)
      h <- tryCatch(harmonise_data(exp, oc, action = 2), error = function(e) NULL)
      if (is.null(h) || !nrow(h) || !h$mr_keep[1]) next
      m <- tryCatch(mr(h, method_list = "mr_wald_ratio"), error = function(e) NULL)
      if (is.null(m) || !nrow(m)) next
      ci <- generate_odds_ratios(m)
      st <- tryCatch(directionality_test(h), error = function(e) NULL)
      ## 单 SNP Wald ratio 下 directionality_test 返回 NaN，改用手算 r2 比较（标准 Steiger 做法）
      maf <- min(lead$eaf, 1 - lead$eaf)
      r2_exp <- 2 * maf * (1 - maf) * lead$beta^2
      r2_out <- 2 * maf * (1 - maf) * om$beta^2
      steiger_manual <- r2_exp > r2_out
      res[[paste(g, t, tag)]] <- data.frame(
        Gene = g, Tissue = t, Outcome = tag, nSNP = 1,
        SNP = paste0("chr", lead$chr, ":", lead$pos, ifelse(is.na(lead$rsids2), "", paste0(" (", lead$rsids2, ")"))),
        eQTL_p = lead$p, OR = ci$or, OR_low = ci$or_lci95, OR_high = ci$or_uci95, p = m$pval,
        F_stat = (lead$beta / lead$se)^2,
        Steiger = steiger_manual, r2_exp = r2_exp, r2_out = r2_out,
        stringsAsFactors = FALSE)
      ## 共定位：该基因在该组织全部 SNP
      if (tag == "GBM_ASTROCYTOMA") {
        d[, `:=`(chr_i = as.integer(chr), pos_i = as.integer(pos))]
        mm <- merge(d[, .(chr_i, pos_i, A1, A2, eaf, beta, se)],
                    O[, .(chr_i = as.integer(chrom), pos_i = as.integer(pos), alt, ref, af_alt, beta2 = beta, se2 = sebeta)],
                    by = c("chr_i", "pos_i"))
        if (nrow(mm) >= 5) {
          cr <- tryCatch({
            r <- coloc.abf(dataset1 = list(beta = mm$beta, varbeta = mm$se^2, snp = paste0(mm$chr_i, ":", mm$pos_i),
                                           N = n1, type = "quant", MAF = pmin(mm$eaf, 1 - mm$eaf)),
                           dataset2 = list(beta = mm$beta2, varbeta = mm$se2^2, snp = paste0(mm$chr_i, ":", mm$pos_i),
                                           N = ma$ncase + ma$nctrl, type = "cc", s = ma$ncase / (ma$ncase + ma$nctrl)))
            s <- r$summary
            data.frame(Gene = g, Tissue = t, nSNP = nrow(mm), PP_H3 = s[["PP.H3.abf"]], PP_H4 = s[["PP.H4.abf"]],
                       stringsAsFactors = FALSE)
          }, error = function(e) NULL)
          if (!is.null(cr)) coloc[[paste(g, t)]] <- cr
        }
      }
    }
  }
}

mr_tab <- rbindlist(res)
mr_tab$FDR <- p.adjust(mr_tab$p, "BH")
setorder(mr_tab, p)
fwrite(mr_tab, "04_results/MR/mr_results.txt", sep = "\t", quote = FALSE)
cat("\n===== MR 结果（主结局）=====\n")
print(mr_tab[Outcome == "GBM_ASTROCYTOMA"][order(p)][1:min(20, .N)], row.names = FALSE, digits = 3)
cat("\n主结局检验数:", nrow(mr_tab[Outcome == "GBM_ASTROCYTOMA"]),
    "| 名义 P<0.05:", nrow(mr_tab[Outcome == "GBM_ASTROCYTOMA" & p < 0.05]),
    "| FDR<0.05:", nrow(mr_tab[Outcome == "GBM_ASTROCYTOMA" & FDR < 0.05]), "\n")
cat("敏感性结局 BRAIN_EXALLC 检验数:", nrow(mr_tab[Outcome == "BRAIN_EXALLC"]),
    "| 名义 P<0.05:", nrow(mr_tab[Outcome == "BRAIN_EXALLC" & p < 0.05]), "\n")

if (length(coloc)) {
  ct <- rbindlist(coloc); setorder(ct, -PP_H4)
  fwrite(ct, "04_results/MR/coloc_results.txt", sep = "\t", quote = FALSE)
  cat("\n===== 共定位 PP.H4 > 0.5 =====\n"); print(ct[PP_H4 > 0.5], row.names = FALSE, digits = 3)
  cat("PP.H4 最大值:", round(max(ct$PP_H4), 3), "（评价基因数:", nrow(ct), "）\n")
}

## 森林图
sig <- mr_tab[Outcome == "GBM_ASTROCYTOMA"][order(p)][1:min(20, .N)][order(OR)]
if (nrow(sig)) {
  lab <- paste0(sig$Gene, " | ", gsub("Brain_", "", sig$Tissue))
  tab <- data.frame(Gene = lab,
                    `OR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", sig$OR, sig$OR_low, sig$OR_high),
                    `P value` = ifelse(sig$p < 0.001, "<0.001", sprintf("%.3f", sig$p)),
                    ` ` = paste(rep(" ", 22), collapse = ""), check.names = FALSE)
  tb <- forest_theme(base_size = 11, refline_gp = gpar(col = "grey45", lty = 2),
                     ci_col = "#0072B2", ci_lwd = 1.5, ci_Theight = 0.22,
                     core = list(bg_params = list(fill = c("white", "grey97"))))
  xmax <- max(sig$OR_high) * 1.1
  p_f <- forest(tab, est = sig$OR, lower = sig$OR_low, upper = sig$OR_high, sizes = 0.75,
                ci_column = 4, ref_line = 1, xlim = c(0, xmax), ticks_at = c(0, round(xmax / 2, 2), round(xmax, 2)), theme = tb)
  hgt <- max(3.0, 0.32 * nrow(sig) + 2.4)
  grDevices::cairo_pdf("03_figures/Fig10A_mr_forest.pdf", width = 7.4, height = hgt, family = "sans"); print(p_f); dev.off()
  grDevices::png("03_figures/Fig10A_mr_forest.png", width = 7.4, height = hgt, units = "in", res = 200, type = "cairo"); print(p_f); dev.off()
  cat("saved: Fig10A_mr_forest.pdf\n")
}
if (length(coloc)) {
  ct2 <- rbindlist(coloc)[order(-PP_H4)][1:min(12, .N)]
  ct2$Label <- factor(paste0(ct2$Gene, " | ", gsub("Brain_", "", ct2$Tissue)), levels = rev(paste0(ct2$Gene, " | ", gsub("Brain_", "", ct2$Tissue))))
  pc <- ggplot(ct2, aes(Label, PP_H4, fill = PP_H4)) +
    geom_col(width = 0.62, colour = "grey25", linewidth = 0.3) +
    geom_hline(yintercept = 0.8, linetype = "dashed", colour = "#D55E00", linewidth = 0.35) +
    scale_fill_gradient(low = "#CCE5F5", high = "#D55E00", guide = "none") +
    labs(x = NULL, y = "PP.H4 (shared causal variant)") + coord_flip()
  save_pdf(pc, "03_figures/Fig10B_mr_coloc.pdf", 6.4, 4.6)
}
cat("MR DONE\n")
