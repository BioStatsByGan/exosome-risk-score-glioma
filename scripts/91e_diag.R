## 诊断：统计每个基因×组织在 MR 流程各步骤的通过数
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(TwoSampleMR)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
inst <- fread("04_results/MR/instruments_fixed.txt")
inst[, chr_i := as.integer(chr)]; inst[, pos_i := as.integer(pos)]
o <- fread("04_results/MR/outcome_GBM_ASTROCYTOMA.txt"); setnames(o, "#chrom", "chrom", skip_absent = TRUE)
o[, chrom := as.integer(gsub("chr", "", chrom))]; o[, pos := as.integer(pos)]
cat("inst:", nrow(inst), "| outcome:", nrow(o), "\n")

stat <- data.table()
for (g in unique(inst$external_gene_name)) {
  for (t in unique(inst$tissue)) {
    d <- inst[external_gene_name == g & tissue == t]
    if (!nrow(d)) next
    lead <- d[which.min(p)]
    om <- o[chrom == lead$chr_i & pos == lead$pos_i]
    s <- data.table(gene = g, tissue = t, matched = nrow(om) > 0, mrkeep = NA, mrok = NA)
    if (nrow(om)) {
      om1 <- om[1]
      exp <- data.frame(SNP = "s1", effect_allele.exposure = lead$A1, other_allele.exposure = lead$A2,
                        eaf.exposure = lead$eaf, beta.exposure = lead$beta, se.exposure = lead$se,
                        pval.exposure = lead$p, exposure = g, stringsAsFactors = FALSE)
      oc <- data.frame(SNP = "s1", effect_allele.outcome = toupper(om1$alt), other_allele.outcome = toupper(om1$ref),
                       eaf.outcome = om1$af_alt, beta.outcome = om1$beta, se.outcome = om1$sebeta,
                       pval.outcome = om1$pval, outcome = "GBM", stringsAsFactors = FALSE)
      h <- tryCatch(harmonise_data(exp, oc, action = 2), error = function(e) {cat('harmonise 报错:', conditionMessage(e), '\n'); NULL})
      s$mrkeep <- if (is.null(h)) NA else h$mr_keep[1]
      if (!is.null(h) && nrow(h) && h$mr_keep[1]) {
        m <- tryCatch(mr(h, method_list = "mr_wald_ratio"), error = function(e) NULL)
        s$mrok <- !is.null(m) && nrow(m) > 0
      } else s$mrok <- FALSE
    }
    stat <- rbind(stat, s)
  }
}
cat("\n=== 各步骤通过数 ===\n")
cat("基因组织组合总数:", nrow(stat), "\n")
cat("结局匹配成功:", sum(stat$matched), "\n")
cat("harmonise 后 mr_keep=TRUE:", sum(stat$mrkeep, na.rm = TRUE), "\n")
cat("mr() 成功:", sum(stat$mrok, na.rm = TRUE), "\n\n")
print(stat[matched == TRUE][1:min(15, .N)])
cat("\n=== mr_keep=FALSE 的位点检查（前 5 个）===\n")
bad <- stat[matched == TRUE & mrkeep == FALSE][1:min(5, .N)]
for (i in seq_len(nrow(bad))) {
  g <- bad$gene[i]; t <- bad$tissue[i]
  d <- inst[external_gene_name == g & tissue == t][which.min(p)]
  om <- o[chrom == d$chr_i & pos == d$pos_i][1]
  cat(sprintf("%s/%s: eQTL A1=%s A2=%s | FinnGen alt=%s ref=%s | eaf %s vs %s\n",
              g, t, d$A1, d$A2, om$alt, om$ref, round(d$eaf, 3), round(om$af_alt, 3)))
}
