## 调试 2：单基因单组织，跳过 clump 直接用 lead SNP 走 harmonise + MR
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(TwoSampleMR)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
inst <- fread("04_results/MR/instruments_all.txt")
out <- fread("04_results/MR/outcome_GBM_ASTROCYTOMA.txt")
setnames(out, "#chrom", "chrom", skip_absent = TRUE)

d <- inst[external_gene_name == "RAB27A" & tissue == "Brain_Cortex"]
d <- d[!is.na(rsids2) & p < 1e-5]
d <- d[which.min(d$p), ]
cat("lead SNP:", d$rsids2, " p =", d$p, " A1 =", d$A1, "\n")
cat("该 SNP 是否在结局文件:", d$rsids2 %in% out$rsids, "\n")
cat("结局文件里该基因附近任意 SNP:", sum(out$rsids %in% inst[external_gene_name == "RAB27A"]$rsids2), "\n")
cat("\n结局文件中前 5 个 rsID:", paste(head(out$rsids, 5), collapse = ", "), "\n")
cat("工具变量前 5 个 rsID:", paste(head(d$rsids2, 5), collapse = ", "), "\n")
cat("\n全部基因的工具变量与结局交集行数:", sum(inst$rsids2 %in% out$rsids), "\n")
tab <- inst[rsids2 %in% out$rsids, .N, by = external_gene_name]
print(tab)

o <- out[out$rsids %in% d$rsids2, ]
if (!nrow(o)) { cat("\n==> 结局文件缺该 SNP，MR 无法进行（需扩大结局 SNP 抽取范围）\n") } else {
  exp <- data.frame(SNP = d$rsids2, effect_allele.exposure = d$A1, other_allele.exposure = d$A2,
                    eaf.exposure = d$eaf, beta.exposure = d$beta, se.exposure = d$se,
                    pval.exposure = d$p, exposure = "RAB27A", id.exposure = "RAB27A", stringsAsFactors = FALSE)
  O <- data.frame(SNP = o$rsids, effect_allele.outcome = toupper(o$alt), other_allele.outcome = toupper(o$ref),
                  eaf.outcome = o$af_alt, beta.outcome = o$beta, se.outcome = o$sebeta,
                  pval.outcome = o$pval, outcome = "GBM", stringsAsFactors = FALSE)
  h <- harmonise_data(exp, O, action = 2)
  cat("\nharmonise 结果:\n"); print(h[, c("SNP", "effect_allele.exposure", "effect_allele.outcome", "beta.exposure", "beta.outcome", "mr_keep")])
  if (sum(h$mr_keep) > 0) {
    r <- mr(h[h$mr_keep, ], method_list = "mr_wald_ratio")
    print(r); print(generate_odds_ratios(r))
  }
}
