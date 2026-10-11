## 调试：单基因单组织走一遍 MR 流程，定位哪一步出问题
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(TwoSampleMR); library(ieugwasr); library(dplyr)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
PLINK <- "E:/01dta/本地clump/bin/plink.exe"; BFILE <- "E:/01dta/本地clump/bin/EUR/EUR"

inst <- fread("04_results/MR/instruments_all.txt")
out <- fread("04_results/MR/outcome_GBM_ASTROCYTOMA.txt")
cat("inst 行:", nrow(inst), "| 列:", paste(names(inst), collapse=","), "\n")
cat("out 行:", nrow(out), "| 列:", paste(names(out), collapse=","), "\n")
cat("p<5e-8 的行数:", sum(inst$p < 5e-8), "\n")

d <- inst[external_gene_name == "RAB27A" & tissue == "Brain_Cortex"]
cat("\nRAB27A/Brain_Cortex 原始行:", nrow(d), " p 范围:", range(d$p), "\n")
d <- d[!is.na(rsids2) & p < 1e-5]
cat("过滤后:", nrow(d), "\n")
cl <- ld_clump(dplyr::tibble(rsid = d$rsids2, pval = d$p), plink_bin = PLINK, bfile = BFILE,
               clump_kb = 10000, clump_r2 = 0.001, clump_p = 1)
cat("clump 后 SNP:", nrow(cl), paste(cl$rsid, collapse = ","), "\n")
d2 <- d[rsids2 %in% cl$rsid]
cat("保留行:", nrow(d2), "\n")

exp <- data.frame(SNP = d2$rsids2, effect_allele.exposure = d2$A1, other_allele.exposure = d2$A2,
                  eaf.exposure = d2$eaf, beta.exposure = d2$beta, se.exposure = d2$se,
                  pval.exposure = d2$p, exposure = "test", id.exposure = "test",
                  samplesize.exposure = 205, stringsAsFactors = FALSE)
o <- out[out$rsids %in% exp$SNP, ]
cat("结局命中 SNP:", nrow(o), "\n")
O <- data.frame(SNP = o$rsids, effect_allele.outcome = toupper(o$alt), other_allele.outcome = toupper(o$ref),
                eaf.outcome = o$af_alt, beta.outcome = o$beta, se.outcome = o$sebeta,
                pval.outcome = o$pval, outcome = "GBM", samplesize.outcome = 544 + 378749,
                stringsAsFactors = FALSE)
h <- harmonise_data(exp, O, action = 2)
cat("harmonise 行:", nrow(h), "| mr_keep TRUE:", sum(h$mr_keep), "\n")
print(h[, c("SNP", "effect_allele.exposure", "effect_allele.outcome", "beta.exposure", "beta.outcome", "mr_keep")])
if (sum(h$mr_keep) > 0) {
  h2 <- h[h$mr_keep, ]
  mrres <- mr(h2, method_list = c("mr_wald_ratio", "mr_ivw"))
  print(mrres)
  ci <- generate_odds_ratios(mrres)
  print(ci)
  st <- directionality_test(h2); print(st)
}
