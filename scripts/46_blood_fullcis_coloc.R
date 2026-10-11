## 步骤 46：血液端「全 cis 窗」共定位（eQTLGen 100kb 全 SNP × FinnGen R12 全量 GWAS）
## 背景：本地 GTEx「lite」besd 实测只含显著对（Brain_Cortex 仅 6 866 个探针），无法做全窗共定位；
##       eQTLGen 100kb cis 为每 SNP 的完整 beta/se（n=30 895），可做规范的全窗 coloc。
## 输出：04_results/audit/blood_fullcis_coloc.txt / _log.txt
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(coloc)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
OUTD <- "04_results/audit"; dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
LOG <- character(); say <- function(...) { s <- sprintf(...); cat(s, "\n"); LOG <<- c(LOG, s); invisible(s) }

## ---------- 1) 暴露：eQTLGen 100kb 全 cis ----------
eq <- readRDS("I:/05-eqtlgen_IEU_Source_Data/06_cis_colocation_data/eqtlgen_ieu_100kbcis.rds")
say("eQTLGen 100kb cis: %d 行, %d 基因", nrow(eq), uniqueN(eq$gene))
genes <- c("PPIA","UTP20","RAB27A","MSN","EZH2","GNG12","ADGRE5","LRRK2","AQP1","ACY3","MYC","EZR","CCT3","SNX18","EIF2S1")
e <- eq[gene %in% genes]
setnames(e, c("chr.exposure","pos.exposure","effect_allele.exposure","other_allele.exposure","beta.exposure","se.exposure","eaf.exposure"),
            c("chr","pos","ea","oa","beta","se","eaf"), skip_absent = TRUE)
say("模型基因在 eQTLGen 的 cis 行数: %d（基因 %d 个：%s）", nrow(e), uniqueN(e$gene), paste(sort(unique(e$gene)), collapse = ", "))
say("（未出现在 eQTLGen 100kb cis 中的基因：%s）", paste(setdiff(genes, unique(e$gene)), collapse = ", "))
print(e[, .(nSNP = .N, chr = chr[1], start = min(pos), end = max(pos)), by = gene][order(gene)])

## ---------- 2) 结局：FinnGen R12 全量 ----------
fg <- fread("00_data/finngen_R12_C3_GBM_ASTROCYTOMA_EXALLC.gz", showProgress = FALSE)
setnames(fg, gsub("^#", "", names(fg)))
say("FinnGen 行数: %d；列: %s", nrow(fg), paste(names(fg), collapse = ", "))
fg <- fg[chrom %in% unique(e$chr)]
say("结局端仅保留相关染色体后: %d 行", nrow(fg))
NCASE <- 544; NCTRL <- 378749
say("结局：%d 例 / %d 对照（%s）", NCASE, NCTRL, "FinnGen R12 C3_GBM_ASTROCYTOMA_EXALLC")

## ---------- 3) 逐基因 coloc（全窗，按 chr:pos + 等位基因一致） ----------
res <- list()
for (g in sort(unique(e$gene))) {
  ee <- e[gene == g]
  ch <- ee$chr[1]; lo <- min(ee$pos) - 2e5; hi <- max(ee$pos) + 2e5
  ff <- fg[chrom == ch & pos >= lo & pos <= hi]
  if (nrow(ff) == 0) { say("%s：结局端该窗口内无 SNP，跳过", g); next }
  ## 注意：eQTLGen(2019) 为 GRCh37、FinnGen R12 为 GRCh38，坐标不可直接比；两边都带 rsID，故按 rsID 合并
  ff[, rsid1 := sub(",.*$", "", rsids)]
  ff <- ff[!is.na(rsid1) & rsid1 != ""]
  m <- merge(ee[, .(SNP, ea, oa, beta, se, eaf)],
             ff[, .(rsid1, ref, alt, beta_o = beta, se_o = sebeta, af = af_alt)],
             by.x = "SNP", by.y = "rsid1")
  m <- m[toupper(ea) == toupper(alt) & toupper(oa) == toupper(ref)]     # 等位基因方向一致
  m <- m[!is.na(beta) & !is.na(se) & !is.na(beta_o) & !is.na(se_o) & is.finite(beta) & is.finite(se) &
         is.finite(beta_o) & is.finite(se_o) & se > 0 & se_o > 0]
  m <- m[!is.na(eaf) & eaf > 0 & eaf < 1]
  m <- m[, maf := pmin(eaf, 1 - eaf)][maf > 0.001 & maf < 0.5]
  m <- m[!duplicated(SNP)]
  if (nrow(m) < 30) { say("%s：等位基因一致且可用的 SNP 仅 %d 个，跳过 coloc", g, nrow(m)); next }
  say("%s：待 coloc SNP %d 个（MAF 中位 %.3f）", g, nrow(m), median(m$maf))
  r <- tryCatch(coloc.abf(dataset1 = list(beta = m$beta, varbeta = m$se^2, snp = m$SNP,
                                          N = 30895, type = "quant", MAF = pmin(m$eaf, 1 - m$eaf)),
                          dataset2 = list(beta = m$beta_o, varbeta = m$se_o^2, snp = m$SNP,
                                          N = NCASE + NCTRL, type = "cc", s = NCASE / (NCASE + NCTRL))),
                 error = function(e) { say("%s：coloc 报错 → %s", g, conditionMessage(e)); NULL })
  if (is.null(r)) next
  s <- r$summary
  res[[g]] <- data.frame(Gene = g, nSNP = nrow(m),
                         PP_H0 = s[["PP.H0.abf"]], PP_H1 = s[["PP.H1.abf"]], PP_H2 = s[["PP.H2.abf"]],
                         PP_H3 = s[["PP.H3.abf"]], PP_H4 = s[["PP.H4.abf"]])
  say("%s：nSNP=%d，PP.H4=%.3f，PP.H3=%.3f，PP.H0=%.3f", g, nrow(m), s[["PP.H4.abf"]], s[["PP.H3.abf"]], s[["PP.H0.abf"]])
}
out <- rbindlist(res)
if (nrow(out) > 0) setorder(out, -PP_H4) else say("没有任何基因通过 nSNP>=30 的门槛")
fwrite(out, file.path(OUTD, "blood_fullcis_coloc.txt"), sep = "\t")
if (nrow(out) > 0) say("\n（血液）全窗共定位：%d 个基因，PP.H4 最大 %.3f（%s）", nrow(out), max(out$PP_H4), out$Gene[which.max(out$PP_H4)])
writeLines(c(sprintf("== 血液端全 cis 窗共定位（eQTLGen 100kb × FinnGen R12，%s）==", format(Sys.time(), "%Y-%m-%d %H:%M")), LOG),
           file.path(OUTD, "blood_fullcis_coloc_log.txt"))
print(out, digits = 3, row.names = FALSE)