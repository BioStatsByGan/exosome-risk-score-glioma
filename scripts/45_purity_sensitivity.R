## 步骤 45：纯度校正敏感性分析（回应外部意见"纯度相关敏感性分析"）
## 纯度来源：① ESTIMATE（胶质瘤 gct，804 例，cos 公式换算）；② GBM/LGG clinicalSE 的 paper_ABSOLUTE.purity（按 12 位 barcode 匹配）
## 输出：04_results/audit/purity_table.txt、purity_sensitivity.txt、immune_all_cells_adjusted_purity.txt
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
OUTD <- "04_results/audit"; dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
LOG <- character(); say <- function(...) { s <- sprintf(...); cat(s, "\n"); LOG <<- c(LOG, s); invisible(s) }

## ---------- 纯度 ----------
est_file <- "F:/0000_Bioinfo/219_棕榈酰化_TCGA_胶质瘤/31.estimate/estimateScore.gct"
g <- readLines(est_file, n = 3)
hdr <- strsplit(g[3], "\t")[[1]]
samples <- gsub("\\.", "-", gsub('"', "", hdr[-(1:2)]))
sc <- read.table(est_file, skip = 3, sep = "\t", header = FALSE, quote = "", check.names = FALSE)
rownames(sc) <- sc[[1]]
est <- as.numeric(sc["ESTIMATEScore", -(1:2)])
names(est) <- samples
pur_est <- cos(0.6049872018 + 0.0001467884 * est)      # Yoshihara 2013 的 ESTIMATE 纯度换算
say("ESTIMATE 纯度：%d 例（中位 %.3f，范围 %.3f–%.3f）", length(pur_est), median(pur_est),
    min(pur_est), max(pur_est))

pur_abs <- numeric(0)
for (co in c("GBM", "LGG")) {
  f <- sprintf("F:/0000_Bioinfo/223_外泌体_泛癌预后/01_pancancer/%s_clinicalSE.csv", co)
  d <- fread(f, data.table = FALSE, check.names = FALSE)
  col <- grep("purity", names(d), ignore.case = TRUE, value = TRUE)[1]
  v <- suppressWarnings(as.numeric(d[[col]]))
  id <- substr(d[[1]], 1, 12)
  keep <- !is.na(v)
  pur_abs <- c(pur_abs, setNames(v[keep], id[keep]))
  say("%s：%s 可用 %d 例", co, col, sum(keep))
}
pur_abs <- pur_abs[!duplicated(names(pur_abs))]
say("ABSOLUTE 纯度（唯一样本）：%d 例", length(pur_abs))

rk <- read.table("04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
cl <- read.table("02_inputs/clinical_glioma.txt", header = TRUE, sep = "\t", check.names = FALSE)
ids <- rownames(rk)
say("风险表 %d 例：ESTIMATE 覆盖 %d，ABSOLUTE 覆盖 %d", length(ids), sum(ids %in% names(pur_est)), sum(ids %in% names(pur_abs)))
pt <- data.frame(patient = ids,
                 purity_estimate = pur_est[ids],
                 purity_absolute = pur_abs[ids])
fwrite(pt, file.path(OUTD, "purity_table.txt"), sep = "\t")

## ---------- 敏感性分析 ----------
d <- merge(data.frame(patient = ids, Risk = rk$Risk, riskScore = rk$riskScore), pt, by = "patient")
d <- merge(d, cl[, c("patient", "grade", "idh")], by = "patient")
d$Grade <- factor(d$grade); d$IDH <- factor(d$idh)
d$Risk <- factor(d$Risk, levels = c("low", "high"))
p_ok <- !is.na(d$purity_estimate)
say("可用于纯度校正的样本：%d（ESTIMATE 完整）", sum(p_ok))

## (a) TMB
tmb <- fread("04_results/31.riskTMB/tmb_risk.txt")
td <- merge(tmb, d[, c("patient", "purity_estimate", "grade", "idh")], by = "patient")
td$Risk <- factor(ifelse(grepl("^[Hh]", td$Risk), "high", "low"), levels = c("low", "high"))   # tmb_risk.txt 里是 "High risk"/"Low risk" 标签
m1 <- summary(lm(log2(TMB + 1) ~ Risk + Grade + IDH, data = transform(td, Grade = factor(grade), IDH = factor(idh))))$coefficients
m2 <- summary(lm(log2(TMB + 1) ~ Risk + Grade + IDH + purity_estimate, data = transform(td, Grade = factor(grade), IDH = factor(idh))))$coefficients
say("\n(a) TMB（log2）~ Risk + 分级 + IDH：Riskhigh 系数 %.4f（P=%.3g）；再加纯度：%.4f（P=%.3g）",
    m1["Riskhigh", 1], m1["Riskhigh", 4], m2["Riskhigh", 1], m2["Riskhigh", 4])

## (b) 22 类 CIBERSORT
imm <- read.table("04_results/27.CIBERSORT/CIBERSORT-Results.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
rownames(imm) <- substr(rownames(imm), 1, 12)
cells <- setdiff(colnames(imm), c("P-value", "Correlation", "RMSE"))
di <- merge(data.frame(patient = rownames(imm), imm[, cells], check.names = FALSE), d, by = "patient")
res <- do.call(rbind, lapply(cells, function(cc) {
  f <- as.formula(sprintf("`%s` ~ Risk + Grade + IDH + purity_estimate", cc))
  s <- summary(lm(f, data = di))$coefficients
  data.frame(Cell = cc, beta_grade_idh_purity = s["Riskhigh", 1], p_grade_idh_purity = s["Riskhigh", 4])
}))
res$FDR <- p.adjust(res$p_grade_idh_purity, "BH")
res <- res[order(res$p_grade_idh_purity), ]
fwrite(res, file.path(OUTD, "immune_all_cells_adjusted_purity.txt"), sep = "\t")
say("(b) 22 类 CIBERSORT 再加纯度协变量后 P<0.05 的细胞数：%d（校正分级+IDH 时为 5）", sum(res$p_grade_idh_purity < 0.05))
print(head(res, 6), digits = 3, row.names = FALSE)

## (c) Hallmark 关键通路
hs <- fread("04_results/37.hallmark/hallmark_scores.txt")
setnames(hs, 1, "Pathway")
key <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
         "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_ANGIOGENESIS", "HALLMARK_INFLAMMATORY_RESPONSE",
         "HALLMARK_IL6_JAK_STAT3_SIGNALING")
hh <- t(as.matrix(hs[Pathway %in% key, -1]))
colnames(hh) <- hs[Pathway %in% key, Pathway]
hd <- merge(data.frame(patient = rownames(hh), hh, check.names = FALSE), d, by = "patient")
say("\n(c) Hallmark 关键通路：Riskhigh 系数（校正分级+IDH → 再加纯度）")
for (k in key) {
  f1 <- as.formula(sprintf("`%s` ~ Risk + Grade + IDH", k))
  f2 <- as.formula(sprintf("`%s` ~ Risk + Grade + IDH + purity_estimate", k))
  s1 <- summary(lm(f1, data = hd))$coefficients["Riskhigh", ]; s2 <- summary(lm(f2, data = hd))$coefficients["Riskhigh", ]
  say("    %-45s %.3f (P=%.2g) → %.3f (P=%.2g)", sub("HALLMARK_", "", k), s1[1], s1[4], s2[1], s2[4])
}

## (d) 药敏（替莫唑胺 / 卡莫司汀 / 顺铂）
dp <- fread("04_results/34.oncoPredict/calcPhenotype_Output/DrugPredictions.csv", data.table = FALSE)
rownames(dp) <- substr(dp[[1]], 1, 12); dp <- dp[, -1]
say("\n(d) 药敏 log2(IC50)：Riskhigh 系数（校正分级+IDH → 再加纯度）")
for (dg in c("Temozolomide_1375", "Carmustine_1807", "Cisplatin_1005")) {
  if (!dg %in% colnames(dp)) { say("    缺 %s", dg); next }
  dd <- merge(data.frame(patient = rownames(dp), ic50 = log2(dp[[dg]] + 1), check.names = FALSE), d, by = "patient")
  s1 <- summary(lm(ic50 ~ Risk + Grade + IDH, data = dd))$coefficients["Riskhigh", ]
  s2 <- summary(lm(ic50 ~ Risk + Grade + IDH + purity_estimate, data = dd))$coefficients["Riskhigh", ]
  say("    %-22s %.3f (P=%.2g) → %.3f (P=%.2g)", dg, s1[1], s1[4], s2[1], s2[4])
}

## ---------- 用 ABSOLUTE 纯度（DNA 来源，n≈675）复核关键结果，避免"用 ESTIMATE 校正免疫比例"的循环性 ----------
say("\n===== ABSOLUTE 纯度（DNA，独立于表达）复核 =====")
dA <- d[!is.na(d$purity_absolute), ]
say("可用样本 %d（ESTIMATE 为 %d）", nrow(dA), sum(p_ok))
## TMB
tdA <- merge(tmb, dA[, c("patient", "purity_absolute", "grade", "idh")], by = "patient")
tdA$Risk <- factor(ifelse(grepl("^[Hh]", tdA$Risk), "high", "low"), levels = c("low", "high"))
sA <- summary(lm(log2(TMB + 1) ~ Risk + Grade + IDH + purity_absolute,
                 data = transform(tdA, Grade = factor(grade), IDH = factor(idh))))$coefficients["Riskhigh", ]
say("(a) TMB + 分级 + IDH + ABSOLUTE 纯度：系数 %.4f（P=%.3g）", sA[1], sA[4])
## CIBERSORT
diA <- merge(data.frame(patient = rownames(imm), imm[, cells], check.names = FALSE), dA, by = "patient")
resA <- do.call(rbind, lapply(cells, function(cc) {
  s <- summary(lm(as.formula(sprintf("`%s` ~ Risk + Grade + IDH + purity_absolute", cc)), data = diA))$coefficients
  data.frame(Cell = cc, beta_abs = s["Riskhigh", 1], p_abs = s["Riskhigh", 4])
}))
resA$FDR <- p.adjust(resA$p_abs, "BH")
mrg <- merge(res, resA, by = "Cell")
fwrite(mrg, file.path(OUTD, "immune_all_cells_adjusted_purity.txt"), sep = "\t")
say("(b) 22 类 CIBERSORT：再加 ESTIMATE 纯度后 %d 类 P<0.05；再加 ABSOLUTE 纯度后 %d 类 P<0.05（校正分级+IDH 为 5 类）",
    sum(res$p_grade_idh_purity < 0.05), sum(resA$p_abs < 0.05))
print(head(mrg[order(mrg$p_abs), c("Cell", "p_grade_idh_purity", "p_abs")], 6), digits = 3, row.names = FALSE)
## Hallmark
hdA <- merge(data.frame(patient = rownames(hh), hh, check.names = FALSE), dA, by = "patient")
say("(c) Hallmark 关键通路 + ABSOLUTE 纯度")
for (k in key) {
  s <- summary(lm(as.formula(sprintf("`%s` ~ Risk + Grade + IDH + purity_absolute", k)), data = hdA))$coefficients["Riskhigh", ]
  say("    %-45s %.3f (P=%.2g)", sub("HALLMARK_", "", k), s[1], s[4])
}
## 药敏
for (dg in c("Temozolomide_1375", "Carmustine_1807", "Cisplatin_1005")) {
  ddA <- merge(data.frame(patient = rownames(dp), ic50 = log2(dp[[dg]] + 1), check.names = FALSE), dA, by = "patient")
  s <- summary(lm(ic50 ~ Risk + Grade + IDH + purity_absolute, data = ddA))$coefficients["Riskhigh", ]
  say("(d) %-22s %.3f (P=%.2g)", dg, s[1], s[4])
}

writeLines(c(sprintf("== 纯度校正敏感性分析（%s）==", format(Sys.time(), "%Y-%m-%d %H:%M")), LOG),
           file.path(OUTD, "purity_sensitivity.txt"))
cat("\nwritten:", file.path(OUTD, "purity_sensitivity.txt"), "\n")