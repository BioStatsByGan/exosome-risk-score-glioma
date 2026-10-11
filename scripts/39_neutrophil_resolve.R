## 步骤 39：中性粒细胞三方法冲突的系统性验证
## 1) 组成性检验：CIBERSORT 比例 CLR 变换后重测
## 2) 分层检验：WHO IV 级 / IDH 野生型内部直接比较
## 3) 标志基因检验：中性粒特异标志基因平均表达（不依赖解卷积）
## 4) 一致性矩阵：五种估计两两 Spearman 相关
## 输出：04_results/39.neutrophil/neutrophil_resolution.txt, method_cor_matrix.txt, FigS5_neutrophil_methods.pdf
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(ggplot2); library(MCPcounter); library(patchwork)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
setwd(ROOT); source("scripts/plot_style.R")
OUTD <- "04_results/39.neutrophil"; dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
sink_out <- character()
say <- function(...) { s <- sprintf(...); cat(s, "\n"); sink_out <<- c(sink_out, s) }

## ---------- 数据 ----------
cb <- fread("04_results/27.CIBERSORT/CIBERSORT-Results.txt"); setnames(cb, 1, "id")
frac_cols <- setdiff(names(cb), c("id", "P-value", "Correlation", "RMSE"))
F <- as.matrix(cb[, ..frac_cols]); rownames(F) <- cb$id

ex <- fread("02_inputs/TCGA_glioma_tumor_log2tpm.csv")
M <- as.matrix(ex[, -1]); rownames(M) <- ex[[1]]; M <- M[!duplicated(rownames(M)), ]

rk <- fread("04_results/20.SHAP/risk.TCGA.txt"); rk <- rk[!is.na(futime) & !is.na(fustat)]
cl <- fread("02_inputs/clinical_glioma.txt")
rk <- merge(rk[, .(id, riskScore, Risk)], cl[, .(patient, grade, idh)],
            by.x = "id", by.y = "patient", all.x = TRUE)
rk[, Risk := factor(ifelse(Risk == "high", "High", "Low"), levels = c("Low", "High"))]
sam <- Reduce(intersect, list(rk$id, rownames(F), colnames(M)))
F <- F[sam, ]; M <- M[, sam]; rk <- rk[match(sam, rk$id)]
say("样本 %d（高危 %d，低危 %d）", length(sam), sum(rk$Risk == "High"), sum(rk$Risk == "Low"))

## ---------- 1) CLR 变换 ----------
clr <- function(x) {
  n <- length(x); d <- 0.65 / (n^2); nz <- sum(x == 0)
  xz <- ifelse(x == 0, d, x * (1 - nz * d))
  log(xz) - mean(log(xz))
}
Fclr <- t(apply(F, 1, clr)); colnames(Fclr) <- colnames(F)

## ---------- 2) 各方法中性粒 ----------
mcp_raw <- fread("04_results/38.MCPcounter/mcp_scores.txt")
mcp <- data.table(id = names(mcp_raw)[-1], MCPcounter = as.numeric(mcp_raw[Population == "Neutrophils", -1]))
ss_raw <- fread("04_results/29.immFunction/ssgsea_scores.txt")
setnames(ss_raw, 1, "Feature")
say("ssGSEA 特征数 %d，含中性粒行: %s", nrow(ss_raw), "Neutrophils" %in% ss_raw$Feature)
ss <- data.table(id = names(ss_raw)[-1], ssGSEA = as.numeric(ss_raw[Feature == "Neutrophils", -1]))
dt <- data.table(id = rownames(F), CIBERSORT_frac = F[, "Neutrophils"], CIBERSORT_clr = Fclr[, "Neutrophils"])
dt <- merge(dt, mcp, by = "id", all.x = TRUE)
dt <- merge(dt, ss, by = "id", all.x = TRUE)

## ---------- 3) 标志基因平均表达 ----------
genes <- tryCatch({
  e <- new.env(); data("MCPcounter.genes", package = "MCPcounter", envir = e)
  unlist(e$MCPcounter.genes[["Neutrophils"]])
}, error = function(e) NULL)
if (is.null(genes)) genes <- c("FCGR3B", "CSF3R", "CXCR2", "FPR1", "AQP9", "S100A12", "MMP9")
genes <- intersect(genes, rownames(M))
say("中性粒标志基因（%d 个）: %s", length(genes), paste(genes, collapse = ", "))
dt[, MarkerScore := colMeans(M[genes, , drop = FALSE])[match(id, colnames(M))]]
dt <- merge(dt, rk[, .(id, Risk, grade, idh, riskScore)], by = "id")

## ---------- 4) 检验 ----------
test_all <- function(v) {
  x <- dt[[v]]
  w <- wilcox.test(x ~ dt$Risk)
  d <- median(x[dt$Risk == "High"]) - median(x[dt$Risk == "Low"])
  pvg <- tryCatch(summary(lm(x ~ dt$Risk + factor(dt$grade)))$coefficients["dt$RiskHigh", 4],
                  error = function(e) NA_real_)
  pv <- tryCatch(summary(lm(x ~ dt$Risk + factor(dt$grade) + factor(dt$idh)))$coefficients["dt$RiskHigh", 4],
                 error = function(e) NA_real_)
  ## 用 which() 取子集，避免逻辑下标里的 NA 混进来（否则 median() 返回 NA，输出表会出现空列）
  g4 <- which(!is.na(dt$grade) & dt$grade == "G4")
  iw <- which(!is.na(dt$idh) & dt$idh == "WT")
  p4 <- tryCatch(wilcox.test(x[g4] ~ droplevels(dt$Risk[g4]))$p.value, error = function(e) NA_real_)
  pidh <- tryCatch(wilcox.test(x[iw] ~ droplevels(dt$Risk[iw]))$p.value, error = function(e) NA_real_)
  dir4 <- median(x[g4][dt$Risk[g4] == "High"]) - median(x[g4][dt$Risk[g4] == "Low"])
  data.table(Method = v, delta = d, P_raw = w$p.value, P_adj_grade = pvg, P_adj_grade_idh = pv,
             delta_G4 = dir4, P_within_G4 = p4, P_within_IDHwt = pidh,
             rho_risk = cor(x, dt$riskScore, method = "spearman"))
}
res <- rbindlist(lapply(c("CIBERSORT_frac", "CIBERSORT_clr", "MCPcounter", "ssGSEA", "MarkerScore"), test_all))
res[, `:=`(delta = round(delta, 3), delta_G4 = round(delta_G4, 3), rho_risk = round(rho_risk, 3))]
say("\n== 五种中性粒估计的组间差异 ==")
print(res)
# 第 6 个估计：脑特异署名集（含小胶质）——读患者级评分，用同一套检验，保证六种口径一致
bs <- tryCatch(fread("04_results/audit/brain_aware_neutrophil_scores.txt"), error = function(e) NULL)
if (!is.null(bs) && nrow(bs) > 0) {
  setnames(bs, c("patient", "BrainAware_ssGSEA"))
  dt <- merge(dt, bs, by.x = "id", by.y = "patient", all.x = TRUE)
  res <- rbindlist(list(res, test_all("BrainAware_ssGSEA")))
} else {
  say("注意：未找到 brain_aware_neutrophil_scores.txt，图S5 只画 5 种估计")
}
fwrite(res, file.path(OUTD, "neutrophil_resolution.txt"), sep = "\t")
say("\n== 统一对比表（含脑特异署名集）==")
for (i in seq_len(nrow(res))) say("%-20s delta=%7.3f  P_raw=%.3g  P_adj(分级)=%.3g  P_adj(分级+IDH)=%.3g  G4内=%.3g  IDHwt内=%.3g",
    res$Method[i], res$delta[i], res$P_raw[i], res$P_adj_grade[i], res$P_adj_grade_idh[i],
    res$P_within_G4[i], res$P_within_IDHwt[i])

## ---------- 5) 一致性矩阵 ----------
mat <- dt[, .(CIBERSORT_frac, CIBERSORT_clr, MCPcounter, ssGSEA, MarkerScore, BrainAware_ssGSEA)]
cm <- cor(as.matrix(mat), method = "spearman", use = "pairwise.complete.obs")
write.table(round(cm, 3), file.path(OUTD, "method_cor_matrix.txt"), sep = "\t", quote = FALSE)
say("\n== 方法间 Spearman 相关 =="); print(round(cm, 3))

## ---------- 6) 结论判定 ----------
f <- res[Method == "CIBERSORT_frac"]; c2 <- res[Method == "CIBERSORT_clr"]
say("\n== 判定 ==")
say("CIBERSORT 比例：原始 P=%.3g，校正分级 P=%.3g，校正分级+IDH P=%.3g", f$P_raw, f$P_adj_grade, f$P_adj_grade_idh)
say("CIBERSORT CLR：原始 P=%.3g，校正分级 P=%.3g，校正分级+IDH P=%.3g", c2$P_raw, c2$P_adj_grade, c2$P_adj_grade_idh)
say("MCP-counter：分级+IDH P=%.3g | ssGSEA：%.3g | 标志基因：%.3g | 脑特异署名：%.3g",
    res[Method == "MCPcounter"]$P_adj_grade_idh, res[Method == "ssGSEA"]$P_adj_grade_idh,
    res[Method == "MarkerScore"]$P_adj_grade_idh, res[Method == "BrainAware_ssGSEA"]$P_adj_grade_idh)
say("G4 内部：CIBERSORTCLR P=%.3g | MCP P=%.3g | ssGSEA P=%.3g | 标志基因 P=%.3g",
    c2$P_within_G4, res[Method == "MCPcounter"]$P_within_G4, res[Method == "ssGSEA"]$P_within_G4,
    res[Method == "MarkerScore"]$P_within_G4)
writeLines(sink_out, file.path(OUTD, "neutrophil_resolution_report.txt"))

## ---------- FigS5 ----------
lg <- melt(dt[, .(id, Risk, CIBERSORT_frac, CIBERSORT_clr, MCPcounter, ssGSEA, MarkerScore, BrainAware_ssGSEA)],
           id.vars = c("id", "Risk"), variable.name = "Method", value.name = "value")
lg[, Method := factor(Method,
                      levels = c("CIBERSORT_frac", "CIBERSORT_clr", "MCPcounter", "ssGSEA", "MarkerScore", "BrainAware_ssGSEA"),
                      labels = c("CIBERSORT fraction", "CIBERSORT CLR", "MCP-counter", "ssGSEA signature", "Marker genes",
                                 "Brain-aware ssGSEA"))]
p1 <- ggplot(lg, aes(Method, value, fill = Risk)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.3, width = 0.7) +
  scale_fill_manual(values = c(Low = OI[1], High = OI[2]), name = "") +
  facet_wrap(~Method, scales = "free_y", nrow = 1, strip.position = "bottom") +
  labs(x = NULL, y = "Neutrophil estimate (scaled within method)") +
  theme_pub() + theme(legend.position = "top", axis.text.x = element_blank(), axis.ticks.x = element_blank())
cmdf <- data.table(expand.grid(Var1 = rownames(cm), Var2 = colnames(cm), stringsAsFactors = FALSE))
cmdf[, N := cm[cbind(match(Var1, rownames(cm)), match(Var2, colnames(cm)))]]
p2 <- ggplot(cmdf, aes(Var1, Var2, fill = N)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", N)), size = 2.4) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, limits = c(-1, 1)) +
  labs(x = NULL, y = NULL, fill = "rho") + theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6.4), axis.text.y = element_text(size = 6.4))
save_pdf(p1 / p2 + plot_layout(heights = c(1.1, 1)), "03_figures/FigS5_neutrophil_methods.pdf", 11.4, 6.8)
say("\nFigS5 完成")
