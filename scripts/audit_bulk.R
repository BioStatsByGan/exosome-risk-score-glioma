## 投稿前数据漏洞审计（bulk 层）：模型乐观度、混杂校正、离群值敏感性
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(survival); library(glmnet); library(ggplot2)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
OUT <- "04_results/audit"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
## 控制台输出同时落盘（C 指数增量、R2、TMB、免疫校正后 P 值均在此日志内）
sink(file.path(OUT, "audit_bulk_log.txt"), split = TRUE)
risk <- read.table("04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t",
                   check.names = FALSE, row.names = 1)
genes <- setdiff(colnames(risk), c("futime", "fustat", "riskScore", "Risk"))

cat("========== 1. 模型乐观度：10 折交叉验证 C 指数 ==========\n")
cat("训练集 C 指数（表观，有乐观偏倚）:",
    round(summary(coxph(Surv(futime, fustat) ~ riskScore, data = risk))$concordance[1], 3), "\n")
set.seed(42)
k <- 10; folds <- sample(rep(1:k, length.out = nrow(risk)))
cv <- numeric(k)
for (i in 1:k) {
  tr <- risk[folds != i, ]; te <- risk[folds == i, ]
  m <- coxph(Surv(futime, fustat) ~ ., data = tr[, c(genes, "futime", "fustat")])
  lp <- predict(m, newdata = te)
  cv[i] <- survConcordance(Surv(te$futime, te$fustat) ~ lp)$concordance
}
cat(sprintf("10 折交叉验证 C 指数（模型系数重新拟合）: %.3f (各折 %s)\n",
            mean(cv), paste(round(cv, 2), collapse = ", ")))
cat(sprintf("乐观度 = %.3f\n", summary(coxph(Surv(futime, fustat) ~ riskScore, data = risk))$concordance[1] - mean(cv)))

cat("\n========== 2. 与临床变量比增量价值 ==========\n")
cli <- read.table("02_inputs/clinical_glioma.txt", header = TRUE, sep = "\t", check.names = FALSE)
d <- merge(data.frame(patient = rownames(risk), risk, check.names = FALSE),
           cli[, c("patient", "age", "sex", "grade", "idh")], by = "patient")
d$Grade <- factor(d$grade); d$IDH <- factor(d$idh); d$Sex <- factor(d$sex)
m_cli <- coxph(Surv(futime, fustat) ~ age + Sex + Grade + IDH, data = d)
m_all <- coxph(Surv(futime, fustat) ~ riskScore + age + Sex + Grade + IDH, data = d)
c_cli <- summary(m_cli)$concordance[1]; c_all <- summary(m_all)$concordance[1]
cat(sprintf("临床变量模型 C 指数: %.3f | 加风险评分后: %.3f | 增量 %.3f\n", c_cli, c_all, c_all - c_cli))
cat(sprintf("似然比检验 P = %.3g\n", anova(m_cli, m_all)$`Pr(>|Chi|)`[2]))
r2 <- summary(lm(riskScore ~ as.numeric(Grade) + as.numeric(IDH), data = d))$r.squared
cat(sprintf("风险评分被分级与 IDH 解释的方差比例 R2 = %.3f（越高说明模型越接近分级代理）\n", r2))

cat("\n========== 3. TMB 分析：离群值与分级混杂 ==========\n")
tmb <- fread("04_results/31.riskTMB/tmb_risk.txt")
cat("TMB 最大值:", round(max(tmb$TMB), 1), "（超突变体）| >50 的样本数:", sum(tmb$TMB > 50), "\n")
cat(sprintf("全部样本: 高危 %.2f vs 低危 %.2f, P=%.3g\n",
            median(tmb$TMB[tmb$Risk == "High risk"]), median(tmb$TMB[tmb$Risk == "Low risk"]),
            wilcox.test(TMB ~ Risk, data = tmb)$p.value))
t2 <- tmb[TMB < 50]
cat(sprintf("去掉超突变体后: 高危 %.2f vs 低危 %.2f, P=%.3g\n",
            median(t2$TMB[t2$Risk == "High risk"]), median(t2$TMB[t2$Risk == "Low risk"]),
            wilcox.test(TMB ~ Risk, data = t2)$p.value))
dd <- merge(tmb, cli[, c("patient", "grade", "idh")], by = "patient")
dd$Grade <- factor(dd$grade)
cat("\n校正 WHO 分级后风险评分与 log2(TMB) 的关系:\n")
print(summary(lm(log2(TMB + 1) ~ riskScore + Grade, data = dd))$coefficients["riskScore", ], digits = 3)
cat("仅分级解释 log2(TMB) 的 R2 =", round(summary(lm(log2(TMB + 1) ~ Grade, data = dd))$r.squared, 3), "\n")
cat("分级 + 风险评分 R2 =", round(summary(lm(log2(TMB + 1) ~ riskScore + Grade, data = dd))$r.squared, 3), "\n")

cat("\n========== 4. 免疫浸润差异：校正分级后是否仍在 ==========\n")
imm <- read.table("04_results/27.CIBERSORT/CIBERSORT-Results.txt", header = TRUE, sep = "\t",
                  check.names = FALSE, row.names = 1)
rownames(imm) <- gsub("^(.*?)-(.*?)-(.*?)-.*", "\\1-\\2-\\3", rownames(imm))
cells <- setdiff(colnames(imm), c("P-value", "Correlation", "RMSE"))
di <- merge(data.frame(patient = rownames(imm), imm[, cells], check.names = FALSE),
            cli[, c("patient", "grade")], by = "patient")
di <- merge(di, data.frame(patient = rownames(risk), Risk = risk$Risk), by = "patient")
di$Grade <- factor(di$grade)
key <- c("Macrophages M2", "Macrophages M0", "Neutrophils", "Monocytes", "T cells CD8")
res <- do.call(rbind, lapply(key, function(cc) {
  f1 <- as.formula(sprintf("%s ~ Risk", paste0("`", cc, "`")))
  f2 <- as.formula(sprintf("%s ~ Risk + Grade", paste0("`", cc, "`")))
  p1 <- summary(lm(f1, data = di))$coefficients["Risklow", 4]
  s2 <- summary(lm(f2, data = di))$coefficients
  data.frame(Cell = cc, p_unadj = p1, p_adj_grade = s2["Risklow", 4],
             beta_unadj = summary(lm(f1, data = di))$coefficients["Risklow", 1],
             beta_adj = s2["Risklow", 1])
}))
print(res, row.names = FALSE, digits = 3)
fwrite(res, file.path(OUT, "audit_immune_grade_adjusted.txt"), sep = "\t", quote = FALSE)
sink()
cat("\nwritten:", file.path(OUT, "audit_bulk_log.txt"), "\n")
cat("AUDIT-BULK DONE\n")
