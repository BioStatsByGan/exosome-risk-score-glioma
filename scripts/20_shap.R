## 步骤 20：LASSO 降维 + stepwise 多因素 Cox + SHAP（kernelshap/shapviz）
## 输入：04_results/19.uniCox/TCGA.uniSigExp.txt（单因素显著基因）
## 输出：04_results/20.SHAP/multiCox.txt、risk.TCGA.txt、risk.TCGA_test.txt（内部留出验证）
##       03_figures/Fig2{A,B,C,D}_shap*.pdf（课程原绘图代码：shapviz bar/bee/waterfall/force）
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({
  library(data.table); library(glmnet); library(survival)
  library(kernelshap); library(shapviz); library(ggplot2)
})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
setwd(ROOT)
dir.create("04_results/20.SHAP", recursive = TRUE, showWarnings = FALSE)

rt <- read.table("04_results/19.uniCox/TCGA.uniSigExp.txt", header = TRUE, sep = "\t",
                 check.names = FALSE, row.names = 1)
rt <- na.omit(rt)
cat("单因素显著基因进入建模:", ncol(rt) - 2, " 样本:", nrow(rt), "\n")

## ---------- 1. LASSO-Cox 降维 ----------
set.seed(20260919)
x <- as.matrix(rt[, -(1:2)]); y <- Surv(rt$futime, rt$fustat)
cvfit <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10)
cf <- as.matrix(coef(cvfit, s = "lambda.min"))
sel <- rownames(cf)[cf[, 1] != 0]
cat("LASSO 选出基因:", length(sel), "个:", paste(sel, collapse = ", "), "\n")
fwrite(data.frame(gene = sel, coef = cf[sel, 1]), "04_results/20.SHAP/lasso_genes.txt", sep = "\t", quote = FALSE)

## ---------- 2. stepwise 多因素 Cox ----------
d <- data.frame(futime = rt$futime, fustat = rt$fustat, rt[, sel, drop = FALSE])
multiCox <- coxph(Surv(futime, fustat) ~ ., data = d)
multiCox <- step(multiCox, direction = "both", trace = 0)
s <- summary(multiCox)
outMulti <- data.frame(id = rownames(s$coefficients), coef = s$coefficients[, "coef"],
                       HR = s$conf.int[, "exp(coef)"], HR.95L = s$conf.int[, "lower .95"],
                       HR.95H = s$conf.int[, "upper .95"], pvalue = s$coefficients[, "Pr(>|z|)"])
outMulti <- outMulti[order(outMulti$HR, decreasing = TRUE), ]
fwrite(outMulti, "04_results/20.SHAP/multiCox.txt", sep = "\t", quote = FALSE)
modelGene <- outMulti$id
cat("stepwise 最终模型基因:", length(modelGene), "个\n")
print(outMulti, row.names = FALSE, digits = 4)
cat(sprintf("模型 C-index (训练集): %.3f\n", s$concordance[1]))

## ---------- 3. SHAP ----------
X <- as.data.frame(d[, modelGene, drop = FALSE])
fit <- additive_shap(multiCox, X)
shp <- shapviz(fit, X_pred = X, X = X, interactions = TRUE)
imp <- sort(colMeans(abs(shp$S)), decreasing = TRUE)
cat("SHAP 平均|值| 排序:\n"); print(round(imp, 4))

sv_importance(shp, kind = "bar", show_numbers = TRUE) + theme_bw() + labs(title = NULL)
ggsave("03_figures/Fig2A_shap_bar.pdf", width = 6, height = max(3.5, 0.32 * length(imp) + 1.6), device = cairo_pdf)
sv_importance(shp, kind = "bee", show_numbers = TRUE) + theme_bw() + labs(title = NULL)
ggsave("03_figures/Fig2B_shap_bee.pdf", width = 7, height = max(3.5, 0.32 * length(imp) + 1.6), device = cairo_pdf)
p_wf <- sv_waterfall(shp, row_id = 4); p_wf <- p_wf + labs(title = NULL)
ggsave("03_figures/Fig2C_shap_waterfall.pdf", p_wf, width = 7.5, height = 5, device = cairo_pdf)
p_fc <- sv_force(shp, row_id = 4); p_fc <- p_fc + labs(title = NULL)
ggsave("03_figures/Fig2D_shap_force.pdf", p_fc, width = 9.5, height = 5, device = cairo_pdf)
cat("saved: Fig2A-D SHAP 四联图\n")

## ---------- 4. 风险评分（训练集全量 + 内部 7:3 留出验证）----------
cutoff <- median(fit$predictions)
writeLines(sprintf("%.6f", cutoff), "04_results/20.SHAP/cutoff.txt")
build <- function(score, dd) {
  data.frame(futime = dd$futime, fustat = dd$fustat, dd[, modelGene, drop = FALSE],
             riskScore = as.vector(score), Risk = ifelse(score > cutoff, "high", "low"))
}
outTab <- build(fit$predictions, d)
fwrite(cbind(id = rownames(outTab), outTab), "04_results/20.SHAP/risk.TCGA.txt", sep = "\t", quote = FALSE)

set.seed(42)
idx <- sample(nrow(d), round(0.3 * nrow(d)))
lp_all <- predict(multiCox, newdata = X, type = "lp")
outTest <- build(lp_all[idx], d[idx, ])
fwrite(cbind(id = rownames(outTest), outTest), "04_results/20.SHAP/risk.TCGA_test.txt", sep = "\t", quote = FALSE)
cat(sprintf("风险分组: 训练集 high=%d low=%d | 留出集 high=%d low=%d\n",
            sum(outTab$Risk == "high"), sum(outTab$Risk == "low"),
            sum(outTest$Risk == "high"), sum(outTest$Risk == "low")))

## 公式导出，供外部队列（CGGA）套用
writeLines(paste(modelGene, collapse = " + "), "04_results/20.SHAP/model_formula.txt")
cat(sprintf("模型公式已存; cutoff = %.6f\n", cutoff))
cat("STEP20 DONE\n")
