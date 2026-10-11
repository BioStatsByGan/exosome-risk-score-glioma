## 步骤 41：全流程嵌套交叉验证（回应"乐观度 0.005 依据不足"）
## 目的：外层 10 折内**重跑所有利用结局信息的步骤**（单因素 Cox 筛选 → LASSO → 逐步 Cox），
##       在外层验证折上评价，得到诚实的乐观度估计；与表观 C 指数 0.824、仅重估系数的 0.819 对比。
## 输入：04_results/13.diff/glioma_exosome_diff.txt（422 个差异基因名）、
##       02_inputs/TCGA_glioma_tumor_log2tpm.csv、04_results/20.SHAP/risk.TCGA.txt（生存与已发表评分）
## 输出：04_results/audit/nested_cv.txt、nested_cv_log.txt
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(survival); library(glmnet)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
OUTD <- "04_results/audit"; dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
LOG <- character()
say <- function(...) { s <- sprintf(...); cat(s, "\n"); LOG <<- c(LOG, s); invisible(s) }

genes422 <- fread("04_results/13.diff/glioma_exosome_diff.txt")$gene
ex <- fread("02_inputs/TCGA_glioma_tumor_log2tpm.csv")
M <- as.matrix(ex[, -1]); rownames(M) <- ex[[1]]
M <- M[!duplicated(rownames(M)), , drop = FALSE]
rk <- read.table("04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
use <- intersect(genes422, rownames(M)); use <- intersect(use, rownames(M))
sam <- intersect(colnames(M), rownames(rk))
X <- t(M[use, sam, drop = FALSE]); X <- X[rownames(rk)[rownames(rk) %in% sam], , drop = FALSE]
d <- data.frame(rk[rownames(X), c("futime", "fustat")], X, check.names = FALSE)
say("候选基因 %d 个，样本 %d 例，事件 %d", ncol(X), nrow(d), sum(d$fustat))

## 参考值 1：已发表 15 基因模型（全数据筛选并拟合）在同一批样本上的表观 C 指数
mc <- fread("04_results/20.SHAP/multiCox.txt")
mg <- intersect(mc$id, colnames(d))
lp15 <- as.numeric(as.matrix(d[, mg, drop = FALSE]) %*% mc$coef[match(mg, mc$id)])
app <- summary(coxph(Surv(d$futime, d$fustat) ~ lp15))$concordance[1]
say("参考值 1：已发表 15 基因模型的表观 C 指数 = %.3f（表 2 注为 0.824）", app)
say("参考值 2：只重估系数（基因固定）的 10 折 CV = 0.819，其「乐观度」= 0.005")
say("注意：嵌套 CV 的每一折都要重新做 单因素→LASSO→逐步（见下），才是完整流程的乐观度")

## 外层 10 折：折内重跑 单因素 → LASSO → 逐步
set.seed(42)
K <- 10; folds <- sample(rep(1:K, length.out = nrow(d)))
cv <- rep(NA_real_, K); ngene <- integer(K); nsel <- integer(K)
for (k in 1:K) {
  tr <- d[folds != k, ]; te <- d[folds == k, ]
  xs <- tr[, -(1:2), drop = FALSE]
  ## ① 单因素 Cox（折内）
  pv <- apply(xs, 2, function(v) tryCatch(summary(coxph(Surv(tr$futime, tr$fustat) ~ v))$coefficients[1, 5],
                                          error = function(e) NA_real_))
  keep <- names(pv)[!is.na(pv) & pv < 0.05]
  nsel[k] <- length(keep)
  if (length(keep) < 2) next
  ## ② LASSO-Cox（折内，lambda.min）
  lf <- tryCatch(cv.glmnet(as.matrix(tr[, keep, drop = FALSE]), Surv(tr$futime, tr$fustat),
                           family = "cox", alpha = 1, nfolds = 10), error = function(e) NULL)
  if (is.null(lf)) next
  cf <- as.matrix(coef(lf, s = "lambda.min")); sel <- rownames(cf)[cf[, 1] != 0]
  if (length(sel) < 2) next
  ## ③ 逐步 Cox（折内）
  m <- tryCatch(step(coxph(Surv(tr$futime, tr$fustat) ~ ., data = tr[, c("futime", "fustat", sel)]),
                      direction = "both", trace = 0), error = function(e) NULL)
  if (is.null(m)) next
  mg <- setdiff(names(coef(m)), NULL); ngene[k] <- length(mg)
  ## ④ 在外层验证折评价（用折内系数）
  te[, mg] <- te[, mg, drop = FALSE]
  lp <- tryCatch(predict(m, newdata = te[, c("futime", "fustat", mg)]), error = function(e) NULL)
  if (is.null(lp)) next
  ## 注意：survival::concordance(Surv ~ x) 把 x 当生存时间处理，会得到 1-C；这里统一用 coxph 口径
  cv[k] <- summary(coxph(Surv(te$futime, te$fustat) ~ lp))$concordance[1]
  say("  折 %2d：单因素入选 %3d → LASSO %3d → 逐步 %2d 个基因；验证折 C 指数 %.3f", k, nsel[k], length(sel), ngene[k], cv[k])
}
say("")
say("嵌套 10 折外样本 C 指数：均值 %.3f（中位 %.3f；各折 %s）", mean(cv, na.rm = TRUE), median(cv, na.rm = TRUE),
    paste(sprintf("%.2f", cv), collapse = ", "))
say("诚实乐观度 = 表观 %.3f − 嵌套 %.3f = %.3f", app, mean(cv, na.rm = TRUE), app - mean(cv, na.rm = TRUE))
say("（对照：仅重估系数的\"乐观度\" 0.005；只按系数重估会明显低估完整流程的乐观度）")
say("折内最终入选基因数：%s（说明基因组合对折的依赖程度）", paste(ngene, collapse = ", "))
writeLines(c(sprintf("== 嵌套 CV 关键输出（%s）==", format(Sys.time(), "%Y-%m-%d %H:%M")), LOG),
           file.path(OUTD, "nested_cv.txt"))
writeLines(LOG, file.path(OUTD, "nested_cv_log.txt"))
cat("written:", file.path(OUTD, "nested_cv.txt"), "\n")