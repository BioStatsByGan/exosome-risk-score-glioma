## 步骤 42：两个敏感性分析
## (A) 完全跳过「TCGA 肿瘤 vs GTEx 正常」差异筛选，直接在 TCGA 生存数据里走 单因素→LASSO→逐步，
##     看模型与已发表 15 基因评分是否一致（回应审稿意见第 1 条）。
## (B) 13 基因缩减模型（CGGA 缺失 ADGRE5/AQP1）在 TCGA 内的表现，用于说明外部验证所用版本的代价（回应第 3 条）。
## 输入：02_inputs/exosome_genes.txt、02_inputs/TCGA_glioma_tumor_log2tpm.csv、
##       04_results/20.SHAP/{risk.TCGA.txt, multiCox.txt}、04_results/13.diff/glioma_exosome_diff.txt
## 输出：04_results/audit/screen_free_model.txt / _log.txt
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(survival); library(glmnet)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
OUTD <- "04_results/audit"; dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
LOG <- character(); say <- function(...) { s <- sprintf(...); cat(s, "\n"); LOG <<- c(LOG, s); invisible(s) }

## ---------- 数据 ----------
rk <- read.table("04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
ex <- fread("02_inputs/TCGA_glioma_tumor_log2tpm.csv"); M <- as.matrix(ex[, -1]); rownames(M) <- ex[[1]]
M <- M[!duplicated(rownames(M)), , drop = FALSE]
eg <- unique(readLines("02_inputs/exosome_genes.txt")); eg <- eg[nzchar(eg)]
pool <- intersect(eg, rownames(M)); say("外泌体基因集中可检出: %d 个", length(pool))
sam <- intersect(colnames(M), rownames(rk)); say("有生存信息的肿瘤样本: %d 例（事件 %d）", length(sam), sum(rk[sam, "fustat"]))
Mp <- M[pool, sam, drop = FALSE]
Mp <- Mp[rowMeans(Mp) > 0.5, , drop = FALSE]; say("经平均表达过滤后候选特征: %d 个", nrow(Mp))
X <- t(Mp); d0 <- data.frame(rk[sam, c("futime", "fustat")], X, check.names = FALSE)
d0 <- na.omit(d0); say("建模集: %d 例", nrow(d0))

## ---------- (A) 无差异筛选的建模 ----------
say("\n===== (A) 跳过肿瘤-正常差异筛选，直接在 TCGA 生存数据建模 =====")
pv <- apply(d0[, -(1:2)], 2, function(v) tryCatch(summary(coxph(Surv(d0$futime, d0$fustat) ~ v))$coefficients[1, 5],
                                                 error = function(e) NA_real_))
keep <- names(pv)[!is.na(pv) & pv < 0.05]; say("单因素 P<0.05: %d 个", length(keep))
set.seed(20260919)
lf <- cv.glmnet(as.matrix(d0[, keep, drop = FALSE]), Surv(d0$futime, d0$fustat), family = "cox", alpha = 1, nfolds = 10)
cf <- as.matrix(coef(lf, s = "lambda.min")); sel <- rownames(cf)[cf[, 1] != 0]
say("LASSO(lambda.min) 保留: %d 个", length(sel))
m <- step(coxph(Surv(futime, fustat) ~ ., data = d0[, c("futime", "fustat", sel)]), direction = "both", trace = 0)
mg <- names(coef(m)); say("逐步回归最终模型: %d 个基因: %s", length(mg), paste(mg, collapse = ", "))
lp_new <- predict(m, type = "lp")
c_new <- summary(coxph(Surv(d0$futime, d0$fustat) ~ lp_new))$concordance[1]
pub <- fread("04_results/20.SHAP/multiCox.txt")$id
lp_pub <- as.numeric(as.matrix(d0[, pub, drop = FALSE]) %*% fread("04_results/20.SHAP/multiCox.txt")$coef)
c_pub <- summary(coxph(Surv(d0$futime, d0$fustat) ~ lp_pub))$concordance[1]
say("新模型 C 指数 = %.3f；已发表 15 基因评分 C 指数 = %.3f", c_new, c_pub)
say("基因重叠: %d/%d（新模型 ∩ 已发表 15 基因）= %s", length(intersect(mg, pub)), length(mg),
    paste(intersect(mg, pub), collapse = ", "))
say("两个评分的 Spearman 相关 = %.3f", cor(lp_new, lp_pub, method = "spearman"))

## ---------- (B) 13 基因缩减模型在 TCGA 内的表现 ----------
say("\n===== (B) 13 基因缩减模型（缺 ADGRE5/AQP1）在 TCGA 内的表现 =====")
drop13 <- c("ADGRE5", "AQP1")                      # CGGA 队列中缺失的两个基因
d13 <- setdiff(pub, drop13); say("缩减为 %d 个基因（剔除: %s）", length(d13), paste(drop13, collapse = ", "))
co <- as.numeric(coef(coxph(Surv(d0$futime, d0$fustat) ~ ., data = d0[, c("futime", "fustat", d13)])))
lp13 <- as.numeric(as.matrix(d0[, d13, drop = FALSE]) %*% co)
c13 <- summary(coxph(Surv(d0$futime, d0$fustat) ~ lp13))$concordance[1]
say("13 基因模型在 TCGA 重拟合的 C 指数 = %.3f（15 基因模型 = %.3f）", c13, c_pub)
say("13 基因评分与 15 基因评分相关 = %.3f", cor(lp13, lp_pub, method = "spearman"))

writeLines(c(sprintf("== 敏感性分析关键输出（%s）==", format(Sys.time(), "%Y-%m-%d %H:%M")), LOG),
           file.path(OUTD, "screen_free_model.txt"))
writeLines(LOG, file.path(OUTD, "screen_free_model_log.txt"))
cat("written:", file.path(OUTD, "screen_free_model.txt"), "\n")