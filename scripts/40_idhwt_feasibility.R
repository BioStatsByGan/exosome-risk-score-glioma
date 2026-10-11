## 步骤 40：可行性验证 —— IDH 野生型专用签名能否超过现有评分
## 目的：若在 IDH-wt 内部重建的签名能显著超过现有 15 基因评分（TCGA 内 C 0.608），
##       则主线可从"外泌体签名"改为"IDH 野生型专用分层"，把增量不足的弱点变成研究动机。
## 评估：10 折 × 5 次重复交叉验证的 Harrell C（折外预测汇总），对照
##       1) 现有 15 基因评分  2) 临床模型（年龄+性别+MGMT+分级）  3) 新签名
## 输出：04_results/40.idhwt/idhwt_feasibility.txt
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(survival); library(glmnet)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
OUTD <- "04_results/40.idhwt"; dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
say <- function(...) cat(sprintf(...), "\n")

## ---- 数据 ----
ex <- fread("02_inputs/TCGA_glioma_tumor_log2tpm.csv")
M <- as.matrix(ex[, -1]); rownames(M) <- ex[[1]]; M <- M[!duplicated(rownames(M)), ]
ck <- fread("02_inputs/clinical_glioma.txt")
rk <- fread("04_results/20.SHAP/risk.TCGA.txt")
d <- merge(rk[, .(id, futime, fustat, riskScore)], ck[, .(patient, age, sex, grade, idh, mgmt)],
           by.x = "id", by.y = "patient")
d <- d[!is.na(futime) & !is.na(fustat)]
d <- d[idh == "WT"]
say("IDH 野生型样本 %d，死亡 %d", nrow(d), sum(d$fustat))

X_all <- t(M[, d$id])
ex_genes <- intersect(fread("02_inputs/exosome_genes.txt")[[1]], colnames(X_all))
X <- X_all[, ex_genes]
say("候选特征（外泌体相关基因，可检出）: %d", ncol(X))

## ---- 1) 现有 15 基因评分的 C 指数 ----
c_exist <- summary(coxph(Surv(futime, fustat) ~ riskScore, data = d))$concordance[1]   # 与临床模型同口径（风险评分越大风险越高）
say("现有 15 基因评分 C 指数: %.3f", c_exist)

## ---- 2) 临床模型 C 指数 ----
cli <- data.frame(futime = d$futime, fustat = d$fustat, age = d$age,
                  sex = factor(d$sex), grade = factor(d$grade), mgmt = factor(d$mgmt))
fit_cli <- coxph(Surv(futime, fustat) ~ ., data = cli)
c_cli <- summary(fit_cli)$concordance[1]
say("临床模型（年龄+性别+分级+MGMT）C 指数: %.3f", c_cli)

## ---- 3) 折外交叉验证的 LASSO 新签名 ----
cv_cindex <- function(X, time, status, nfold = 10, nrep = 5, seed = 20260926) {
  cs <- c()
  for (r in seq_len(nrep)) {
    set.seed(seed + r)
    preds <- rep(NA_real_, length(time))
    fold <- sample(rep(1:nfold, length.out = length(time)))
    for (k in 1:nfold) {
      tr <- fold != k; te <- fold == k
      cvfit <- cv.glmnet(X[tr, ], Surv(time[tr], status[tr]), family = "cox", alpha = 1,
                         nfolds = 5, standardize = TRUE)
      preds[te] <- as.numeric(predict(cvfit, newx = X[te, ], s = "lambda.min", type = "link"))
    }
    cs <- c(cs, survival::concordance(Surv(time, status) ~ preds)$concordance)
  }
  c(mean = mean(cs), sd = sd(cs))
}
res <- cv_cindex(X, d$futime, d$fustat)
say("新签名（外泌体基因，折外 CV）C 指数: %.3f (SD %.3f)", res["mean"], res["sd"])

## ---- 4) 全基因组特征对照 ----
Xg <- X_all[, setdiff(colnames(X_all), ex_genes)]
sel <- apply(Xg, 2, function(v) sd(v) > 0.1)
Xg <- Xg[, sel]
say("全基因组特征数: %d", ncol(Xg))
resg <- tryCatch(cv_cindex(Xg, d$futime, d$fustat, nrep = 3), error = function(e) c(mean = NA, sd = NA))
say("全基因组签名 C 指数: %.3f", resg["mean"])

out <- c(sprintf("IDH-wt n=%d events=%d", nrow(d), sum(d$fustat)),
         sprintf("现有 15 基因评分      C = %.3f", c_exist),
         sprintf("临床模型              C = %.3f", c_cli),
         sprintf("新签名(外泌体基因)    C = %.3f (SD %.3f, 5x10 折)", res["mean"], res["sd"]),
         sprintf("新签名(全基因组)      C = %.3f", resg["mean"]))
writeLines(out, file.path(OUTD, "idhwt_feasibility.txt"))
say("\n== 结论 =="); for (o in out) say("%s", o)
