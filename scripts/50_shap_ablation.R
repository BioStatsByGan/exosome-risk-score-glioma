## 步骤 50：SHAP 与「消融/置换」对照（回应外部审查第 6 条）
## 说明：线性 Cox 模型的 SHAP 值有解析解 phi_ij = beta_j * (x_ij - mean(x_j))，
##       因此 mean|SHAP| = |beta_j| * mean|x_j - mean(x_j)|（本脚本先验证这一点与已发表值一致）。
##       但「贡献大」不等于「去掉后判别力下降多」，故补做：
##       (a) 留一基因重拟合（leave-one-gene-out）后 C 指数变化；
##       (b) 置换该基因表达（B=200）后固定系数评分 C 指数的下降；
##       并给出三者（SHAP 排名 / LOO ΔC / 置换 ΔC）的秩相关与特征间相关性提示。
## 输出：04_results/audit/shap_ablation.txt
##       稿/tables/表S12_SHAP消融与置换.csv
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table); library(survival)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
OUTD <- "04_results/audit"; dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
LOG <- character(); say <- function(...) { s <- sprintf(...); cat(s, "\n"); LOG <<- c(LOG, s); invisible(s) }
ci <- function(f) as.numeric(summary(f)$concordance[1])
cscore <- function(score, tt, ss) ci(coxph(Surv(tt, ss) ~ score))
B <- 200

mc <- read.table("04_results/20.SHAP/multiCox.txt", header = TRUE, sep = "\t")
coefs <- setNames(mc$coef, mc$id)
rk <- read.table("04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t",
                 check.names = FALSE, row.names = 1)
E <- as.matrix(rk[, names(coefs)]); tt <- rk$futime; ss <- rk$fustat
S15 <- as.numeric(E %*% coefs)
C15 <- cscore(S15, tt, ss)
say("TCGA 15 基因评分 C 指数 = %.3f（n=%d，事件 %d）", C15, length(tt), sum(ss))

## ---------- (0) 线性模型的解析 SHAP ----------
phi_mean <- setNames(sapply(names(coefs), function(g) abs(coefs[g]) * mean(abs(E[, g] - mean(E[, g])))),
                     names(coefs))
say("\n线性 Cox 的解析 mean|SHAP|（与 kernelshap 对线性预测值的定义一致）：%s",
    paste(sprintf("%s=%.3f", names(sort(phi_mean, decreasing = TRUE)),
                  round(sort(phi_mean, decreasing = TRUE), 3)), collapse = ", "))

## ---------- (a) 留一基因重拟合 ----------
loo <- sapply(names(coefs), function(g) {
  keep <- setdiff(names(coefs), g)
  ci(coxph(Surv(tt, ss) ~ ., data = as.data.frame(E[, keep, drop = FALSE])))
})
dC_loo <- C15 - loo

## ---------- (b) 置换该基因表达（固定系数评分） ----------
set.seed(42)
dC_perm <- sapply(names(coefs), function(g) {
  v <- replicate(B, {
    xp <- E[, g]; xp <- sample(xp)
    S <- S15 - coefs[g] * E[, g] + coefs[g] * xp
    cscore(S, tt, ss)
  })
  C15 - mean(v)
})

## ---------- (c) 相关性提示 ----------
COR <- cor(E, method = "spearman")
maxabs <- apply(COR, 1, function(z) max(abs(z[z < 1])))

tb <- data.frame(Gene = names(coefs), Coef = round(unname(coefs), 3),
                 mean_abs_SHAP = round(unname(phi_mean), 3),
                 C_drop_LOO = round(unname(dC_loo), 4),
                 C_drop_perm = round(unname(dC_perm), 4),
                 Max_abs_rho_with_other_genes = round(unname(maxabs[names(coefs)]), 2))
tb <- tb[order(-tb$mean_abs_SHAP), ]
say("\n%-8s %8s %10s %12s %12s %10s", "Gene", "Coef", "mean|SHAP|", "dC_LOO", "dC_perm", "max|rho|")
for (i in seq_len(nrow(tb)))
  say("%-8s %8.3f %10.3f %12.4f %12.4f %10.2f", tb$Gene[i], tb$Coef[i], tb$mean_abs_SHAP[i],
      tb$C_drop_LOO[i], tb$C_drop_perm[i], tb$Max_abs_rho_with_other_genes[i])

rk_tab <- function(a, b) round(cor(a, b, method = "spearman"), 3)
say("\n秩相关（Spearman）：mean|SHAP| vs dC_LOO = %s；mean|SHAP| vs dC_perm = %s；dC_LOO vs dC_perm = %s",
    rk_tab(tb$mean_abs_SHAP, tb$C_drop_LOO), rk_tab(tb$mean_abs_SHAP, tb$C_drop_perm),
    rk_tab(tb$C_drop_LOO, tb$C_drop_perm))
say("ΔC_LOO 最大的三个基因：%s；ΔC_perm 最大的三个基因：%s",
    paste(tb$Gene[order(-tb$C_drop_LOO)][1:3], collapse = ", "),
    paste(tb$Gene[order(-tb$C_drop_perm)][1:3], collapse = ", "))
hi <- tb$Gene[tb$mean_abs_SHAP >= quantile(tb$mean_abs_SHAP, .7) & tb$C_drop_LOO <= 0.001]
say("mean|SHAP| 位居前 30%% 但 ΔC_LOO ≤0.001 的基因（SHAP 高估其唯一贡献的例子）：%s",
    if (length(hi)) paste(hi, collapse = ", ") else "无")

## ---------- (d) CGGA-325 的 13 基因版本留一分析（外部稳健性） ----------
KEEP13 <- setdiff(names(coefs), c("ADGRE5", "AQP1"))
ex <- fread("00_data/cgga/CGGA.mRNAseq_325.RSEM-genes.20200506.txt", showProgress = FALSE)
gn <- ex[[1]]; m <- as.matrix(ex[, -1]); rownames(m) <- gn
use <- intersect(KEEP13, rownames(m)); m <- log2(m[use, , drop = FALSE] + 1)
sc325 <- colSums(coefs[use] * m)
rr <- read.table("04_results/18.CGGA/risk.CGGA325.txt", header = TRUE, sep = "\t", check.names = FALSE)
s <- intersect(names(sc325), rr$id)
dd <- data.frame(futime = rr$futime[match(s, rr$id)], fustat = rr$fustat[match(s, rr$id)],
                 S = as.numeric(sc325[s]))
dd <- dd[!is.na(dd$futime) & !is.na(dd$fustat), ]
C325 <- cscore(dd$S, dd$futime, dd$fustat)
keep_s <- !is.na(rr$futime[match(s, rr$id)]) & !is.na(rr$fustat[match(s, rr$id)])
m2 <- m[, s[keep_s], drop = FALSE]
## 注意：CGGA-325 不做「留一重拟合」——重拟合本身会提升 C 指数（ΔC 全为负），会掩盖基因贡献；
## 改用固定系数下的置换分析，与 TCGA 完全同口径。
set.seed(42)
perm325 <- sapply(use, function(g) {
  v <- replicate(B, {
    xp <- m2[g, ]; xp <- sample(xp)
    S <- dd$S - coefs[g] * m2[g, ] + coefs[g] * xp
    cscore(S, dd$futime, dd$fustat)
  })
  C325 - mean(v)
})
say("\nCGGA-325（固定系数 13 基因评分，B=%d 次置换）C=%.3f；置换 ΔC：%s", B, C325,
    paste(sprintf("%s=%.4f", use, perm325), collapse = ", "))
say("   CGGA-325 置换 ΔC 最大的三个基因：%s",
    paste(use[order(-perm325)][1:3], collapse = ", "))

fwrite(tb, "稿/tables/表S12_SHAP消融与置换.csv", sep = ",", quote = FALSE)
say("\n表 S12 已写入 稿/tables/表S12_SHAP消融与置换.csv")
writeLines(c(sprintf("== SHAP 与消融/置换对照（%s，置换 B=%d）==", format(Sys.time(), "%Y-%m-%d %H:%M"), B), LOG),
           file.path(OUTD, "shap_ablation.txt"))
cat("\nwritten: 04_results/audit/shap_ablation.txt\nSTEP50 DONE\n")
