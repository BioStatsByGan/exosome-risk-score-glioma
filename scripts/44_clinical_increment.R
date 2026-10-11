## 步骤 44：外部队列的「临床模型 vs 临床模型+评分」增量分析（含 ΔC 的 bootstrap 95%CI）
## 回应外部意见第 4 条：需要在外部队列比较相同患者中的临床模型与临床模型加评分，并报告 ΔC 的置信区间。
## 输出：04_results/audit/clinical_increment.txt / _log.txt
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(survival); library(data.table)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
OUTD <- "04_results/audit"; dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
LOG <- character(); say <- function(...) { s <- sprintf(...); cat(s, "\n"); LOG <<- c(LOG, s); invisible(s) }
B <- 1000

increment <- function(d, f_cli, label) {
  d <- d[complete.cases(d[, all.vars(as.formula(f_cli))]), ]
  f_all <- as.formula(paste(f_cli, "+ riskScore"))
  fit <- function(dat, form) coxph(as.formula(form), data = dat)
  c_cli <- summary(fit(d, f_cli))$concordance[1]
  c_all <- summary(fit(d, f_all))$concordance[1]
  set.seed(7)
  deltas <- numeric(B); ok <- 0
  for (b in seq_len(B)) {
    idx <- sample(nrow(d), replace = TRUE)
    db <- d[idx, ]
    if (length(unique(db$fustat)) < 2) next
    cc <- tryCatch(summary(fit(db, f_cli))$concordance[1], error = function(e) NA)
    ca <- tryCatch(summary(fit(db, f_all))$concordance[1], error = function(e) NA)
    ok <- ok + 1; deltas[ok] <- ca - cc
  }
  deltas <- deltas[1:ok]
  ci <- quantile(deltas, c(0.025, 0.975), na.rm = TRUE)
  say("%s（n=%d，事件 %d）：临床模型 C=%.3f → 加评分 C=%.3f，ΔC=%.3f（bootstrap 95%%CI %.3f 至 %.3f，B=%d）",
      label, nrow(d), sum(d$fustat), c_cli, c_all, c_all - c_cli, ci[1], ci[2], ok)
  data.frame(Cohort = label, n = nrow(d), Events = sum(d$fustat),
             C_clinical = c_cli, C_plus_score = c_all, deltaC = c_all - c_cli,
             CI_low = ci[1], CI_high = ci[2], B = ok)
}

## ---------- TCGA ----------
rk <- read.table("04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
cl <- read.table("02_inputs/clinical_glioma.txt", header = TRUE, sep = "\t", check.names = FALSE)
dt <- merge(data.frame(id = rownames(rk), futime = rk$futime, fustat = rk$fustat, riskScore = rk$riskScore),
            cl[, c("patient", "age", "sex", "grade", "idh")], by.x = "id", by.y = "patient")
dt$Grade <- factor(dt$grade); dt$IDH <- factor(dt$idh); dt$Sex <- factor(dt$sex)
r1 <- increment(dt, "Surv(futime, fustat) ~ age + Sex + Grade + IDH", "TCGA（年龄+性别+分级+IDH）")

## ---------- CGGA-325 / 693 ----------
out <- r1
for (nm in c("325", "693")) {
  cg <- read.table(sprintf("04_results/18.CGGA/risk.CGGA%s.txt", nm), header = TRUE, sep = "\t", check.names = FALSE)
  cgc <- fread(sprintf("00_data/cgga/CGGA.mRNAseq_%s_clinical.20200506.txt", nm), data.table = FALSE)
  names(cgc) <- gsub(" .*", "", names(cgc))
  d <- merge(cg[, c("id", "futime", "fustat", "riskScore")],
             cgc[, c("CGGA_ID", "Age", "Gender", "Grade", "IDH_mutation_status", "PRS_type")],
             by.x = "id", by.y = "CGGA_ID")
  d <- d[!is.na(d$futime) & !is.na(d$fustat), ]
  d$Grade <- factor(d$Grade, levels = c("WHO II", "WHO III", "WHO IV"))
  d$IDH <- factor(ifelse(grepl("^mut", d$IDH_mutation_status, ignore.case = TRUE), "Mutant", "Wildtype"),
                  levels = c("Mutant", "Wildtype"))
  d$Age <- as.numeric(d$Age); d$Sex <- factor(d$Gender)
  say("\nCGGA-%s：可分析 %d 例；PRS_type 构成：%s", nm, nrow(d),
      paste(sprintf("%s=%d", names(table(d$PRS_type)), as.integer(table(d$PRS_type))), collapse = " "))
  out <- rbind(out, increment(d, "Surv(futime, fustat) ~ Age + Sex + Grade", sprintf("CGGA-%s（年龄+性别+分级）", nm)))
  out <- rbind(out, increment(d, "Surv(futime, fustat) ~ Age + Sex + Grade + IDH", sprintf("CGGA-%s（年龄+性别+分级+IDH）", nm)))
}
fwrite(out, file.path(OUTD, "clinical_increment.txt"), sep = "\t", quote = FALSE)
writeLines(c(sprintf("== 临床模型 vs 加评分：ΔC 与 bootstrap 95%%CI（%s）==", format(Sys.time(), "%Y-%m-%d %H:%M")), LOG),
           file.path(OUTD, "clinical_increment_log.txt"))
cat("\n"); print(out, digits = 3, row.names = FALSE)