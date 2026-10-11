## 胶质瘤项目 输入准备：TCGA-GBM/LGG（本地成表）+ GTEx 脑（TOIL）→ 课程流水线输入格式
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table)})
base <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
dat  <- file.path(base, "00_data"); inp <- file.path(base, "02_inputs")
dir.create(inp, recursive = TRUE, showWarnings = FALSE)
setwd(base)

rd <- function(f) fread(f, data.table = FALSE, check.names = FALSE)

## ---------- 1 本地 TCGA-GBM / LGG 表达（行=基因符号，列=barcode） ----------
merge_tcga <- function(csv) {
  m <- rd(file.path(dat, csv))
  rn <- m[[1]]; m <- as.matrix(m[, -1]); rownames(m) <- rn
  m <- m[!duplicated(rownames(m)), , drop = FALSE]
  m
}
gbm <- merge_tcga("TCGA-GBM_mrna_expr_tpm.csv")
lgg <- merge_tcga("TCGA-LGG_mrna_expr_tpm.csv")
cat("GBM:", dim(gbm), " LGG:", dim(lgg), "\n")
gene <- intersect(rownames(gbm), rownames(lgg))
tcga <- cbind(gbm[gene, , drop = FALSE], lgg[gene, , drop = FALSE])
cat("合并后:", dim(tcga), "\n")

## ---------- 2 样本类型与患者级折叠 ----------
p4 <- sapply(strsplit(colnames(tcga), "-"), `[`, 4)
type <- substr(p4, 1, 2)
pat  <- sapply(strsplit(colnames(tcga), "-"), function(x) paste(x[1:3], collapse = "-"))
isT <- type %in% c("01", "02", "03", "06")     # 原发/复发肿瘤
isN <- type == "11"                            # 癌旁正常
cat("肿瘤样本:", sum(isT), " 癌旁正常:", sum(isN), "\n")

collapse_patient <- function(m, keep) {
  m <- m[, keep, drop = FALSE]
  ids <- sapply(strsplit(colnames(m), "-"), function(x) paste(x[1:3], collapse = "-"))
  agg <- rowsum(t(m), group = ids)                 # 样本×基因 按患者求和
  cnt <- as.vector(table(ids)[rownames(agg)])
  t(agg / cnt)
}
tum <- collapse_patient(tcga, isT)
cat("患者级肿瘤矩阵:", dim(tum), "\n")

## ---------- 3 GTEx 脑（TOIL，ENSG→symbol） ----------
phe <- rd(file.path(dat, "TcgaTargetGTEX_phenotype.txt.gz"))
brain <- phe[phe[["_primary_site"]] == "Brain" & phe[["_study"]] == "GTEX", "sample"]
cat("GTEx 脑样本:", length(brain), "\n")
hdr <- fread(file.path(dat, "TcgaTargetGtex_rsem_gene_tpm.gz"), nrows = 0)
cols <- colnames(hdr)
idx <- match(brain, cols); idx <- idx[!is.na(idx)]
sel <- c(1L, idx)
toil <- fread(file.path(dat, "TcgaTargetGtex_rsem_gene_tpm.gz"), select = sel, data.table = FALSE,
              check.names = FALSE, showProgress = FALSE)
cat("TOIL 子集:", dim(toil), "\n")
gi <- rd(file.path(dat, "TCGA-GBM_gene_info.csv"))
map <- setNames(gi$gene_name, sub("\\..*$", "", gi$gene_id))     # 去掉 ENSG 版本号再匹配
sym <- map[sub("\\..*$", "", toil[[1]])]
keep <- !is.na(sym) & sym != "" & !duplicated(sym)
gx <- as.matrix(toil[keep, -1]); rownames(gx) <- sym[keep]
gx <- gx[intersect(rownames(gx), rownames(tum)), , drop = FALSE]
cat("GTEx 脑 × 共有基因:", dim(gx), "\n")

## ---------- 4 生存与临床（clinicalSE，患者级） ----------
mk_cli <- function(cs) {
  d <- rd(file.path(dat, cs))
  d$patient <- sub("^(TCGA-[^-]+-[^-]+).*$", "\\1", d$barcode)
  num <- function(x) suppressWarnings(as.numeric(x))
  d$os_days <- ifelse(!is.na(num(d$days_to_death)), num(d$days_to_death), num(d$days_to_last_follow_up))
  d$os_stat <- ifelse(d$vital_status == "Dead", 1, 0)
  d$age <- num(d$age_at_diagnosis) / 365.25
  d$sex <- ifelse(grepl("female", d$gender, ignore.case = TRUE), "Female",
                  ifelse(grepl("male", d$gender, ignore.case = TRUE), "Male", NA))
  d$grade <- d[["paper_Grade"]]
  d$idh <- d[["paper_IDH.status"]]
  d$mgmt <- d[["paper_MGMT.promoter.status"]]
  d$cohort <- ifelse(grepl("^TCGA-GBM", cs), "GBM", "LGG")
  keep <- c("patient", "cohort", "os_days", "os_stat", "age", "sex", "grade", "idh", "mgmt")
  d <- d[!is.na(d$patient), keep]
  d <- d[!duplicated(d$patient), ]
  d
}
cli <- rbind(mk_cli("TCGA-GBM_clinicalSE.csv"), mk_cli("TCGA-LGG_clinicalSE.csv"))
cli <- cli[!is.na(cli$os_days) & cli$os_days > 0, ]
cli$futime <- cli$os_days / 365.25
cli$fustat <- cli$os_stat
time <- data.frame(Id = cli$patient, futime = round(cli$os_days, 1), fustat = cli$fustat)
cat("有生存的患者:", nrow(cli), " 事件数:", sum(cli$fustat, na.rm = TRUE), "\n")

## ---------- 5 输出 ----------
write.csv(round(log2(tum + 1), 4), file.path(inp, "TCGA_glioma_tumor_log2tpm.csv"))
write.csv(round(gx, 4), file.path(inp, "GTEx_brain_log2tpm.csv"))
write.table(time, file.path(inp, "time.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cli, file.path(inp, "clinical_glioma.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
file.copy(file.path(base, "../222_tcgaSHAP_胃癌外泌体预后模型/12.GeneCards/gene.txt"),
          file.path(inp, "exosome_genes.txt"), overwrite = TRUE)
cat("\n样本类型统计\n"); print(table(ifelse(isT, "Tumor", ifelse(isN, "Normal", "Other"))))
cat("\n队列统计\n"); print(table(cli$cohort, useNA = "ifany"))
cat("\nGrade:", paste(names(table(cli$grade)), table(cli$grade), sep = "=", collapse = ", "), "\n")
cat("IDH:", paste(names(table(cli$idh)), table(cli$idh), sep = "=", collapse = ", "), "\n")
cat("MGMT:", paste(names(table(cli$mgmt)), table(cli$mgmt), sep = "=", collapse = ", "), "\n")
cat("DONE ->", inp, "\n")
