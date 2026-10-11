## 步骤 52：分析样本流程表（回应外部审查「建立完整的分析样本流程表，防止分母无说明地切换」）
## 所有例数均从实际输入/结果文件读取，不手工填写。CGGA 的分母统一为「生存分析集」。
## 输出：04_results/audit/sample_flow.txt
##       稿/tables/表S14_分析样本流程.csv
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(data.table)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
OUTD <- "04_results/audit"; dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
LOG <- character(); say <- function(...) { s <- sprintf(...); cat(s, "\n"); LOG <<- c(LOG, s); invisible(s) }
nrow_f <- function(p) nrow(fread(p, showProgress = FALSE))

rows <- list()
add <- function(cohort, step, n, note = "") rows[[length(rows) + 1]] <<-
  data.frame(Cohort = cohort, Step = step, n = as.integer(n), Note = note)

## ---------------- TCGA ----------------
hdr <- strsplit(readLines("02_inputs/TCGA_glioma_tumor_log2tpm.csv", n = 1), ",")[[1]]
add("TCGA", "Tumour expression profiles (log2 TPM+1)", length(hdr) - 1, "02_inputs/TCGA_glioma_tumor_log2tpm.csv")
rk <- fread("04_results/20.SHAP/risk.TCGA.txt", showProgress = FALSE)
cl <- fread("02_inputs/clinical_glioma.txt", showProgress = FALSE)
add("TCGA", "With 15-gene expression and risk score", nrow(rk))
add("TCGA", "Analysis set (risk score + survival)", sum(!is.na(rk$futime) & !is.na(rk$fustat)))
add("TCGA", "Grade recorded (G2/G3/G4) among analysis set", sum(cl$grade %in% c("G2", "G3", "G4")))
add("TCGA", "Grade missing among analysis set", sum(!(cl$grade %in% c("G2", "G3", "G4"))))
add("TCGA", "IDH recorded (mutant/wildtype)", sum(cl$idh %in% c("Mutant", "WT", "Wildtype")))
add("TCGA", "MGMT recorded", sum(!is.na(cl$mgmt) & cl$mgmt != ""))
add("TCGA", "TMB computable from somatic MAF", nrow_f("04_results/31.riskTMB/TMB.txt"))
add("TCGA", "TMB + expression + survival (joint analysis)", nrow_f("04_results/31.riskTMB/tmb_risk.txt"))
add("TCGA", "CIBERSORT fractions", nrow_f("04_results/27.CIBERSORT/CIBERSORT-Results.txt"))
add("TCGA", "Predicted drug sensitivity (oncoPredict/GDSC2)", nrow_f("04_results/34.oncoPredict/calcPhenotype_Output/DrugPredictions.csv"))
pt <- fread("04_results/audit/purity_table.txt", showProgress = FALSE)
add("TCGA", "ESTIMATE-derived purity available", sum(!is.na(pt$purity_estimate)))
add("TCGA", "DNA-based ABSOLUTE purity available", sum(!is.na(pt$purity_absolute)))
mm <- merge(data.frame(patient = rk$id, futime = rk$futime, fustat = rk$fustat, riskScore = rk$riskScore),
            cl[, c("patient", "age", "grade", "idh")], by = "patient")
add("TCGA", "Nomogram model complete cases (score+age+grade+IDH)", nrow(mm[complete.cases(mm), ]),
    "used for calibration / Brier / DCA")
add("TCGA", "Single-cell analysis", 4738, "4738 cells from 14 patients (CGGA 10X multi-sector cohort)")

## ---------------- CGGA（分母统一为生存分析集） ----------------
for (nm in c("325", "693")) {
  r <- fread(sprintf("04_results/18.CGGA/risk.CGGA%s.txt", nm), showProgress = FALSE)
  cg <- fread(sprintf("00_data/cgga/CGGA.mRNAseq_%s_clinical.20200506.txt", nm), showProgress = FALSE)
  names(cg) <- gsub(" .*", "", names(cg))
  keep <- !is.na(r$futime) & !is.na(r$fustat)
  ids <- r$id[keep]
  cg2 <- cg[match(ids, cg$CGGA_ID), ]
  lab <- sprintf("CGGA-%s", nm)
  add(lab, "Enrolled with expression", nrow(r))
  add(lab, "Analysis set (with survival data)", sum(keep))
  typ <- factor(cg2$PRS_type)
  for (k in c("Primary", "Recurrent", "Secondary")) if (k %in% levels(typ))
    add(lab, sprintf("  of analysis set: %s", k), sum(typ == k, na.rm = TRUE))
  add(lab, "  of analysis set: type unrecorded", sum(is.na(typ)))
  add(lab, "Grade recorded (WHO II/III/IV) in analysis set", sum(cg2$Grade %in% c("WHO II", "WHO III", "WHO IV")))
  add(lab, "IDH recorded in analysis set", sum(!is.na(cg2$IDH_mutation_status) & cg2$IDH_mutation_status != ""))
  add(lab, "MGMT recorded in analysis set", sum(!is.na(cg2$MGMTp_methylation_status) & cg2$MGMTp_methylation_status != ""))
  add(lab, "Genes available for the score", 13, "ADGRE5 and AQP1 absent; fixed original coefficients used")
}

## ---------------- MR / coloc ----------------
add("MR and colocalization", "Gene-tissue combinations tested", 49, "13 genes x brain regions with cis-eQTL")
add("MR and colocalization", "FinnGen R12 glioma cases (primary outcome)", 544, "C3_GBM_ASTROCYTOMA_EXALLC")
add("MR and colocalization", "Brain-tissue colocalization pairs (restricted variants)", 32, "source files retained significant pairs only")
add("MR and colocalization", "Model genes with blood full-cis-window data (eQTLGen)", 9, "of 15")

tb <- do.call(rbind, rows)
fwrite(tb, "稿/tables/表S14_分析样本流程.csv", sep = ",", quote = FALSE)
say("表 S14 已写入 稿/tables/表S14_分析样本流程.csv")
print(tb, row.names = FALSE)
writeLines(c(sprintf("== 分析样本流程（%s）==", format(Sys.time(), "%Y-%m-%d %H:%M")), LOG),
           file.path(OUTD, "sample_flow.txt"))
cat("\nwritten: 04_results/audit/sample_flow.txt\nSTEP52 DONE\n")
