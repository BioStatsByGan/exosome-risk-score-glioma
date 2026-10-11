suppressMessages(library(data.table)); setwd("F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤")
cl <- fread("02_inputs/clinical_glioma.txt"); rk <- fread("04_results/20.SHAP/risk.TCGA.txt")
tc <- merge(rk[, .(id, futime, fustat, riskScore, Risk)],
            cl[, .(patient, age, sex, grade, idh, mgmt)], by.x = "id", by.y = "patient")
tc <- tc[!(is.na(futime) | is.na(fustat))]
tc[, Grade := paste0("WHO grade ", sub("^G", "", grade))]
tc[, IDH := ifelse(grepl("mut", idh, ignore.case = TRUE), "IDH-mutant", "IDH-wildtype")]
cgcl <- fread("00_data/cgga/CGGA.mRNAseq_325_clinical.20200506.txt")
setnames(cgcl, gsub(" .*", "", names(cgcl)))
cg <- merge(fread("04_results/18.CGGA/risk.CGGA325.txt")[, .(id, futime, fustat, riskScore, Risk, Grade, IDH)],
            cgcl[, .(CGGA_ID, Radio_status, `1p19q_codeletion_status`)], by.x = "id", by.y = "CGGA_ID")
cg <- cg[!(is.na(futime) | is.na(fustat))]
cg[, Codel := ifelse(grepl("^codel", `1p19q_codeletion_status`, ignore.case = TRUE), "1p/19q codeleted", "1p/19q intact")]
cg[, Radio := ifelse(Radio_status == 1, "Radiotherapy treated", "Radiotherapy untreated")]
cat("tc IDH-wt:", nrow(tc[IDH == "IDH-wildtype"]), "| tc IDH-mut:", nrow(tc[IDH == "IDH-mutant"]),
    "| tc G3:", nrow(tc[Grade == "WHO grade III"]), "\n")
cat("cg IDH-wt:", nrow(cg[IDH == "IDH-wildtype"]), "| cg IDH-mut:", nrow(cg[IDH == "IDH-mutant"]),
    "| cg G4:", nrow(cg[Grade == "WHO IV"]), "| cg codel-intact:", nrow(cg[Codel == "1p/19q intact"]),
    "| cg radio-tx:", nrow(cg[Radio == "Radiotherapy treated"]), "\n")
cat("cg Grade:", paste(unique(cg$Grade), collapse = "/"), "\n")
cat("cg IDH:", paste(unique(cg$IDH), collapse = "/"), "\n")
cat("cg 1p19q:", paste(unique(cg$`1p19q_codeletion_status`), collapse = "/"), "\n")
