suppressMessages(library(data.table))
setwd("F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤")
cg <- fread("04_results/18.CGGA/risk.CGGA325.txt")
cat("risk:", nrow(cg), "| NA fustat:", sum(is.na(cg$fustat)), "| NA futime:", sum(is.na(cg$futime)), "\n")
cgcl <- fread("00_data/cgga/CGGA.mRNAseq_325_clinical.20200506.txt")
setnames(cgcl, gsub(" .*", "", names(cgcl)))
cat("临床列:", paste(names(cgcl), collapse=" | "), "\n")
cat("临床行数:", nrow(cgcl), "\n")
m <- merge(cg[, .(id, futime, fustat, riskScore, Risk, Grade, IDH)],
           cgcl[, .(CGGA_ID, Grade_cl = Grade, Gender, Age, Radio_status, Chemo_status,
                    IDH_mutation_status, `1p19q_codeletion_status`, MGMTp_methylation_status)],
           by.x = "id", by.y = "CGGA_ID")
cat("merged:", nrow(m), "| NA fustat:", sum(is.na(m$fustat)), "\n")
print(head(m, 3))
