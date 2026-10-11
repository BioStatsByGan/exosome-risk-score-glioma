options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(survival); library(rms); library(pec); library(data.table)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT)
risk <- read.table("04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
cli <- read.table("02_inputs/clinical_glioma.txt", header = TRUE, sep = "\t", check.names = FALSE)
d <- merge(data.frame(patient = rownames(risk), risk, check.names = FALSE), cli[, c("patient","age","sex","grade","idh")], by = "patient")
d$Grade <- factor(d$grade); d$IDH <- factor(d$idh); d$Sex <- factor(d$sex); units(d$futime) <- "Year"
dd <- datadist(d); options(datadist = "dd")
f <- cph(Surv(futime, fustat) ~ riskScore + age + Grade + IDH, data = d, x = TRUE, y = TRUE, surv = TRUE)
set.seed(2)
bs <- pec(list(`Full model` = f), data = d, formula = Surv(futime, fustat) ~ 1,
          times = c(0, 1, 3, 5), splitMethod = "BootConv", B = 50, verbose = FALSE)
cat("AppErr 结构:\n"); print(dim(bs$AppErr)); print(rownames(bs$AppErr)); print(colnames(bs$AppErr))
print(bs$AppErr)
