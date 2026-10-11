## 步骤 30：TCGA 免疫分型（Thorsson 六型）与风险分层的关联
## 输入：00_data/TCGA_immune_MAE.rda、04_results/20.SHAP/risk.TCGA.txt
## 输出：04_results/30.immSubtype/immune_subtype_risk.txt、03_figures/Fig8A_immSubtype.pdf
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({library(MultiAssayExperiment); library(SummarizedExperiment); library(ggplot2); library(data.table)})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT); source("scripts/plot_style.R")
dir.create("04_results/30.immSubtype", recursive = TRUE, showWarnings = FALSE)

e <- new.env(); load("00_data/TCGA_immune_MAE.rda", envir = e)
obj <- get(ls(e)[1], envir = e)
cd <- as.data.frame(SummarizedExperiment::colData(obj))
cat("colData 列:", paste(colnames(cd), collapse = " | "), "\n")
cand <- grep("Subtype|subtype|Immune|immune", colnames(cd), value = TRUE)
cat("候选列:", paste(cand, collapse = ", "), "\n")
col <- cand[grep("Subtype", cand)][1]
if (is.na(col)) col <- cand[1]
sub <- data.frame(patient = substr(rownames(cd), 1, 12), subtype = cd[[col]])
cat("使用列:", col, "| 样本:", nrow(sub), "| 分型分布:\n"); print(table(sub$subtype, useNA = "ifany"))

risk <- read.table("04_results/20.SHAP/risk.TCGA.txt", header = TRUE, sep = "\t",
                   check.names = FALSE, row.names = 1)
m <- merge(sub, data.frame(patient = rownames(risk), Risk = risk$Risk), by = "patient")
m <- m[!is.na(m$subtype) & m$subtype != "", ]
cat("匹配样本:", nrow(m), "\n")
tb <- table(m$subtype, m$Risk)
print(tb)
tb <- tb[rowSums(tb) > 0, , drop = FALSE]
pv <- tryCatch(chisq.test(tb, simulate.p.value = TRUE, B = 5000)$p.value, error = function(x) NA)
cat("卡方 p =", signif(pv, 3), "\n")
fwrite(m, "04_results/30.immSubtype/immune_subtype_risk.txt", sep = "\t", quote = FALSE)

## 堆叠比例图
m$Risk <- factor(m$Risk, levels = c("low", "high"), labels = c("Low risk", "High risk"))
pd <- as.data.frame(prop.table(table(m$subtype, m$Risk), 2) * 100)
colnames(pd) <- c("Subtype", "Risk", "Pct")
pd <- pd[!is.na(pd$Subtype) & pd$Subtype != "", ]
fills <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#7F7F7F")
p <- ggplot(pd, aes(Risk, Pct, fill = Subtype)) +
  geom_col(width = 0.62, colour = "grey25", linewidth = 0.3) +
  geom_text(aes(label = ifelse(Pct >= 5, sprintf("%.0f", Pct), "")),
            position = position_stack(vjust = 0.5), size = 2.8) +
  scale_fill_manual(values = fills, name = "") +
  annotate("text", x = 0.62, y = 102, hjust = 0, size = 3,
           label = ifelse(is.na(pv), "", sprintf("Chi-square P = %s", ifelse(pv < 0.001, "<0.001", sprintf("%.3f", pv))))) +
  labs(x = NULL, y = "Proportion of patients (%)") +
  ylim(0, 108) +
  theme(legend.position = "right")
save_pdf(p, "03_figures/Fig8A_immSubtype.pdf", 6.6, 4.6)
cat("STEP30 DONE\n")
