## 单细胞层：CGGA 10X 胶质瘤单细胞（6148 细胞）中 15 基因风险模块的定位
## 输入：00_data/cgga/CGGA.scRNAseq_6148.count.matrix.tsv、CGGA.scRNA_6148_clinical.20220104.txt
## 输出：04_results/SC/（对象与表格）、03_figures/Fig9A 至 Fig9D
options(stringsAsFactors = FALSE, warn = -1)
suppressMessages({
  library(Seurat); library(Matrix); library(data.table); library(ggplot2)
  library(patchwork); library(harmony)
})
ROOT <- "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤"; setwd(ROOT); source("scripts/plot_style.R")
dir.create("04_results/SC", recursive = TRUE, showWarnings = FALSE)
set.seed(20260920)

GENES <- c("PPIA", "UTP20", "RAB27A", "MSN", "EZH2", "GNG12", "ADGRE5", "LRRK2",
           "AQP1", "ACY3", "MYC", "EZR", "CCT3", "SNX18", "EIF2S1")

RDS <- "04_results/SC/cgga_scRNA_seurat.rds"
if (file.exists(RDS)) {
  cat("读取已建好的 Seurat 对象（跳过 QC 与聚类）\n")
  obj <- readRDS(RDS)
  cat("细胞数:", ncol(obj), " 聚类数:", length(levels(obj$seurat_clusters)), "\n")
} else {
cat("读取表达矩阵...\n")
m <- fread("00_data/cgga/CGGA.scRNAseq_6148.count.matrix.tsv", showProgress = FALSE)
gn <- m[[1]]; mat <- as.matrix(m[, -1]); rownames(mat) <- gn; rm(m)
cat("矩阵:", nrow(mat), "基因 x", ncol(mat), "细胞\n")
obj <- CreateSeuratObject(counts = mat, project = "CGGA_scRNA", min.cells = 3, min.features = 200)
rm(mat); gc()

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
cat("QC 前细胞数:", ncol(obj), "\n")
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
obj$sample <- sub("_SC.*", "", colnames(obj))
cat("QC 后细胞数:", ncol(obj), " 样本数:", length(unique(obj$sample)), "\n")

obj <- NormalizeData(obj) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData() %>%
  RunPCA(npcs = 30, verbose = FALSE) %>% RunHarmony("sample", verbose = FALSE) %>%
  RunUMAP(reduction = "harmony", dims = 1:20, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)
cat("聚类数:", length(levels(obj$seurat_clusters)), "\n")
saveRDS(obj, "04_results/SC/cgga_scRNA_seurat.rds")
}

## 细胞类型注释（谱系模块评分取最高）
MARK <- list(
  `Tumour cell` = c("GFAP", "PTPRZ1", "SOX2", "OLIG2", "NES", "VIM", "S100B"),
  Microglia = c("P2RY12", "CX3CR1", "TMEM119", "SALL1", "CSF1R"),
  Macrophage = c("CD163", "AIF1", "LYZ", "CD68", "MSR1", "MRC1"),
  `T cell` = c("CD3D", "CD3E", "CD2", "IL7R", "CD8A", "CD4"),
  `NK cell` = c("NKG7", "KLRD1", "GNLY", "PRF1"),
  `B cell` = c("MS4A1", "CD79A", "CD79B"),
  Oligodendrocyte = c("MBP", "PLP1", "MAG", "MOBP", "MOG"),
  OPC = c("PDGFRA", "CSPG4", "OLIG1"),
  Endothelial = c("CLDN5", "VWF", "PECAM1", "FLT1"),
  Neuron = c("RBFOX3", "SYT1", "SNAP25", "NEFL")
)
for (k in names(MARK)) obj <- AddModuleScore(obj, features = list(MARK[[k]]), name = gsub("[^A-Za-z]", "", k), seed = 1)
scorecols <- paste0(gsub("[^A-Za-z]", "", names(MARK)), "1")   # AddModuleScore 会加序号后缀
avg <- sapply(scorecols, function(cn) tapply(obj@meta.data[[cn]], obj$seurat_clusters, mean))
lab <- colnames(avg)[apply(avg, 1, which.max)]
names(lab) <- rownames(avg)
obj$celltype <- factor(gsub("1$", "", unname(lab[as.character(obj$seurat_clusters)])))   # 去掉 AddModuleScore 的序号后缀
cat("细胞类型分布:\n"); print(table(obj$celltype))
fwrite(data.frame(cell = colnames(obj), sample = obj$sample, cluster = obj$seurat_clusters,
                  celltype = as.character(obj$celltype), risk_module = NA),
       "04_results/SC/cell_annotation.txt", sep = "\t", quote = FALSE)

## 三种细胞级评分：
##   ① RiskModule1  15 基因无权重模块分（AddModuleScore；原文所用口径，保留作对照）
##   ② RiskWeighted 沿用 bulk Cox 系数加权（Σ coef × log-normalized 表达），是 bulk 风险评分的细胞级类比
##   ③ RiskPos / RiskNeg 正、负系数基因各成一个模块，检查方向是否被无权重平均抵消
obj <- AddModuleScore(obj, features = list(GENES), name = "RiskModule", seed = 1)
mc15 <- fread("04_results/20.SHAP/multiCox.txt")
w <- setNames(mc15$coef, mc15$id); gp <- mc15$id[mc15$coef > 0]; gn <- mc15$id[mc15$coef < 0]
obj <- AddModuleScore(obj, features = list(intersect(gp, rownames(obj))), name = "RiskPos", seed = 1)
obj <- AddModuleScore(obj, features = list(intersect(gn, rownames(obj))), name = "RiskNeg", seed = 1)
expr <- tryCatch(GetAssayData(obj, assay = "RNA", layer = "data"),
                 error = function(e) GetAssayData(obj, assay = "RNA", slot = "data"))
g15 <- intersect(names(w), rownames(expr))
obj$RiskWeighted <- as.numeric(t(w[g15]) %*% as.matrix(expr[g15, , drop = FALSE]))
## 矩阵类型核对：原始计数是否为整数（判断该数据是 UMI/count 还是已标准化）
raw <- tryCatch(GetAssayData(obj, assay = "RNA", layer = "counts"),
                error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts"))
sub_raw <- as.matrix(raw[1:min(200, nrow(raw)), 1:min(500, ncol(raw)), drop = FALSE])
cat(sprintf("矩阵类型检查：抽样内原始计数全为整数 = %s（总计数 %.3g）\n", all(sub_raw %% 1 == 0), sum(raw)))
obj$patient <- sub("P.*$", "", obj$sample)   # 细胞名前缀 S2P1 → 患者 S2
cat("按细胞名推得的患者数:", length(unique(obj$patient)), "\n")
md <- obj@meta.data
md$cell <- rownames(md)
fwrite(md[, c("cell", "sample", "patient", "seurat_clusters", "celltype", "RiskModule1", "RiskWeighted",
              "RiskPos1", "RiskNeg1", "percent.mt", "nFeature_RNA")],
       "04_results/SC/risk_module_scores.txt", sep = "\t", quote = FALSE)

ctb <- tapply(md$RiskWeighted, md$celltype, mean)
ctab <- data.frame(celltype = names(ctb),
                   mean_module = as.numeric(tapply(md$RiskModule1, md$celltype, mean)),
                   median_module = as.numeric(tapply(md$RiskModule1, md$celltype, median)),
                   mean_weighted = as.numeric(ctb),
                   median_weighted = as.numeric(tapply(md$RiskWeighted, md$celltype, median)),
                   mean_riskpos = as.numeric(tapply(md$RiskPos1, md$celltype, mean)),
                   mean_riskneg = as.numeric(tapply(md$RiskNeg1, md$celltype, mean)),
                   n_cells = as.numeric(table(md$celltype)[names(ctb)]))
write.csv(ctab, "04_results/SC/risk_module_by_celltype.csv", row.names = FALSE)
cat("各细胞类型评分均值（无权重模块 / Cox 加权）:\n"); print(ctab, digits = 3, row.names = FALSE)

## ---------- 图 9A UMAP ----------
pA <- DimPlot(obj, group.by = "celltype", label = TRUE, repel = TRUE, pt.size = 0.35,
              cols = c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9",
                       "#8C564B", "#7F7F7F", "#BCBD22", "#17BECF")) +
  labs(x = "UMAP 1", y = "UMAP 2") + theme(legend.position = "right", legend.text = element_text(size = 8))
save_pdf(pA, "03_figures/Fig9A_scRNA_umap.pdf", 6.8, 5.4)

## ---------- 图 9B 15 基因在各细胞类型表达（点图）----------
pB <- DotPlot(obj, features = GENES, group.by = "celltype", cols = c("grey90", "#D55E00")) +
  labs(x = NULL, y = NULL) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                                   axis.text.y = element_text(size = 8))
save_pdf(pB, "03_figures/Fig9B_scRNA_model_genes.pdf", 7.2, 4.6)

## ---------- 图 9C 风险模块评分按细胞类型 ----------
ord <- names(sort(tapply(md$RiskWeighted, md$celltype, mean), decreasing = TRUE))
md$celltype2 <- factor(md$celltype, levels = ord)
pC <- ggplot(md, aes(celltype2, RiskWeighted, fill = celltype2)) +
  geom_violin(scale = "width", linewidth = 0.25, colour = "grey30") +
  geom_boxplot(width = 0.14, outlier.size = 0.15, linewidth = 0.25, fill = "white") +
  scale_fill_manual(values = rep(OI, length.out = length(ord)), guide = "none") +
  labs(x = NULL, y = "Weighted risk score (Cox coefficients)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
save_pdf(pC, "03_figures/Fig9C_scRNA_risk_module.pdf", 7.2, 4.6)

## ---------- 图 9D 肿瘤细胞风险模块 vs 增殖/髓系 M1M2 标记 ----------
mk <- list(Proliferation = c("MKI67", "TOP2A", "PCNA", "CCNB1"),
           `M1 marker` = c("CD86", "NOS2", "IL1B", "TNF", "CXCL9"),
           `M2 marker` = c("CD163", "MRC1", "MSR1", "IL10", "TGFB1"),
           `T cell exhaustion` = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4"),
           `Antigen presentation` = c("HLA-DRA", "HLA-DRB1", "CD74", "B2M", "TAP1"))
## 注意：名字里必须保留数字，否则 "M1 marker"/"M2 marker" 都被削成 "Mmarker" 互相覆盖
for (k in names(mk)) obj <- AddModuleScore(obj, features = list(mk[[k]]), name = gsub("[^A-Za-z0-9]", "", k), seed = 1)
md2 <- obj@meta.data
md2$cell <- rownames(md2)
sub <- md2[as.character(md2$celltype) == "Tumourcell", ]
sig_names <- c("Proliferation", "M1marker", "M2marker", "Tcellexhaustion", "Antigenpresentation")
cors <- do.call(rbind, lapply(sig_names, function(k) {
  col <- if (paste0(k, "1") %in% colnames(sub)) paste0(k, "1") else if (k %in% colnames(sub)) k else NA
  if (is.na(col)) return(NULL)
  ct <- suppressWarnings(cor.test(sub$RiskWeighted, sub[[col]], method = "spearman"))
  data.frame(Signature = k, rho = as.numeric(ct$estimate), p = ct$p.value, n = nrow(sub))
}))
cors <- as.data.frame(cors)
if (nrow(cors) == 0) stop("相关性计算无结果")
cors <- cors[!is.na(cors$rho), ]
cors$Signature <- c(Proliferation = "Proliferation", M1marker = "M1 marker", M2marker = "M2 marker",
                     Tcellexhaustion = "T cell exhaustion",
                     Antigenpresentation = "Antigen presentation")[cors$Signature]
cors <- cors[!is.na(cors$Signature), ]
fwrite(cors, "04_results/SC/risk_module_correlation.txt", sep = "\t", quote = FALSE)
print(cors, row.names = FALSE, digits = 3)

## ---------- 患者层面分析（肿瘤细胞均值；2190 个细胞来自 14 位患者，细胞级 P 值属伪重复） ----------
pt <- aggregate(cbind(RiskWeighted, RiskModule1, RiskPos1, RiskNeg1, Proliferation1, Antigenpresentation1) ~ patient,
                data = sub, FUN = mean)
pt$n_cells <- as.numeric(table(sub$patient)[pt$patient])
say_pt <- function(a, b) {
  ct <- suppressWarnings(cor.test(pt[[a]], pt[[b]], method = "spearman"))
  sprintf("%s vs %s: rho = %.3f, P = %.3g (n = %d 位患者)", a, b, ct$estimate, ct$p.value, nrow(pt))
}
ptl <- c(sprintf("患者层面汇总（肿瘤细胞均值，%d 位患者）", nrow(pt)),
         say_pt("RiskWeighted", "Proliferation1"), say_pt("RiskWeighted", "Antigenpresentation1"),
         say_pt("RiskModule1", "Proliferation1"),
         if (requireNamespace("lme4", quietly = TRUE)) {
           mm <- tryCatch(lme4::lmer(Proliferation1 ~ RiskWeighted + (1 | patient), data = sub),
                          error = function(e) NULL)
           if (is.null(mm)) "混合模型未收敛，未报告"
           else {
             cs <- summary(mm)$coefficients
             if (nrow(cs) >= 2) {
               tcol <- if ("t value" %in% colnames(cs)) "t value" else colnames(cs)[ncol(cs)]
               sprintf("混合模型 Proliferation ~ RiskWeighted + (1|patient): 斜率 = %.4f, t = %.3f（lme4 不直接给 P 值）",
                       cs[2, 1], cs[2, tcol])
             }
             else "混合模型固定效应秩不足，未报告"
           }
         } else "（未安装 lme4，未做混合模型）")
print(pt, digits = 3, row.names = FALSE)
cat(paste(ptl, collapse = "\n"), "\n")
fwrite(pt, "04_results/SC/risk_score_patient_level.txt", sep = "\t", quote = FALSE)
writeLines(ptl, "04_results/SC/risk_score_patient_level_report.txt")

## ---------- 正/负系数基因子集（回应外部意见：无权重平均会抵消方向） ----------
say_pn <- function(a, b, dat) {
  ct <- suppressWarnings(cor.test(dat[[a]], dat[[b]], method = "spearman"))
  sprintf("%s vs %s: rho = %+.3f, P = %.3g (n = %d)", a, b, ct$estimate, ct$p.value, nrow(dat))
}
pn <- c(sprintf("正系数基因 %d 个: %s", length(gp), paste(gp, collapse = ", ")),
        sprintf("负系数基因 %d 个: %s", length(gn), paste(gn, collapse = ", ")),
        sprintf("RiskWeighted vs RiskPos1: rho = %+.3f", cor(sub$RiskWeighted, sub$RiskPos1, method = "spearman")),
        sprintf("RiskWeighted vs RiskNeg1: rho = %+.3f", cor(sub$RiskWeighted, sub$RiskNeg1, method = "spearman")),
        sprintf("RiskPos1 vs RiskNeg1（肿瘤细胞内）: rho = %+.3f", cor(sub$RiskPos1, sub$RiskNeg1, method = "spearman")),
        "肿瘤细胞内（细胞级）:", say_pn("RiskPos1", "Proliferation1", sub), say_pn("RiskNeg1", "Proliferation1", sub),
        "肿瘤细胞内（患者层面，n=14）:", say_pn("RiskPos1", "Proliferation1", pt), say_pn("RiskNeg1", "Proliferation1", pt))
cat(paste(pn, collapse = "\n"), "\n")
writeLines(pn, "04_results/SC/risk_pos_neg_modules.txt")
pD <- ggplot(cors, aes(reorder(Signature, rho), rho, fill = rho)) +
  geom_col(width = 0.62, colour = "grey25", linewidth = 0.3) +
  geom_text(aes(label = sprintf("rho = %.2f\nP = %.1e", rho, p)), vjust = ifelse(cors$rho > 0, -0.15, 1.15), size = 2.7) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", guide = "none") +
  labs(x = NULL, y = "Spearman rho with the coefficient-weighted risk score (tumour cells)") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))
save_pdf(pD, "03_figures/Fig9D_scRNA_correlation.pdf", 6.6, 4.4)
cat("SCRNA DONE\n")
