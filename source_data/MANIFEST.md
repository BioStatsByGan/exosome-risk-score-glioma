# MANIFEST — 图/表与源数据对应关系

补充源数据 = `04_results/`（分析中间与最终结果表；>5 MB 的中间矩阵、Seurat 对象与原始工具变量未纳入）+ `tables/`（表 1–3、表 S2–S5b 的 CSV）。

投稿包 `Tables/` 与 `tables/` 的关系分三类，逐类说明：

1. **Table S2–S5b**：与 `tables/` 中对应文件**逐字节一致**（直接复制）。
2. **Table 1**：`tables/表1a_基线_TCGA.csv` 与 `tables/表1b_基线_CGGA325.csv` 是两个队列各自的基线表（由 `make_tables.R` 生成）；投稿包内的 `Table1_baseline.csv` 是装配脚本 `assemble_en_v2.py` / `assemble_cn_v2.py` 把两张表**合并为一张**（英文标签、两队列并列）后的版本，数值逐项对应源文件。
3. **Table 2 / Table 3**：`tables/表2_多因素Cox_15基因.csv`、`tables/表3_独立预后分析.csv` 是中文标签的源表；投稿包内的 `Table2_multicox.csv` / `Table3_independence.csv` 是同一批数值的**英文标签版**（如 `riskScore` → `Risk score`、`TCGA univariable` → `TCGA, univariable`），标签与引号格式不同。

其中 **Table S1** 由装配脚本从 `04_results/audit/mr_power_analysis.txt` 生成（两文件按 TSV 制表符读取后转为 markdown 表）。

| 图/表面板 | 源文件 |
|---|---|
| 图 1A 火山图 | `04_results/13.diff/glioma_exosome_all.txt`（691 基因的 limma 输出） |
| 图 1B 热图（|log2FC| 前 50） | 同上（`glioma_exosome_diff.txt` 给出 422 个差异基因） |
| 图 1C GSEA 曲线 / 图 1D 富集网络 | `04_results/13.diff/KEGG.txt`、`GO.txt` |
| 图 2A LASSO-CV | `04_results/20.SHAP/lasso_genes.txt`（29 基因与系数） |
| 图 2B 森林图 / **表 2** | `04_results/20.SHAP/multiCox.txt` → `tables/表2_多因素Cox_15基因.csv`（投稿包 `Tables/Table2_multicox.csv` 为其英文标签版） |
| **表 1** 基线特征（TCGA 与 CGGA-325 按风险组） | `tables/表1a_基线_TCGA.csv`、`tables/表1b_基线_CGGA325.csv`（由 `make_tables.R` 从 `04_results/20.SHAP/risk.TCGA.txt`、`04_results/18.CGGA/risk.CGGA325.txt` 与临床表生成）；投稿包 `Tables/Table1_baseline.csv` 为装配脚本合并后的英文版 |
| **表 3** 单/多因素独立预后分析 | `04_results/25_26/independence_univariable_TCGA.txt`、`independence_multivariable_TCGA.txt`、`independence_multivariable_CGGA325.txt` → `tables/表3_独立预后分析.csv`（`make_tables.R` 生成；投稿包 `Tables/Table3_independence.csv` 为其英文标签版） |
| 图 2C–D SHAP 重要性/蜂群 | `04_results/20.SHAP/risk.TCGA.txt`（含 15 基因表达与 riskScore） |
| 图 S1 瀑布/力图 | 同上（脚本 20_shap.R 以单个患者示例绘制） |
| 图 3A–B 风险分布/生存状态 | `04_results/20.SHAP/risk.TCGA.txt` |
| 图 3C KM / 图 3D 时间 ROC | `04_results/21_24/train_stats.txt`（C 指数、HR、1/3/5 年 AUC） |
| 图 4 外部验证 | `04_results/18.CGGA/risk.CGGA325.txt`、`risk.CGGA693.txt`、`validate_stats.txt` |
| 图 5A/5A2 独立预后 | `04_results/25_26/independence_univariable_TCGA.txt`、`independence_multivariable_TCGA.txt` |
| 图 5B–D 临床相关性/热图/列线图 | `04_results/25_26/model_and_calibration_stats.txt`、`calibration.rds` |
| 图 5E 校准曲线 | `04_results/25_26/calibration.rds`、`model_and_calibration_stats.txt` |
| 图 6A/B 亚层 HR 与层内 KM；**表 S2** | `04_results/36.subgroup/subgroup_continuous_TCGA.txt`、`subgroup_continuous_CGGA325.txt` → `tables/表S2_分层验证.csv`（与投稿包 `Tables/TableS2_subgroup.csv` 逐字节一致） |
| 图 6C/D Hallmark 通路；**表 S3** | `04_results/37.hallmark/hallmark_diff.txt` → `tables/表S3_Hallmark通路.csv`（与投稿包 `Tables/TableS3_hallmark.csv` 逐字节一致） |
| 图 7A/B 29 项 ssGSEA | `04_results/29.immFunction/ssgsea_diff.txt` |
| 图 7C CIBERSORT 22 细胞；表 S5 输入 | `04_results/27.CIBERSORT/CIBERSORT-Results.txt`、`28.immDiff/imm_cell_diff.txt` |
| 图 7D 免疫分型 | `04_results/30.immSubtype/immune_subtype_risk.txt` |
| 图 7E MCP-counter；**表 S4** | `04_results/38.MCPcounter/mcp_diff.txt` → `tables/表S4_MCPcounter.csv`（与投稿包 `Tables/TableS4_MCPcounter.csv` 逐字节一致） |
| 图 8A/B TMB 与四组生存 | `04_results/31.riskTMB/TMB.txt`、`tmb_risk.txt`、`tmb_surv.txt`、`run.log` |
| 图 8C 药敏火山 / 图 S4 药敏箱线 | `04_results/35.drugSensitivity/drug_diff.txt`、`34.oncoPredict/calcPhenotype_Output/DrugPredictions.csv` |
| 图 8D 决策曲线 | `04_results/audit/dca_3year.txt`；Brier/IPA：`audit/brier_scores.txt` |
| 图 9A–D 单细胞 | `04_results/SC/cell_annotation.txt`、`risk_module_scores.txt`、`risk_module_by_celltype.csv`、`risk_module_correlation.txt`、`audit/scrna_module_null_test.txt` |
| 图 S2 TIDE | `04_results/33.TIDE/tide_diff.txt`、`TIDE.csv` |
| 图 S3 MR 与共定位 | `04_results/MR/mr_results.txt`、`coloc_results.txt`、`eqtl_brain_summary.txt`；跨组织：`audit/cross_tissue_mr_blood.txt` |
| **表 S1** MR 把握度（装配脚本按 TSV 读入后生成） | `04_results/audit/mr_power_analysis.txt`、`mr_power_min_or.txt` |
| **图 S5** 与 **表 S5 / 表 S5b** 中性粒六法（表与投稿包同名 CSV 逐字节一致） | `04_results/39.neutrophil/neutrophil_resolution.txt`、`method_cor_matrix.txt`；脑特异署名：`audit/brain_immune_ssgsea_sensitivity.txt`、`audit/brain_aware_neutrophil_scores.txt` |
| 审计类数字（C 指数增量/R²/TMB 分级校正/免疫校正） | `04_results/audit/audit_bulk_log.txt` |
| 审计类数字（DCA 例数/Brier/切点敏感性） | `04_results/audit/audit_b_group_log.txt` |
| **表 S6** 22 类 CIBERSORT 三种校正口径 + 免疫含量校正 | `04_results/audit/immune_all_cells_adjusted.txt`（投稿包 `Tables/TableS6_immune_cells_adjusted.csv` 与其一致） |
| **表 S8** 血液端完整 cis 窗共定位（eQTLGen 100 kb 全位点 × FinnGen R12） | `04_results/audit/blood_fullcis_coloc.txt`、`blood_fullcis_coloc_log.txt` |
| 纯度校正敏感性分析（ESTIMATE 799 例；DNA 来源 ABSOLUTE 675 例；TMB/免疫/通路/药敏四组模型） | `04_results/audit/purity_table.txt`、`purity_sensitivity.txt`、`immune_all_cells_adjusted_purity.txt` |
| **表 S7** 三队列「临床 vs 临床+评分」C 指数与 ΔC 的 bootstrap 95%CI | `04_results/audit/clinical_increment.txt`（投稿包 `Tables/TableS7_clinical_increment.csv` 与其一致） |
| **图 S6** 与 **表 S9** 外部队列校准与预测误差（三个队列；TCGA 一列与正文 §3.10 同值） | `04_results/audit/external_calibration.txt` → `tables/表S9_外部队列校准与预测误差.csv`（投稿包 `Tables/TableS9_external_calibration.csv` 与其一致；脚本 `47_ext_calibration.R`） |
| 单细胞正/负系数基因子集（RiskPos/RiskNeg） | `04_results/SC/risk_pos_neg_modules.txt`、`risk_module_by_celltype.csv` |
| 参考文献内容支持性核对（Europe PMC 摘要/开放获取全文；脚本见 `refcheck_*.py`） | `04_results/audit/ref_supportability.txt`、`ref_supportability_detail.txt`、`ref_supportability_numbers.txt`、`ref_verify.json` |
| **表 S10** 固定系数 13 基因评分对照（TCGA/两 CGGA；含 TPM vs FPKM 尺度） | `04_results/audit/gene13_fixed.txt` → `tables/表S10_13基因固定系数对照.csv` |
| **表 S11** CGGA 原发肿瘤限定分析 | `04_results/audit/primary_only.txt` → `tables/表S11_原发肿瘤限定分析.csv` |
| **表 S12** 各基因平均 |SHAP| 与留一/置换 ΔC | `04_results/audit/shap_ablation.txt` → `tables/表S12_SHAP消融与置换.csv` |
| **图 S7** 与 **表 S13** IPCW 决策曲线、校准截距、共线性与子集选择检查 | `04_results/audit/dca_ipcw.txt`、`calibration_diag.txt`、`collinearity.txt` → `tables/表S13_校准与DCA诊断.csv` |
| **表 S14** 分析样本流程（每个分析的分母） | `04_results/audit/sample_flow.txt` → `tables/表S14_分析样本流程.csv` |
| **表 S15** TMB 文件选择规则敏感性 | `04_results/audit/tmb_rule_sensitivity.txt` → `tables/表S15_TMB文件规则敏感性.csv` |
| 预测生存率量纲核对（Brier 修正依据） | `04_results/audit/prediction_scale_check.txt` |
| 图 8D 决策曲线（含仅临床变量模型） | `04_results/audit/dca_3year.txt`（列：model/clinical/all/none） |
| 图 9C/9D 与患者层面统计（加权评分、14 位患者） | `04_results/SC/risk_module_scores.txt`、`risk_module_by_celltype.csv`、`risk_module_correlation.txt`、`risk_score_patient_level.txt`、`risk_score_patient_level_report.txt` |
| 嵌套交叉验证（全流程乐观度 0.022，各折 11–29 基因） | `04_results/audit/nested_cv.txt`、`nested_cv_log.txt` |
| 跳过肿瘤–正常差异筛选的敏感性建模（22 基因/C 0.843/rho 0.923）与 13 基因缩减模型（C 0.821） | `04_results/audit/screen_free_model.txt`、`screen_free_model_log.txt` |
| 模型与校准（TCGA/CGGA-325 多因素、校准斜率） | `04_results/25_26/model_and_calibration_stats.txt` |

注：`04_results/` 内另有若干 >5 MB 的中间文件未纳入（CIBERSORT 的 uniq.symbol 矩阵、Seurat 对象 `cgga_scRNA_seurat.rds`、eQTLGen 原始工具变量、GDSC2 表达矩阵、deg 缓存），它们由公开数据与脚本可完整重建。
