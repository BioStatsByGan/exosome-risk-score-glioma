# Supplementary Code — 222.tcgaSHAP 胶质瘤

本目录为论文的分析脚本（与运行顺序说明），对应 Supplementary Code。

## 环境

- R **4.4.2**（唯一验证版本；`C:\Program Files\R\R-4.4.2\bin\Rscript.exe`）。主要依赖：data.table、survival、survminer、glmnet、rms、pec、forestploter、ggplot2、patchwork、ComplexHeatmap、circlize、GSVA、GSEABase、limma、sva、clusterProfiler、org.Hs.eg.db、enrichplot、Seurat v5、harmony、MCPcounter、timeROC、oncoPredict、shapviz/kernelshap、coloc、TwoSampleMR/EQLDtools。
- Python 3.11+，依赖：pandas、numpy、python-docx、PyMuPDF(fitz)、requests/urllib（标准库）。
- 部分脚本（下载类）内写的是 WSL 路径 `/mnt/f/...`；R 脚本一律使用 Windows 路径 `F:/...`。

## 运行前必读（两个坑）

1. **Windows 下 R 的中文路径**：若环境里存在 `LANG/LC_ALL/LC_CTYPE=C.UTF-8`（部分终端/工具会注入），R 的 locale 会退化为 `C`，`setwd("F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤")` 直接失败、`list.files()` 返回 0。运行前清空：
   ```powershell
   $env:LC_ALL=$null; $env:LC_CTYPE=$null; $env:LC_COLLATE=$null; $env:LC_TIME=$null; $env:LC_MONETARY=$null; $env:LANG=$null
   ```
2. **下载类步骤要核对「清单数 vs 落地数」**：`dl_gdc_maf.py` 下载完成后解包，若某个 tar 分片被截断，日志里下载仍记为完成而解包失败，后续按目录扫文件不会报错，会静默少样本。补齐脚本：`dl_gdc_maf_missing.py`（比对 `file_list.json` 与 `00_data/maf/extracted/`，只重下缺失文件）。

## 数据准备（可跳过：数据均可从公开资源获取）

| 脚本 | 作用 |
|---|---|
| `dl_stage1.sh` | TCGA/GTEx/FinnGen 等 Stage-1 下载（curl，WSL 下运行） |
| `dl_stage2_cgga.sh` | CGGA 验证队列、正常对照与单细胞数据下载 |
| `dl_gdc_maf.py` / `dl_gdc_maf_missing.py` | GDC masked somatic mutation MAF 下载 / 缺失补齐（TMB 输入） |
| `mr_extract_eqtl.py`、`mr_get_outcome.py`、`mr_fix_coords.py` | MR 暴露（GTEx 脑 cis-eQTL）与结局（FinnGen R12）取数与坐标核对 |
| `fetch_refs.py`、`fetch_refs2.py`、`fetch_refs3.py`、`splice_refs.py`、`verify_refs_doi.py`、`check_doi_crossref.py` | 参考文献抓取与 DOI 核验 |

## 分析主链（按顺序）

| 顺序 | 脚本 | 产出 |
|---|---|---|
| 0 | `00_prep_inputs.R` | `02_inputs/`（表达、临床、基因集） |
| 1 | `13_diff_figs.R` | 差异表达与富集 → 图 1；`04_results/13.diff/` |
| 2 | `17_19_mergeTime_uniCox.R` | 生存合并与单因素 Cox → `19.uniCox/` |
| 3 | `20_shap.R` | LASSO-Cox + 逐步 Cox + SHAP → 图 2、`20.SHAP/` |
| 4 | `18_cgga_validate.R` | CGGA 外部验证 → 图 4、`18.CGGA/`（含 `validate_stats.txt`） |
| 5 | `21_24_risk_km_roc.R` | 训练队列风险分层与 ROC → 图 3 |
| 6 | `25_26_independence_clinical.R`、`26b_calibration.R` | 独立预后/临床相关性/列线图/校准 → 图 5、`25_26/` |
| 7 | `27_cibersort.R`…`35_boxplot.R` | 免疫浸润、免疫功能、免疫分型、TMB、TIDE、药敏 → 图 7、图 8、图 S2、图 S4 |
| 8 | `36_subgroup_survival.R` | 亚层分析与通路活性 → 图 6 |
| 9 | `37_hallmark_pathway.R`、`38_mcpcounter.R`、`39_neutrophil_resolve.R`、`40_idhwt_feasibility.R` | 通路、MCP-counter、中性粒方法一致性（图 S5）、IDH 野生型可行性 |
| 10 | `90_scrna.R` | 单细胞定位 → 图 9 |
| 11 | `91_mr.R` | MR 与共定位 → 图 S3 |
| 12 | `audit_bulk.R`、`audit_b_group.R`、`audit_brain_immune.R`、`audit_mr_sc.R`、`audit_cross_tissue_mr.py` 等 | 投稿前审计（C 指数增量、DCA/Brier、脑特异署名敏感性、跨组织 MR 等）→ `04_results/audit/` |
| 13 | `make_tables.R` | 表 1–3（CSV） |
| 14 | `assemble_figures.py`、`assemble_figures_v2.py` | 面板拼版（图 1–9、图 S1–S5；图 S6 由 `47_ext_calibration.R` 直接输出） |
| 15 | `assemble_en_v2.py`、`assemble_cn_v2.py` | 生成中/英文正文 markdown（含文献编号与表注） |
| 16 | `md2docx_en_figs.py`、`md2docx_cn_figs.py`、`md2docx_cover.py` | 生成 docx；`build_submission_pkg.py` 组装投稿包 |
| 17 | `41_nested_cv.R`、`42_screen_free_model.R`、`43_immune_adjusted.R`、`44_clinical_increment.R`、`45_purity_sensitivity.R`、`46_blood_fullcis_coloc.R`、`47_ext_calibration.R` | 审稿补充分析：嵌套交叉验证、跳过差异筛选重建模型、22 类免疫四口径校正、外部队列 ΔC、肿瘤纯度校正、血液端完整 cis 窗共定位、外部队列校准与预测误差 → 表 S6–S9、图 S6 |

## 复现要点

- 随机种子：`20_shap.R`(42)、`26b_calibration.R`(7)、`90_scrna.R`(20260920)、`40_idhwt_feasibility.R`(20260926)、校准曲线 bootstrap(200)、timeROC CI(n.sim=500)。
- 主要数值与脚本日志的对应关系见 `Supplementary_SourceData/MANIFEST.md`。
- 脚本以 UTF-8（无 BOM）保存；R 与 Python 脚本使用 LF 换行。
