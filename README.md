# exosome-risk-score-glioma

Analysis code and derived result tables for the manuscript

> **An exosome-related gene score for glioma: reproducible risk stratification, SHAP-based dissection of its drivers, and its myeloid-enriched tumour-state phenotype**
>
> 外泌体相关基因评分在胶质瘤中的应用：可复现的风险分层、SHAP 驱动的判别力来源解析及其髓系富集肿瘤状态表型

---

## What is in this repository / 本仓库内容

| Path | Contents |
|---|---|
| `scripts/` | 84 analysis scripts (R / Python / shell) for the whole pipeline, plus `scripts/README.md` with the detailed run order, environment and package list |
| `source_data/04_results/` | Derived result tables behind every figure and table (`04_results/` in the manuscript's Supplementary Source Data) |
| `source_data/tables/` | Table 1-3 and Table S1-S9 as CSV (ASCII file names; see `Tables_map.csv` below) |
| `source_data/MANIFEST.md` | Panel-by-panel mapping between figures/tables and their source files |

The 15-gene score is built in TCGA-GBM/LGG, validated in CGGA-325 and CGGA-693, and interpreted with SHAP;
downstream analyses cover Hallmark pathway activity, immune deconvolution, tumour mutational burden,
drug-sensitivity prediction, single-cell mapping, and germline Mendelian randomization with colocalization.

## Requirements

- **R 4.4.2** (the only validated version) - data.table, survival, survminer, glmnet, rms, pec, timeROC,
  ggplot2, patchwork, ComplexHeatmap, circlize, GSVA, GSEABase, limma, sva, clusterProfiler, org.Hs.eg.db,
  enrichplot, Seurat v5, harmony, MCPcounter, oncoPredict, kernelshap/shapviz, coloc, TwoSampleMR, EQLDtools
- **Python 3.11+** - pandas, numpy, python-docx, PyMuPDF, requests

```r
install.packages(c("data.table","survival","survminer","glmnet","rms","pec","timeROC","ggplot2",
                   "patchwork","ComplexHeatmap","circlize","GSVA","limma","sva","coloc"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GSEABase","clusterProfiler","org.Hs.eg.db","enrichplot"))
```

## How to run

`scripts/README.md` lists the full run order. In short: `00_prep_inputs.R` -> differential expression and
enrichment (`13_diff_figs.R`) -> survival merge and univariable Cox (`17_19_mergeTime_uniCox.R`) ->
LASSO-Cox + stepwise Cox + SHAP (`20_shap.R`) -> external validation (`18_cgga_validate.R`) ->
training-cohort stratification (`21_24_risk_km_roc.R`) -> independence/nomogram/calibration
(`25_26_independence_clinical.R`, `26b_calibration.R`) -> immune / TMB / TIDE / drug-sensitivity
(`27_cibersort.R` ... `35_boxplot.R`) -> subgroups and pathways (`36_subgroup_survival.R`,
`37_hallmark_pathway.R`, `38_mcpcounter.R`, `39_neutrophil_resolve.R`) -> single cell (`90_scrna.R`) ->
Mendelian randomization (`91_mr.R`) -> reviewer-requested supplementary analyses
(`41_nested_cv.R` ... `47_ext_calibration.R`), with audit scripts under `audit_*.R`.

## Reproducibility notes (four traps we hit)

1. **Windows + R + non-ASCII paths.** If `LANG`/`LC_ALL`/`LC_CTYPE` is set to `C.UTF-8`, R degrades its
   locale to `C` and `setwd()` into a path containing Chinese characters fails silently
   (`list.files()` returns 0). Clear them first:
   `$env:LC_ALL=$null; $env:LC_CTYPE=$null; $env:LC_COLLATE=$null; $env:LANG=$null`
2. **`survival::concordance(Surv ~ x)`** treats `x` as the survival time. Always use
   `summary(coxph(...))$concordance[1]`.
3. **`rms::cph()`** resolves `datadist` by name, so the `datadist()` object must live in the global
   environment (assign with `<<-`), otherwise `cph()` fails with "dataset dd not found".
4. **R 4.4 `sprintf()`** raises `too few arguments` when the format string contains non-ASCII characters
   *and* vector arguments are passed; keep such format strings ASCII.

## Data

No individual-level data are redistributed here; all datasets are public and are downloaded by the scripts
in `scripts/`. See `DATA.md` for sources, versions and access terms. The download steps must be run before
the analysis chain; derived tables in `source_data/` let you check the figures without re-running everything.

## Key results

| Analysis | Result |
|---|---|
| Training discrimination (TCGA) | C-index 0.824; multivariable HR 2.15 (1.72-2.68) |
| External validation | CGGA-325 C-index 0.731; CGGA-693 0.722 |
| Increment over a clinical model | dC 0.010 (TCGA), 0.032 (CGGA-325), 0.044 (CGGA-693), all bootstrap 95% CI excluding 0 |
| Calibration / prediction error | slopes 1.05-1.19 (1-year) and 1.46-1.65 (3-5 years) in CGGA; Brier/IPA in `source_data/tables/TableS9_external_calibration.csv` |
| Germline causality | no significant MR estimate (49 tests); max PP.H4 0.06 (brain, restricted) and 0.203 (blood, full cis window) |

## Citation

If you use this code, please cite the manuscript (see `CITATION.cff`). Fill in the author list and the DOI
after acceptance.

## License

MIT - see `LICENSE`.
