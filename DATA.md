# Data sources

None of the individual-level data are redistributed in this repository. Everything below is public; the
scripts in `scripts/` download or reconstruct the required inputs. Please follow each provider's access terms.

| Dataset | Used for | Source | Notes |
|---|---|---|---|
| TCGA-GBM / TCGA-LGG | discovery cohort (n=804 with expression; 799 with survival) | GDC / UCSC Xena (log2(TPM+1) and clinical) | scripts `dl_stage1.sh`, `00_prep_inputs.R` |
| GDC masked somatic mutation MAF | tumour mutational burden | GDC via `dl_gdc_maf.py` / `dl_gdc_maf_missing.py` | 994 files, 887 patients; one file per patient retained |
| GTEx v8 brain | normal-brain comparator for differential expression | GTEx portal | exploratory only: source and group are fully confounded, no batch correction applied |
| CGGA mRNAseq_325 / mRNAseq_693 | external validation | CGGA data portal | RSEM-FPKM; downloaded by `dl_stage2_cgga.sh` |
| CGGA scRNA-seq (6148 cells) | single-cell mapping of the score | CGGA single-cell resource | multi-sector biopsies, 14 patients / 75 sample points |
| GTEx v8 brain cis-eQTL | MR instruments (primary) | GTEx / SMR-format summary data | instruments = lead cis-eQTL per gene and tissue |
| eQTLGen full cis (100 kb) | MR cross-tissue replication and full-cis-window colocalization | eQTLGen (n=30,895) | 19,156,973 rows; GRCh37 - merge with FinnGen by rsID, not by position |
| FinnGen R12 | MR outcome (C3_GBM_ASTROCYTOMA_EXALLC; malignant neoplasm of brain) | FinnGen R12 public summary statistics | 544 cases / 378,749 controls (primary) |
| CIBERSORT LM22 | immune fractions | CIBERSORT | 22 fractions, 1,000 permutations |
| MCP-counter, 29-immune-feature ssGSEA, brain-aware signatures | immune deconvolution Comparison | package data in `scripts/` | see `38_mcpcounter.R`, `29_immFunction.R`, `audit_brain_immune.R` |
| MSigDB Hallmark | pathway activity | MSigDB | 50 gene sets, GSVA with Gaussian kernel |
| GDSC2 via oncoPredict | predicted drug sensitivity | GDSC2 | 805 cell lines, 198 drugs |
| TIDEpy | predicted immune checkpoint blockade response | TIDE / TIDEpy | binary response uninformative in glioma; only continuous components reported |
| ESTIMATE | tumour purity sensitivity analysis | `estimate` R package / ESTIMATE score | cosine-transformed purity |
| TCGA pan-cancer ABSOLUTE purity | DNA-based purity covariate | TCGA pan-cancer atlas annotation | available for 675 of the 799 tumours |

## Reference data that must be supplied by the user

Path `H:/...` and `I:/...` strings inside a few scripts are the authors' local copies of large public
resources (GTEx besd files, eQTLGen cis matrices, FinnGen summary statistics). Replace them with your own
local paths; the file names are given in `scripts/README.md`.
