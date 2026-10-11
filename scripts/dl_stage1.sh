#!/bin/bash
# 胶质瘤项目 Stage 1 数据下载（幂等：已存在且大小>0 则跳过；支持断点续传）
set -u
export https_proxy=http://127.0.0.1:7897 http_proxy=http://127.0.0.1:7897
OUT="/mnt/f/0000_Bioinfo/222_tcgaSHAP_胶质瘤/00_data"
mkdir -p "$OUT"
cd "$OUT"
LOG="$OUT/download.log"
echo "=== START $(date) ===" >> "$LOG"

get() {  # get <url> <目标文件名>
  local url="$1" out="$2"
  if [ -s "$out" ]; then echo "[skip] $out 已存在 $(du -h "$out" | cut -f1)" | tee -a "$LOG"; return 0; fi
  echo "[get ] $out  <- $url" | tee -a "$LOG"
  curl -L -C - --retry 5 --retry-delay 5 --max-time 7200 -o "$out.part" "$url" >> "$LOG" 2>&1
  if [ -s "$out.part" ]; then mv "$out.part" "$out"; echo "[done] $out $(du -h "$out" | cut -f1)" | tee -a "$LOG";
  else echo "[FAIL] $out" | tee -a "$LOG"; fi
}

# 1) TCGA + GTEx 统一处理矩阵（TOIL RSEM TPM，含肿瘤与正常脑）
get "https://toil.xenahubs.net/download/TcgaTargetGtex_rsem_gene_tpm.gz"        "TcgaTargetGtex_rsem_gene_tpm.gz"
get "https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz"        "TcgaTargetGTEX_phenotype.txt.gz"
get "https://toil.xenahubs.net/download/TCGA_survival_data"                     "TCGA_survival_data"

# 2) 胶质瘤 GWAS（FinnGen R12，结局端）
get "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_C3_GBM_ASTROCYTOMA_EXALLC.gz" "finngen_R12_C3_GBM_ASTROCYTOMA_EXALLC.gz"
get "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_C3_BRAIN_EXALLC.gz"            "finngen_R12_C3_BRAIN_EXALLC.gz"

# 3) 突变数据（TMB 与 IDH 状态）
get "https://toil.xenahubs.net/download/mc3.v0.2.8.PUBLIC.xena.gz"               "mc3.v0.2.8.PUBLIC.xena.gz"

echo "=== END $(date) ===" >> "$LOG"
du -sh "$OUT" >> "$LOG"
echo "ALL DOWNLOADS FINISHED" | tee -a "$LOG"
