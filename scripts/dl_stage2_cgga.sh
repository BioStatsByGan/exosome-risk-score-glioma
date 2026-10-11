#!/bin/bash
# 胶质瘤项目 Stage 2 下载：CGGA（验证队列 + 正常对照 + 单细胞）
# CGGA 下载链接是查询式：/download?file=download/<日期>/<文件名>&type=<类型>&time=<日期>
set -u
OUT="/mnt/f/0000_Bioinfo/222_tcgaSHAP_胶质瘤/00_data"
mkdir -p "$OUT"; cd "$OUT"
LOG="$OUT/download_stage2.log"
echo "=== STAGE2 START $(date) ===" >> "$LOG"
UA="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 Chrome/120 Safari/537.36"
REF="http://www.cgga.org.cn/download.jsp"

get() {
  local url="$1" out="$2"
  if [ -s "$out" ] && [ "$(stat -c%s "$out")" -gt 10000 ]; then
    echo "[skip] $out $(du -h "$out" | cut -f1)" >> "$LOG"; return 0
  fi
  echo "[get ] $out" >> "$LOG"
  curl -sL -C - --retry 5 --retry-delay 5 --max-time 7200 \
       -A "$UA" -e "$REF" -o "$out.part" "$url" >> "$LOG" 2>&1
  if [ -s "$out.part" ] && [ "$(stat -c%s "$out.part" 2>/dev/null || echo 0)" -gt 10000 ]; then
    mv "$out.part" "$out"; echo "[done] $out $(du -h "$out" | cut -f1)" >> "$LOG"
  else
    echo "[FAIL] $out  size=$(stat -c%s "$out.part" 2>/dev/null || echo 0)" >> "$LOG"; rm -f "$out.part"
  fi
}

B="http://www.cgga.org.cn/download?file=download"

# 验证队列：mRNAseq_325 (FPKM) + 临床；693 做敏感性分析
get "$B/20200506/CGGA.mRNAseq_325.RSEM-genes.20200506.txt.zip&type=mRNAseq_325&time=20200506"                 "CGGA_mRNAseq_325_FPKM.zip"
get "$B/20200506/CGGA.mRNAseq_325_clinical.20200506.txt.zip&type=mRNAseq_325_clinical&time=20200506"           "CGGA_mRNAseq_325_clinical.zip"
get "$B/20200506/CGGA.mRNAseq_693.RSEM-genes.20200506.txt.zip&type=mRNAseq_693&time=20200506"                  "CGGA_mRNAseq_693_FPKM.zip"
get "$B/20200506/CGGA.mRNAseq_693_clinical.20200506.txt.zip&type=mRNAseq_693_clinical&time=20200506"           "CGGA_mRNAseq_693_clinical.zip"
# 正常脑对照 20 例
get "$B/20200614/CGGA_RNAseq_Control_20.txt.zip&type=CGGA_RNAseq_Control_20&time=20200614"                     "CGGA_normal20_FPKM.zip"
# 单细胞（16 例 IDH 突变胶质瘤，6148 细胞 UMI 矩阵 + 临床注释）
get "$B/20220104/CGGA.scRNAseq_6148.count.matrix.tsv.zip&type=scRNAseq_6148&time=20220104"                     "CGGA_scRNAseq_6148_counts.zip"
get "$B/20220104/CGGA.scRNAseq_6148_clinical.20220104.txt.zip&type=scRNAseq_6148_clinical&time=20220104"       "CGGA_scRNAseq_6148_clinical.zip"
# 新增一版（2025 发布，IDH 突变队列，做第二验证）
get "$B/20250915/CGGA_IDH_A_RNAseq_RSEM_20250915.zip&type=RNAseq_RSEM_counts&time=20250915"                    "CGGA_IDH_A_RNAseq_RSEM.zip"
get "$B/20250915/CGGA_IDH_A_Clinical_Information.zip&type=Clinical_Information_counts&time=20250915"           "CGGA_IDH_A_clinical.zip"

echo "=== STAGE2 END $(date) ===" >> "$LOG"
du -sh "$OUT" >> "$LOG"
echo "STAGE2 FINISHED" >> "$LOG"
