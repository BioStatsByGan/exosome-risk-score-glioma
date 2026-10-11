# -*- coding: utf-8 -*-
"""补齐 batch_000.tar.gz 截断导致的缺失 MAF（2026-10-09 修复）
背景：dl_gdc_maf.py 下载后解包时 batch_000.tar.gz 报
      "Compressed file ended before the end-of-stream marker was reached"，
      994 个清单文件仅落地 926 个。本脚本按 00_data/maf/file_list.json 与已解包目录比对，
      只重新下载缺失的文件并解包到同一目录结构。
输出：00_data/maf/batch_missing.tar.gz、解包到 00_data/maf/extracted/、日志 download_missing.log
"""
import json, os, tarfile, time, urllib.request

ROOT = r"F:\0000_Bioinfo\222_tcgaSHAP_胶质瘤\00_data\maf"
EXTR = os.path.join(ROOT, "extracted")
LOG = open(os.path.join(ROOT, "download_missing.log"), "a", encoding="utf-8")

def log(m):
    LOG.write(time.strftime("[%H:%M:%S] ") + m + "\n"); LOG.flush(); print(m, flush=True)

j = json.load(open(os.path.join(ROOT, "file_list.json"), encoding="utf-8"))
ids = j["ids"]
have = {d for d in os.listdir(EXTR) if os.path.isdir(os.path.join(EXTR, d))}
missing = [i for i in ids if i not in have]
log(f"清单 {len(ids)} 个，已解包 {len(have)} 个，缺失 {len(missing)} 个")
if missing:
    out = os.path.join(ROOT, "batch_missing.tar.gz")
    data = json.dumps({"ids": missing}).encode()
    req = urllib.request.Request("https://api.gdc.cancer.gov/data", data=data,
                                headers={"Content-Type": "application/json",
                                         "User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=1800) as r, open(out, "wb") as f:
        total = 0
        while True:
            b = r.read(1 << 20)
            if not b:
                break
            f.write(b); total += len(b)
    log(f"[done] {os.path.basename(out)} {total/1048576:.2f} MB")
    with tarfile.open(out) as t:
        t.extractall(EXTR)
    log("解包完成")
have2 = {d for d in os.listdir(EXTR) if os.path.isdir(os.path.join(EXTR, d))}
log(f"补齐后已解包 {len(have2)} 个（清单 {len(ids)} 个，仍缺 {len(set(ids) - have2)} 个）")
log("DONE")