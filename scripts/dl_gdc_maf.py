"""从 GDC 批量下载 TCGA-GBM/LGG 的 masked somatic mutation MAF（用于 TMB 计算）。
输出：00_data/maf/<file>.maf.gz，日志 00_data/maf/download.log
"""
import json, os, sys, tarfile, time, urllib.parse, urllib.request

ROOT = "/mnt/f/0000_Bioinfo/222_tcgaSHAP_胶质瘤/00_data/maf"
os.makedirs(ROOT, exist_ok=True)
LOG = open(os.path.join(ROOT, "download.log"), "a", encoding="utf-8")


def log(msg):
    LOG.write(time.strftime("[%H:%M:%S] ") + msg + "\n")
    LOG.flush()
    print(msg, flush=True)


FILT = {"op": "and", "content": [
    {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-GBM", "TCGA-LGG"]}},
    {"op": "in", "content": {"field": "data_category", "value": ["Simple Nucleotide Variation"]}},
    {"op": "in", "content": {"field": "data_type", "value": ["Masked Somatic Mutation"]}},
    {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["Aliquot Ensemble Somatic Variant Merging and Masking"]}},
    {"op": "in", "content": {"field": "cases.samples.sample_type", "value": ["Primary Tumor", "Recurrent Tumor"]}},
]}

ids, samples = [], []
frm = 0
while True:
    qs = urllib.parse.urlencode({"filters": json.dumps(FILT), "fields": "file_id,file_name,cases.samples.submitter_id",
                                 "format": "JSON", "size": "500", "from": str(frm)})
    with urllib.request.urlopen(f"https://api.gdc.cancer.gov/files?{qs}", timeout=90) as r:
        d = json.load(r)
    hits = d["data"]["hits"]
    total = d["data"]["pagination"]["total"]
    for h in hits:
        ids.append(h["file_id"])
        try:
            samples.append(h["cases"][0]["samples"][0]["submitter_id"])
        except Exception:
            samples.append("")
    frm += len(hits)
    log(f"已列出 {frm}/{total}")
    if frm >= total or not hits:
        break

json.dump({"ids": ids, "samples": samples}, open(os.path.join(ROOT, "file_list.json"), "w"))
log(f"MAF 文件总数 {len(ids)}")

BATCH = 120
for i in range(0, len(ids), BATCH):
    chunk = ids[i:i + BATCH]
    out = os.path.join(ROOT, f"batch_{i//BATCH:03d}.tar.gz")
    if os.path.exists(out) and os.path.getsize(out) > 10000:
        log(f"[skip] {os.path.basename(out)}")
        continue
    data = json.dumps({"ids": chunk}).encode()
    req = urllib.request.Request("https://api.gdc.cancer.gov/data", data=data,
                                 headers={"Content-Type": "application/json",
                                          "User-Agent": "Mozilla/5.0"})
    try:
        with urllib.request.urlopen(req, timeout=1800) as r, open(out, "wb") as f:
            while True:
                b = r.read(1 << 20)
                if not b:
                    break
                f.write(b)
        log(f"[done] {os.path.basename(out)} {os.path.getsize(out)/1048576:.1f} MB")
    except Exception as e:
        log(f"[FAIL] {os.path.basename(out)} {type(e).__name__} {str(e)[:80]}")

for f in sorted(os.listdir(ROOT)):
    if f.endswith(".tar.gz"):
        try:
            with tarfile.open(os.path.join(ROOT, f)) as t:
                t.extractall(os.path.join(ROOT, "extracted"))
            log(f"解包 {f}")
        except Exception as e:
            log(f"解包失败 {f} {e}")
log("MAF DOWNLOAD FINISHED")
