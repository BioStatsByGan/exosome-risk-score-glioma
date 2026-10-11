"""把正文里的 [REF:PMID] 占位符按首次出现顺序编号，并生成 GB/T 7714 风格参考文献清单。"""
import json, os, re, sys, urllib.parse, urllib.request

ROOT = "/mnt/f/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
MD = os.path.join(ROOT, "稿", "手稿_中文_v1.md")
OUT_MD = os.path.join(ROOT, "稿", "手稿_中文_v2_带文献.md")
OUT_REF = os.path.join(ROOT, "稿", "参考文献_引文清单.md")

pool = {}
for f in ["参考文献_原始.json", "参考文献_精查.json", "参考文献_扩充.json"]:
    j = json.load(open(os.path.join(ROOT, "稿", f), encoding="utf-8"))
    if isinstance(j, list):
        items = j
    else:
        items = j.get("landmarks", [])
        for v in j.get("topics", {}).values():
            items += v
    for it in items:
        if it.get("pmid"):
            pool.setdefault(str(it["pmid"]), it)

text = open(MD, encoding="utf-8").read()
pmids = re.findall(r"\[REF:(\d+)\]", text)
order = []
for p in pmids:
    if p not in order:
        order.append(p)

missing = [p for p in order if p not in pool]
if missing:
    print("警告：以下 PMID 不在本地文献池，尝试在线补取")
    API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    for p in missing:
        try:
            url = f"{API}?query=EXT_ID:{p}%20AND%20(SRC:MED)&format=json&resultType=core"
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=40) as r:
                res = json.load(r).get("resultList", {}).get("result", [])
            if res:
                a = res[0]
                auth = a.get("authorString", "").rstrip(".")
                parts = [x.strip() for x in auth.split(",") if x.strip()]
                if len(parts) > 3: auth = ", ".join(parts[:3]) + ", et al"
                ji = a.get("journalInfo", {})
                j = ji.get("journal", {}).get("title", ""); vol = ji.get("volume", ""); pages = a.get("pageInfo", "")
                pool[p] = {"pmid": p, "doi": a.get("doi", ""), "title": a.get("title", "").strip(),
                           "journal": j, "year": a.get("pubYear", ""), "volume": vol, "pages": pages, "authors": auth,
                           "ref": f"{auth}. {a.get('title','').strip().rstrip('.')}. {j}, {a.get('pubYear','')}"
                                  + (f", {vol}" if vol else "") + (f": {pages}" if pages else "")
                                  + (f". doi: {a.get('doi')}" if a.get("doi") else "") + f". PMID: {p}"}
                print("  补取成功:", p)
            else:
                print("  仍缺失:", p)
        except Exception as e:
            print("  补取失败:", p, e)

num = {p: i + 1 for i, p in enumerate(order)}
def repl(m):
    p = m.group(1)
    return f"[{num[p]}]" if p in num else m.group(0)
new = re.sub(r"\[REF:(\d+)\]", repl, text)

lines = ["# 参考文献", ""]
for p, i in sorted(num.items(), key=lambda kv: kv[1]):
    it = pool.get(p)
    if not it:
        lines.append(f"{i}. （未取到元数据，PMID {p}）")
        continue
    lines.append(f"{i}. {it['ref']}")
refmd = "\n".join(lines) + "\n"

new = new.replace("<!--REFS-->", refmd)
open(OUT_MD, "w", encoding="utf-8").write(new)
open(OUT_REF, "w", encoding="utf-8").write(refmd)

print("引用文献总数:", len(num))
print("正文剩余未解析占位符:", len(re.findall(r"\[REF:", new)))
print("\n--- 前 8 条 ---")
print("\n".join(lines[2:10]))
