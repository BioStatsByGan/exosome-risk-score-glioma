"""用 Crossref API 并发核验 DOI 是否存在（200=真实，404=查无此文），并比对期刊名。"""
import re, sys, json, urllib.request, urllib.parse, concurrent.futures as cf

refs = open("稿/参考文献_英文清单_v2.md", encoding="utf-8").read().splitlines()
items = []
for ln in refs:
    m = re.match(r"^(\d+)\.\s+(.*)$", ln.strip())
    if not m: continue
    body = m.group(2)
    dm = re.search(r"doi:(\S+?)\.?$", body)
    # 期刊名 = 倒数第二个句号段
    segs = [x.strip() for x in body.split(".") if x.strip()]
    jour = segs[-3] if len(segs) >= 3 else ""
    items.append({"n": int(m.group(1)), "doi": dm.group(1).rstrip(".") if dm else None, "jour": jour, "body": body})

def check(it):
    if not it["doi"]:
        return (it["n"], "NO_DOI", "", it["jour"])
    u = "https://api.crossref.org/works/" + urllib.parse.quote(it["doi"])
    try:
        req = urllib.request.Request(u, headers={"User-Agent": "refcheck/1.0 (mailto:chk@example.org)"})
        with urllib.request.urlopen(req, timeout=25) as r:
            d = json.load(r)["message"]
        ct = (d.get("container-title") or [""])[0]
        return (it["n"], "OK", ct, it["jour"])
    except urllib.error.HTTPError as e:
        return (it["n"], f"HTTP {e.code}", "", it["jour"])
    except Exception as e:
        return (it["n"], f"ERR {type(e).__name__}", "", it["jour"])

with cf.ThreadPoolExecutor(max_workers=8) as ex:
    res = list(ex.map(check, items))
bad = [r for r in res if r[1] != "OK"]
print(f"Crossref 核验 {len(res)} 条 → 存在 {len(res)-len(bad)}，异常 {len(bad)}")
for b in bad: print("  #%-3d %-12s doi-journal=%-38s | 稿中期刊=%s" % (b[0], b[1], b[2][:38], b[3][:34]))
json.dump(res, open("04_results/audit/ref_crossref_check.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
