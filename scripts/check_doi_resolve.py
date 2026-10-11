"""并发用 doi.org 解析核验每条 DOI 是否真实（HEAD 请求，302=存在）。"""
import re, sys, urllib.request, concurrent.futures as cf

refs = open("稿/参考文献_英文清单_v2.md", encoding="utf-8").read().splitlines()
items = []
for ln in refs:
    m = re.match(r"^(\d+)\.\s+(.*)$", ln.strip())
    if not m: continue
    body = m.group(2)
    dm = re.search(r"doi:(\S+?)\.?$", body)
    items.append((int(m.group(1)), dm.group(1).rstrip(".") if dm else None, body[:80]))

def check(it):
    n, doi, _ = it
    if not doi: return (n, "NO_DOI", doi)
    url = "https://doi.org/" + doi
    for method in ("HEAD", "GET"):
        try:
            req = urllib.request.Request(url, method=method, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=25) as r:
                return (n, f"OK {r.status}", doi)
        except urllib.error.HTTPError as e:
            if e.code in (403, 406) and method == "HEAD":   # 部分出版社拒绝 HEAD，换 GET 再试
                continue
            return (n, f"HTTP {e.code}", doi)
        except Exception as e:
            if method == "HEAD": continue
            return (n, f"ERR {type(e).__name__}", doi)
    return (n, "UNKNOWN", doi)

with cf.ThreadPoolExecutor(max_workers=8) as ex:
    res = list(ex.map(check, items))
bad = [r for r in res if not r[1].startswith("OK")]
print(f"核验 {len(res)} 条，通过 {len(res)-len(bad)} 条")
for b in bad: print("  #%-3d %-22s %s" % b)
sys.stdout.flush()
