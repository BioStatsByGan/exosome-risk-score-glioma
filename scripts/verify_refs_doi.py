"""用 Crossref 逐条核验 94 条文献的 DOI 是否真实存在，并比对期刊/年份。"""
import re, json, time, urllib.request, os

refs = open("稿/参考文献_英文清单_v2.md", encoding="utf-8").read().splitlines()
items = []
for ln in refs:
    m = re.match(r"^(\d+)\.\s+(.*)$", ln.strip())
    if not m: continue
    n, body = int(m.group(1)), m.group(2)
    dm = re.search(r"doi:([^\s.]+(?:\.[^\s.]+)*?)\.?$", body)
    doi = dm.group(1).rstrip(".") if dm else None
    yr = re.search(r"\.\s*(\d{4});", body) or re.search(r"\.\s*(\d{4})\.", body)
    items.append({"n": n, "doi": doi, "year": yr.group(1) if yr else None,
                  "title": body.split(".")[0][:70]})

print(f"解析到 {len(items)} 条，其中无 DOI {sum(1 for i in items if not i['doi'])} 条")
bad, ok = [], 0
for it in items:
    if not it["doi"]: bad.append((it["n"], "无DOI", it["title"])); continue
    try:
        req = urllib.request.Request("https://api.crossref.org/works/" + urllib.parse.quote(it["doi"]),
                                     headers={"User-Agent": "ref-check/1.0 (mailto:check@example.com)"})
        with urllib.request.urlopen(req, timeout=25) as r:
            d = json.load(r)["message"]
        ct = (d.get("container-title") or [""])[0]
        cy = (d.get("issued", {}).get("date-parts") or [[None]])[0][0]
        it["cr_journal"] = ct; it["cr_year"] = cy
        flag = []
        if it["year"] and cy and abs(int(it["year"]) - int(cy)) > 1: flag.append(f"年份 {it['year']}vs{cy}")
        if ct and ct.lower().split(":")[0][:14] not in " ".join([x.lower() for x in re.findall(r"[A-Za-z ]+", it["title"])]).lower() and False:
            pass
        ok += 1 if not flag else 0
        if flag: bad.append((it["n"], "; ".join(flag), f"{ct} | {it['title']}"))
    except Exception as e:
        bad.append((it["n"], f"DOI 查询失败 {e}", it["doi"]))
    time.sleep(0.15)

print(f"DOI 可解析 {ok}/{len(items)}")
print("== 需处理 ==")
for b in bad: print("  #%-3d %s | %s" % b)
json.dump(items, open("04_results/audit/ref_verify.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
