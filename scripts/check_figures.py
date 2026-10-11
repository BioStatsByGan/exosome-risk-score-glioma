"""出图后验墨迹：检查单张 panel 与拼版主图是否有空白页 / 空白 panel。
用法：python3 scripts/check_figures.py            # 检查 03_figures 全部 pdf
      python3 scripts/check_figures.py 稿/figures # 检查指定目录
"""
import glob, os, sys

import fitz
import numpy as np

ROOT = "/mnt/f/0000_Bioinfo/222_tcgaSHAP_胶质瘤"


def ink(page, clip=None):
    pix = page.get_pixmap(dpi=60, clip=clip)
    a = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)
    return (a[:, :, :3].min(axis=2) < 230).mean() * 100


def check_dir(d, grid=None):
    bad = []
    for f in sorted(glob.glob(os.path.join(d, "*.pdf"))):
        doc = fitz.open(f)
        pg = doc[0]
        text = len(pg.get_text().strip())
        v = ink(pg)
        tag = ""
        if v < 0.3 and text == 0:
            tag = "  <== 空白页!"
            bad.append(os.path.basename(f))
        print(f"{os.path.basename(f):40s} 墨迹 {v:5.2f}%  文字 {text:5d}{tag}")
        if grid:
            cols, rows = grid
            r = pg.rect
            for i in range(cols * rows):
                rr, c = divmod(i, cols)
                clip = fitz.Rect(c * r.width / cols, rr * r.height / rows,
                                 (c + 1) * r.width / cols, (rr + 1) * r.height / rows)
                pv = ink(pg, clip)
                if pv < 0.3:
                    print(f"      panel {chr(65 + i)} 墨迹 {pv:5.2f}%  <== 空白 panel!")
                    bad.append(f"{os.path.basename(f)}#{chr(65 + i)}")
    return bad


dirs = sys.argv[1:] or ["03_figures"]
grids = {"稿/figures": (2, 3)}
bad = []
for d in dirs:
    p = d if os.path.isabs(d) else os.path.join(ROOT, d)
    print(f"===== {p} =====")
    bad += check_dir(p, grids.get(d))
print()
print("结论:", "全部有内容" if not bad else f"发现 {len(bad)} 处空白 -> {bad}")
