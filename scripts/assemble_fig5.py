"""图5 专用拼版：A、B 两个森林图各占整行（横向放大），C 至 F 两两并排。
用法：python3 scripts/assemble_fig5.py
"""
import os

import fitz

ROOT = "/mnt/f/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
FIG = os.path.join(ROOT, "03_figures")
OUT = os.path.join(ROOT, "稿", "figures")
os.makedirs(OUT, exist_ok=True)

PW, PH = 8.6, 12.4           # 英寸
ROWS = [                     # (行高权重, [panel 文件...]）
    (1.00, ["Fig5A_independence_TCGA.pdf"]),
    (1.05, ["Fig5A2_independence_multi_TCGA.pdf"]),
    (1.75, ["Fig5B_risk_by_clinical.pdf", "Fig5C_clinical_heatmap.pdf"]),
    (1.45, ["Fig5D_nomogram.pdf", "Fig5E_calibration.pdf"]),
]


def content_clip(page, pad=6):
    xs, ys, xe, ye = [], [], [], []
    for w in page.get_text("words"):
        xs.append(w[0]); ys.append(w[1]); xe.append(w[2]); ye.append(w[3])
    for d in page.get_drawings():
        r = d["rect"]
        if r.is_empty or (r.width > page.rect.width * 0.98 and r.height > page.rect.height * 0.98):
            continue
        xs.append(r.x0); ys.append(r.y0); xe.append(r.x1); ye.append(r.y1)
    for im in page.get_image_info():
        r = fitz.Rect(im["bbox"]); xs.append(r.x0); ys.append(r.y0); xe.append(r.x1); ye.append(r.y1)
    if not xs:
        return page.rect
    return (fitz.Rect(min(xs) - pad, min(ys) - pad, max(xe) + pad, max(ye) + pad)) & page.rect


doc = fitz.open()
page = doc.new_page(width=PW * 72, height=PH * 72)
tot = sum(w for w, _ in ROWS)
gap = 0.40
usable_h = PH - 0.5 - gap * (len(ROWS) - 1)
y = 0.5
idx = 0
for weight, panels in ROWS:
    rh = usable_h * weight / tot
    n = len(panels)
    for j, pf in enumerate(panels):
        src = fitz.open(os.path.join(FIG, pf))
        clip = content_clip(src[0])
        sw, sh = clip.width / 72, clip.height / 72
        cellw = (PW - 0.6) / n
        avail_w, avail_h = cellw - 0.25, rh - 0.30
        k = min(avail_w / sw, avail_h / sh)
        w, h = sw * k, sh * k
        x0 = (0.3 + j * cellw + (cellw - w) / 2) * 72
        y0 = (y + 0.26 + (avail_h - h) / 2) * 72
        page.show_pdf_page(fitz.Rect(x0, y0, x0 + w * 72, y0 + h * 72), src, 0, clip=clip)
        page.insert_text(fitz.Point((0.3 + j * cellw) * 72, (y + 0.18) * 72),
                         chr(ord("A") + idx), fontsize=14, fontname="hebo")
        idx += 1
        src.close()
    y += rh + gap

pdf = os.path.join(OUT, "图5_独立预后与临床.pdf")
doc.save(pdf)
pix = page.get_pixmap(dpi=170)
pix.save(os.path.join(OUT, "图5_独立预后与临床.png"))
doc.close()
print("saved:", pdf, round(os.path.getsize(pdf) / 1024), "KB")
