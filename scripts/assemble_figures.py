"""把单张 panel 的矢量 PDF 拼成主图（保持矢量），并输出 PNG 预览。
用法：python3 scripts/assemble_figures.py
"""
import os, fitz

ROOT = "/mnt/f/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
FIG = os.path.join(ROOT, "03_figures")
OUT = os.path.join(ROOT, "稿", "figures")
os.makedirs(OUT, exist_ok=True)

LAYOUTS = {
    "图1_差异与富集": dict(panels=["Fig1A_volcano.pdf", "Fig1B_heatmap.pdf",
                                    "Fig1C_gsea_curves.pdf", "Fig1D_cnetplot.pdf"],
                        cols=2, rows=2, page=(7.4, 7.4)),
    "图2_模型构建与SHAP": dict(panels=["Fig2A_lasso_cv.pdf", "Fig2B_multicox_forest.pdf",
                                      "Fig2A_shap_bar.pdf", "Fig2B_shap_bee.pdf"],
                            cols=2, rows=2, page=(8.6, 7.4)),
    "图3_训练队列风险分层": dict(panels=["Fig4A_riskscore.pdf", "Fig4B_survstatus.pdf",
                                        "Fig4C_km_TCGA.pdf", "Fig4D_roc_TCGA.pdf"],
                              cols=2, rows=2, page=(8.0, 8.0)),
    "图5_独立预后与临床": dict(panels=["Fig5A_independence_TCGA.pdf", "Fig5A2_independence_multi_TCGA.pdf",
                                     "Fig5B_risk_by_clinical.pdf", "Fig5C_clinical_heatmap.pdf",
                                     "Fig5D_nomogram.pdf", "Fig5E_calibration.pdf"],
                            cols=2, rows=3, page=(8.6, 11.0)),
}
# 图4 外部验证本身就是 4 panel 拼好的
SINGLE = {"图4_外部验证": "Fig3_external_validation.pdf"}

def content_clip(page, pad=6):
    """取 pad=6 的源页内容外框（文字 + 矢量图元），用于裁掉四周空白，让 panel 填满格子。"""
    xs, ys, xe, ye = [], [], [], []
    for w in page.get_text("words"):
        xs.append(w[0]); ys.append(w[1]); xe.append(w[2]); ye.append(w[3])
    for d in page.get_drawings():
        r = d["rect"]
        if r.is_empty or r.width > page.rect.width * 0.98 and r.height > page.rect.height * 0.98:
            continue
        xs.append(r.x0); ys.append(r.y0); xe.append(r.x1); ye.append(r.y1)
    for im in page.get_image_info():
        r = fitz.Rect(im["bbox"])
        xs.append(r.x0); ys.append(r.y0); xe.append(r.x1); ye.append(r.y1)
    if not xs:
        return page.rect
    r = fitz.Rect(min(xs), min(ys), max(xe), max(ye))
    r = fitz.Rect(r.x0 - pad, r.y0 - pad, r.x1 + pad, r.y1 + pad)
    return r & page.rect


def compose(name, spec):
    pw, ph = spec["page"]
    cols, rows = spec["cols"], spec["rows"]
    doc = fitz.open()
    page = doc.new_page(width=pw * 72, height=ph * 72)
    cw, ch = pw / cols, ph / rows
    for i, pf in enumerate(spec["panels"]):
        p = os.path.join(FIG, pf)
        if not os.path.exists(p):
            print("  [缺] ", pf); continue
        src = fitz.open(p)
        r, c = divmod(i, cols)
        clip = content_clip(src[0])
        sw, sh = clip.width / 72, clip.height / 72
        avail_w, avail_h = cw - 0.45, ch - 0.58
        k = min(avail_w / sw, avail_h / sh)
        w, h = sw * k, sh * k
        x0 = (c * cw + (cw - w) / 2) * 72
        y0 = (r * ch + 0.48 + (avail_h - h) / 2) * 72
        page.show_pdf_page(fitz.Rect(x0, y0, x0 + w * 72, y0 + h * 72), src, 0, clip=clip)
        lab = chr(ord("A") + i)
        page.insert_text(fitz.Point(c * cw * 72 + 14, (r * ch + 0.34) * 72),
                         lab, fontsize=13, fontname="hebo")
        src.close()
    pdf = os.path.join(OUT, name + ".pdf")
    doc.save(pdf)
    pix = page.get_pixmap(dpi=150)
    pix.save(os.path.join(OUT, name + ".png"))
    doc.close()
    print("saved:", pdf, round(os.path.getsize(pdf) / 1024), "KB")

for nm, spec in LAYOUTS.items():
    compose(nm, spec)

for nm, f in SINGLE.items():
    src = fitz.open(os.path.join(FIG, f))
    dst = fitz.open(); dst.insert_pdf(src)
    pdf = os.path.join(OUT, nm + ".pdf"); dst.save(pdf)
    dst[0].get_pixmap(dpi=150).save(os.path.join(OUT, nm + ".png"))
    print("saved:", pdf, round(os.path.getsize(pdf) / 1024), "KB")

print("done ->", OUT)
