"""v2 图件拼版：主图 6 至 8 + 补充图 S1 至 S4（矢量 PDF + PNG 预览）"""
import os

import fitz

ROOT = "/mnt/f/0000_Bioinfo/222_tcgaSHAP_胶质瘤"
FIG = os.path.join(ROOT, "03_figures")
OUT = os.path.join(ROOT, "稿", "figures_v2")
os.makedirs(OUT, exist_ok=True)

LAYOUTS = {
    "图6_亚组与通路": dict(panels=["Fig6A_subgroup_forest.pdf", "Fig6B_km_strata.pdf",
                               "Fig6C_hallmark_heatmap.pdf", "Fig6D_hallmark_violin.pdf"],
                       cols=2, rows=2, page=(11.6, 10.6)),
    "图7_免疫景观": dict(panels=["Fig6A_ssgsea_immune.pdf", "Fig6B_ssgsea_immune2.pdf",
                             "Fig6C_imm_cells.pdf", "Fig8A_immSubtype.pdf",
                             "Fig6E_mcpcounter.pdf"],
                      cols=2, rows=3, page=(11.6, 11.6)),
    "图8_TMB与药敏": dict(panels=["Fig8B_tmb.pdf", "Fig8C_tmb_surv.pdf",
                                "Fig11B_drug_volcano.pdf", "Fig12_dca.pdf"],
                       cols=2, rows=2, page=(10.4, 8.4)),
    "图9_单细胞": dict(panels=["Fig9A_scRNA_umap.pdf", "Fig9B_scRNA_model_genes.pdf",
                             "Fig9C_scRNA_risk_module.pdf", "Fig9D_scRNA_correlation.pdf"],
                     cols=2, rows=2, page=(10.4, 8.6)),
    "图S1_SHAP个体解释": dict(panels=["Fig2C_shap_waterfall.pdf", "Fig2D_shap_force.pdf"],
                          cols=1, rows=2, page=(8.0, 8.0)),
    "图S2_TIDE免疫治疗预测": dict(panels=["Fig7A_tide.pdf", "Fig7B_tide_response.pdf"],
                            cols=1, rows=2, page=(9.2, 8.6)),
    "图S3_MR与共定位": dict(panels=["Fig10A_mr_forest.pdf", "Fig10B_mr_coloc.pdf"],
                         cols=2, rows=1, page=(12.0, 5.2)),
    "图S4_药敏箱线": dict(panels=["Fig11A_drug_box.pdf"], cols=1, rows=1, page=(11.0, 5.6)),
}


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


for name, spec in LAYOUTS.items():
    pw, ph = spec["page"]; cols, rows = spec["cols"], spec["rows"]
    doc = fitz.open(); page = doc.new_page(width=pw * 72, height=ph * 72)
    cw, ch = pw / cols, ph / rows
    idx = 0
    for i, pf in enumerate(spec["panels"]):
        p = os.path.join(FIG, pf)
        if not os.path.exists(p):
            print("  [缺]", pf); continue
        src = fitz.open(p)
        r, c = divmod(i, cols)
        clip = content_clip(src[0])
        sw, sh = clip.width / 72, clip.height / 72
        avail_w, avail_h = cw - 0.5, ch - 0.62
        k = min(avail_w / sw, avail_h / sh)
        w, h = sw * k, sh * k
        x0 = (c * cw + (cw - w) / 2) * 72
        y0 = (r * ch + 0.52 + (avail_h - h) / 2) * 72
        page.show_pdf_page(fitz.Rect(x0, y0, x0 + w * 72, y0 + h * 72), src, 0, clip=clip)
        page.insert_text(fitz.Point(c * cw * 72 + 16, (r * ch + 0.36) * 72),
                         chr(ord("A") + idx), fontsize=14, fontname="hebo")
        idx += 1
        src.close()
    pdf = os.path.join(OUT, name + ".pdf")
    doc.save(pdf)
    page.get_pixmap(dpi=160).save(os.path.join(OUT, name + ".png"))
    doc.close()
    print("saved:", name, round(os.path.getsize(pdf) / 1024), "KB")
print("done ->", OUT)
