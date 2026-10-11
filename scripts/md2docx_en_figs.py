#!/usr/bin/env python
"""Markdown -> DOCX（英文投稿版）：Times New Roman 12 pt / 全黑 / 双倍行距 / 表格线 / 图片居中 15.5 cm
用法：md2docx_en_figs.py 正文.md 正文.docx
"""
import os, re, sys
from docx import Document
from docx.shared import Pt, RGBColor, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn

src, dst = sys.argv[1], sys.argv[2]
base = os.path.dirname(os.path.abspath(src))
lines = open(src, encoding="utf-8").read().split("\n")
doc = Document()


def set_font(style, size, bold=False):
    f = style.font
    f.name = "Times New Roman"
    f.size = Pt(size)
    f.bold = bold
    f.color.rgb = RGBColor(0, 0, 0)
    rpr = style.element.get_or_add_rPr()
    rf = rpr.find(qn("w:rFonts"))
    if rf is None:
        rf = rpr.makeelement(qn("w:rFonts"), {})
        rpr.append(rf)
    for k in ("w:ascii", "w:hAnsi", "w:cs", "w:eastAsia"):
        rf.set(qn(k), "Times New Roman")


st = doc.styles["Normal"]
set_font(st, 12)
st.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
st.paragraph_format.space_after = Pt(0)
for s, sz in (("Heading 1", 15), ("Heading 2", 13), ("Heading 3", 12)):
    set_font(doc.styles[s], sz, bold=True)
    doc.styles[s].paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE


def add_runs(par, text):
    for seg in re.split(r"(\*\*[^*]+\*\*|\*[^*]+\*)", text):
        if not seg:
            continue
        if seg.startswith("**") and seg.endswith("**"):
            par.add_run(seg[2:-2]).bold = True
        elif seg.startswith("*") and seg.endswith("*") and len(seg) > 2:
            par.add_run(seg[1:-1]).italic = True
        else:
            par.add_run(seg)


def add_table(rows):
    ncol = max(len(r) for r in rows)
    rows = [r + [""] * (ncol - len(r)) for r in rows]
    t = doc.add_table(rows=0, cols=ncol)
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for ri, row in enumerate(rows):
        cells = t.add_row().cells
        for ci, txt in enumerate(row):
            cells[ci].text = ""
            p = cells[ci].paragraphs[0]
            p.paragraph_format.space_after = Pt(0)
            p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            add_runs(p, txt.replace("<br>", " "))
            for r in p.runs:
                r.font.size = Pt(9.5)
                if ri == 0:
                    r.bold = True


i, in_fm, fm = 0, False, []
while i < len(lines):
    ln = lines[i].rstrip()
    if ln.strip() == "---" and i < 8:
        in_fm = not in_fm
        i += 1
        continue
    if in_fm:
        fm.append(ln)
        i += 1
        continue
    m = re.match(r"^!\[[^\]]*\]\(([^)]+)\)\s*$", ln.strip())
    if m:
        path = m.group(1)
        if not os.path.isabs(path):
            path = os.path.join(base, path)
        if os.path.exists(path):
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            p.paragraph_format.space_before = Pt(6)
            p.paragraph_format.space_after = Pt(2)
            p.add_run().add_picture(path, width=Cm(15.5))
        else:
            doc.add_paragraph(f"[missing image] {path}")
        i += 1
        continue
    if ln.startswith("|") and i + 1 < len(lines) and not lines[i + 1].replace("|", "").replace("-", "").replace(":", "").strip():
        rows = []
        while i < len(lines) and lines[i].startswith("|"):
            cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
            if not all(set(c) <= set("-: ") for c in cells):
                rows.append(cells)
            i += 1
        add_table(rows)
        doc.add_paragraph()
        continue
    if ln.startswith("### "):
        doc.add_heading(ln[4:], level=3)
    elif ln.startswith("## "):
        doc.add_heading(ln[3:], level=2)
    elif ln.startswith("# "):
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        r = p.add_run(ln[2:])
        r.bold = True
        r.font.size = Pt(16)
    elif re.match(r"^\s*[-*] ", ln):
        p = doc.add_paragraph(style="List Bullet")
        add_runs(p, re.sub(r"^\s*[-*] ", "", ln))
    elif ln.strip():
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        add_runs(p, ln.strip())
    i += 1

if fm:
    stamp = []
    p = doc.add_paragraph()
    r = p.add_run("FILE INFORMATION")
    r.bold = True
    r.font.color.rgb = RGBColor(0xA0, 0, 0)
    stamp.append(p)
    for f in fm:
        if f.strip():
            pp = doc.add_paragraph()
            rr = pp.add_run(f.strip())
            rr.font.size = Pt(9)
            rr.font.color.rgb = RGBColor(0xA0, 0, 0)
            stamp.append(pp)
    for el in reversed([x._element for x in stamp]):
        doc.element.body.insert(0, el)

doc.save(dst)
print("DOCX saved:", dst)
