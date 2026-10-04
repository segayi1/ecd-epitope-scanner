"""Render and assemble the PD-1 antibody epitope atlas report."""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt, RGBColor

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "pd1_antibody_epitope_atlas"
DATA = json.loads((OUT / "interface_data.json").read_text(encoding="utf-8"))
FIGURES = OUT / "figures"
FIGURES.mkdir(exist_ok=True)
PYMOL = Path(r"F:\pymol\Scripts\pymol.exe")

EPITOPE_ZONES = [
    ("N端环区", "N-loop/A′/A", "25–43", "PD-1 N 端延伸段及 A′/A 邻区；NBO1a 对该区域有直接接触。"),
    ("前侧β片层", "CC′/C′D", "58–84", "PD-L1结合面周围的前侧片层和环区；多种阻断抗体集中于此。"),
    ("FG环热点", "FG loop", "123–142", "PD-1 FG环及相邻结构；常与前侧β片层构成连续复合表位。"),
    ("膜近端背面", "A/F/G 背侧面", "结构定义，不以单一连续编号段概括", "Clone 19 的背侧、膜近端结合方向；不同于典型配体阻断抗体。"),
    ("侧向 C′D/D–E 面", "C′D + D/E strands", "约 80–115（文献映射）", "UDIZ-007/008 的XL-MS映射表位；无共晶，单列为非晶体证据。"),
]


def pymol_quote(path: Path) -> str:
    return str(path.resolve()).replace("\\", "/").replace('"', '\\"')


def image_filename(antibody):
    return f"{antibody['pdb']}_{antibody['name'].replace('/', '-').replace(' ', '-')}.png"


def render_all():
    if not PYMOL.exists():
        raise FileNotFoundError(f"PyMOL not found: {PYMOL}")
    commands = [
        "reinitialize",
        "bg_color white",
        "set antialias, 2",
        "set ray_shadows, 0",
        "set specular, 0.25",
        "set shininess, 18",
        "set surface_quality, 1",
        "set ambient, 0.28",
        "set ray_opaque_background, 1",
        "set cartoon_smooth_loops, 1",
        "set valence, 0",
        "set ray_trace_mode, 0",
        "set ray_trace_color, grey70",
    ]
    for antibody in DATA["antibodies"]:
        s = antibody.get("structure")
        if not s:
            continue
        pdb = antibody["pdb"]
        pdb_path = OUT / "structures" / f"{pdb}_aligned.pdb"
        img_path = FIGURES / image_filename(antibody)
        obj = f"c_{pdb}"
        pd1_sel = f"{obj} and chain {s['antigen_chain']}"
        fab_sel = " or ".join(f"({obj} and chain {c})" for c in s["binder_chains"])
        commands.extend([
            f"load \"{pymol_quote(pdb_path)}\", {obj}",
            f"hide everything, {obj}",
            f"show surface, {pd1_sel}",
            f"color 0xB9D1D8, {pd1_sel}",
            f"show surface, ({fab_sel})",
            f"color 0x426E9A, ({fab_sel})",
            f"set transparency, 0.02, {pd1_sel}",
            f"set transparency, 0.0, ({fab_sel})",
        ])
        for pos, distance in s["distance_by_residue_A"].items():
            if distance <= 3.5:
                color = "0xB6233C"
            elif distance <= 4.5:
                color = "0xD64C50"
            elif distance <= 6.0:
                color = "0xEE9364"
            else:
                color = "0xF5D896"
            sel = f"{pd1_sel} and resi {pos}"
            commands.append(f"color {color}, {sel}")
        commands.extend([
            f"orient {pd1_sel}",
            "zoom all, 8",
            f"center {pd1_sel}",
            f"png figures/{img_path.name}, 1200, 900, ray=1",
            f"delete {obj}",
        ])
    commands.append("quit")
    pml = OUT / "render_all.pml"
    pml.write_text("\n".join(commands) + "\n", encoding="utf-8")
    result = subprocess.run([str(PYMOL), "-cq", "-d", "@render_all.pml"], cwd=OUT,
                            capture_output=True, text=True, timeout=900)
    if result.returncode != 0:
        raise RuntimeError(f"PyMOL render failed ({result.returncode}):\n{result.stdout}\n{result.stderr}")
    expected = [FIGURES / image_filename(a)
                for a in DATA["antibodies"] if a.get("structure")]
    missing = [str(p) for p in expected if not p.exists() or p.stat().st_size < 20_000]
    if missing:
        raise RuntimeError("Missing or undersized PyMOL image(s): " + ", ".join(missing))
    return expected


def fmt_residues(values):
    vals = sorted(set(values))
    if not vals:
        return "—"
    runs = []
    start = prev = vals[0]
    for n in vals[1:]:
        if n == prev + 1:
            prev = n
        else:
            runs.append(f"{start}–{prev}" if start != prev else str(start))
            start = prev = n
    runs.append(f"{start}–{prev}" if start != prev else str(start))
    return ", ".join(runs)


def hyperlink(paragraph, text, url):
    part = paragraph.part
    rid = part.relate_to(url, "http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink", is_external=True)
    hyperlink = __import__("docx").oxml.shared.OxmlElement("w:hyperlink")
    hyperlink.set(__import__("docx").oxml.shared.qn("r:id"), rid)
    run = __import__("docx").oxml.shared.OxmlElement("w:r")
    props = __import__("docx").oxml.shared.OxmlElement("w:rPr")
    color = __import__("docx").oxml.shared.OxmlElement("w:color")
    color.set(__import__("docx").oxml.shared.qn("w:val"), "1565C0")
    props.append(color)
    run.append(props)
    text_node = __import__("docx").oxml.shared.OxmlElement("w:t")
    text_node.text = text
    run.append(text_node)
    hyperlink.append(run)
    paragraph._p.append(hyperlink)


def set_cell_shading(cell, fill):
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = __import__("docx").oxml.shared.OxmlElement("w:shd")
    shd.set(__import__("docx").oxml.shared.qn("w:fill"), fill)
    tc_pr.append(shd)


def set_repeat_table_header(row):
    tr_pr = row._tr.get_or_add_trPr()
    tbl_header = __import__("docx").oxml.shared.OxmlElement("w:tblHeader")
    tbl_header.set(__import__("docx").oxml.shared.qn("w:val"), "true")
    tr_pr.append(tbl_header)


def build_docx(images):
    doc = Document()
    sec = doc.sections[0]
    sec.top_margin = Inches(0.58)
    sec.bottom_margin = Inches(0.58)
    sec.left_margin = Inches(0.62)
    sec.right_margin = Inches(0.62)
    normal = doc.styles["Normal"]
    normal.font.name = "Microsoft YaHei"
    normal.font.size = Pt(9)
    normal.font.color.rgb = RGBColor(42, 56, 65)
    normal.paragraph_format.space_after = Pt(4)
    for style_name in ["Title", "Heading 1", "Heading 2", "Heading 3"]:
        doc.styles[style_name].font.name = "Microsoft YaHei"
        doc.styles[style_name].font.color.rgb = RGBColor(31, 68, 91)
    doc.styles["Title"].font.size = Pt(27)
    doc.styles["Heading 1"].font.size = Pt(17)
    doc.styles["Heading 2"].font.size = Pt(12)

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("PD‑1 抗体结合表位图谱")
    r.bold = True
    r.font.size = Pt(27)
    r.font.color.rgb = RGBColor(31, 68, 91)
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run("已检索抗体来源、共晶界面与表位分区汇总\n检索核对：2026‑09‑26 · 结构参考：人 PD‑1 UniProt Q15116")
    p.paragraph_format.space_after = Pt(13)

    doc.add_heading("结论速览", level=1)
    for text in [
        "本次归档 18 个抗体条目：16 个有 PD‑1–抗体实验复合物结构；普樂利單抗（prolgolimab）有文献中的计算模型而未找到共晶；UDIZ‑007/008 有交联质谱/表位映射而未找到共晶。",
        "大部分有共晶的阻断抗体接触 PD‑1 前侧 CC′/C′D β片层并延伸至 FG 环；纳武利尤单抗、D12、NBO1a 等更显著涉及 N 端；Clone 19 位于膜近端背面；UDIZ‑007/008 被报道为侧向 C′D/D–E 区域。",
        "结构图统一叠合到 5WT9 的 PD‑1 坐标系。PD‑1 以浅青灰表面显示；Fab 用较深蓝色；PD‑1 表位按距离着色：外围 4.5–8 Å 为浅黄，3.5–4.5 Å 为红橙，≤3.5 Å 为深红。每个小图只表示一个实验复合物，不暗示多抗体可同时结合。",
    ]:
        doc.add_paragraph(text, style="List Bullet")

    doc.add_heading("表位分区", level=1)
    doc.add_paragraph("按 PD‑1 UniProt 前体序号列出主要结构表面区域。连续结构面可跨越不连续序列；实际抗体界面以表 2 的晶体接触残基为准。一个抗体可横跨多个分区，因此分区用于归纳空间邻近表位，而非互斥分类。")
    zone_table = doc.add_table(rows=1, cols=4)
    zone_table.style = "Light Shading Accent 1"
    zone_table.alignment = WD_TABLE_ALIGNMENT.CENTER
    for cell, title in zip(zone_table.rows[0].cells, ["分区名称", "结构区域", "PD‑1 序号", "空间含义"]):
        cell.text = title
    set_repeat_table_header(zone_table.rows[0])
    for name, region, pos, desc in EPITOPE_ZONES:
        cells = zone_table.add_row().cells
        for c, val in zip(cells, [name, region, pos, desc]):
            c.text = val
            c.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER

    groups = {
        "N端环区 + FG环复合表位": ["Nivolumab", "GY-5", "NBO1a", "D12"],
        "前侧β片层（CC′/C′D）": ["Pembrolizumab", "Tislelizumab", "Sasanlimab"],
        "BC/C′D–FG环连续脊": ["Toripalimab", "Camrelizumab", "Cemiplimab", "Dostarlimab", "MW11-h317", "mAb059c", "Serplulimab", "GY-14"],
        "膜近端背面": ["Clone 19"],
        "侧向 C′D/D–E 面（非晶体映射）": ["UDIZ-007 / UDIZ-008"],
        "FG/配体面重叠模型（非晶体）": ["Prolgolimab"],
    }
    doc.add_heading("空间邻近的抗体归组", level=2)
    group_table = doc.add_table(rows=1, cols=2)
    group_table.style = "Light Shading Accent 1"
    group_table.rows[0].cells[0].text = "归纳名称"
    group_table.rows[0].cells[1].text = "对应抗体"
    set_repeat_table_header(group_table.rows[0])
    for name, members in groups.items():
        row = group_table.add_row().cells
        row[0].text, row[1].text = name, "、".join(members)

    doc.add_heading("抗体来源与晶体界面对应表", level=1)
    doc.add_paragraph("PD‑1 界面残基通过各复合物结构的原子间距离计算：直接接触 ≤4.5 Å；近邻范围 ≤8 Å。表中列出直接接触残基（UniProt Q15116 前体编号）。晶体结构中的 author chain 与报告图中 PD‑1/Fab 链均以 RCSB 坐标文件为准。")
    # landscape table on its own section
    sec2 = doc.add_section()
    sec2.orientation = WD_ORIENT.LANDSCAPE
    sec2.page_width, sec2.page_height = Inches(11.7), Inches(8.3)
    sec2.top_margin = Inches(0.5)
    sec2.bottom_margin = Inches(0.5)
    sec2.left_margin = Inches(0.45)
    sec2.right_margin = Inches(0.45)
    table = doc.add_table(rows=1, cols=7)
    table.style = "Table Grid"
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    headers = ["抗体", "来源/开发方", "抗体来源类型", "PDB / 分辨率", "PD‑1直接接触残基 ≤4.5 Å", "表位归组", "证据"]
    for c, h in zip(table.rows[0].cells, headers):
        c.text = h
        set_cell_shading(c, "DCEAF0")
    set_repeat_table_header(table.rows[0])
    for a in DATA["antibodies"]:
        s = a.get("structure")
        contacts = fmt_residues(s["contacts_4_5A"]) if s else "未解析"
        res = a.get("resolution_A")
        pdb_text = f"{a['pdb']} · {res:.2f} Å" if a.get("pdb") and res else (a["pdb"] or "无共晶")
        evidence = "X-ray 共晶" if s else a.get("evidence", "非晶体映射")
        row = table.add_row().cells
        values = [a["name"] + (f"\n{a['zh']}" if a["zh"] != a["name"] else ""), a["source"], a["origin"], pdb_text, contacts, a["group"], evidence]
        for c, val in zip(row, values):
            c.text = val
            c.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            for p in c.paragraphs:
                for run in p.runs:
                    run.font.size = Pt(7.3)

    doc.add_page_break()
    doc.add_heading("PD‑1–抗体共晶结构图", level=1)
    doc.add_paragraph("每幅图为 PDB 坐标中 PD‑1 与其对应 Fab 的结构表面；所有 PD‑1 已叠合到 5WT9 参考结构。图中着色残基按各结构中到抗体任一原子的最近重原子距离分档。")
    key = doc.add_paragraph()
    key.alignment = WD_ALIGN_PARAGRAPH.CENTER
    key.add_run("PD‑1：浅青灰     Fab：深蓝     表位外围 4.5–8 Å：浅黄     直接接触：红橙至深红")

    by_pdb = {a["pdb"]: a for a in DATA["antibodies"] if a.get("structure")}
    for n, a in enumerate(by_pdb.values()):
        if n % 4 == 0:
            if n:
                doc.add_page_break()
            plate = doc.add_table(rows=2, cols=2)
            plate.alignment = WD_TABLE_ALIGNMENT.CENTER
            plate.autofit = False
            plate.columns[0].width = Inches(5.05)
            plate.columns[1].width = Inches(5.05)
        cell = plate.cell((n % 4) // 2, n % 2)
        cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.TOP
        p = cell.paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        run = p.add_run(f"{a['name']}  ·  {a['pdb']}\n")
        run.bold = True
        run.font.size = Pt(10)
        s = a["structure"]
        run = p.add_run(f"PD‑1链 {s['antigen_chain']} | Fab链 {','.join(s['binder_chains'])} | ≤4.5 Å：{fmt_residues(s['contacts_4_5A'])}")
        run.font.size = Pt(7.3)
        p = cell.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        image = FIGURES / image_filename(a)
        p.add_run().add_picture(str(image), width=Inches(4.75))
        p = cell.add_paragraph(f"来源：{a['source']} · {a['group']}")
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        for run in p.runs:
            run.font.size = Pt(7.3)
        if n % 2 == 1:
            plate.rows[(n % 4) // 2].height = Inches(3.2)

    doc.add_page_break()
    doc.add_heading("非共晶条目与证据边界", level=1)
    doc.add_paragraph("普樂利單抗（prolgolimab）：已发表比较研究报告了与 pembrolizumab/nivolumab 对照的计算结构模型和预测表位；该界面属于模型推断，不是抗体–PD‑1 共晶接触列表。")
    doc.add_paragraph("UDIZ‑007/UDIZ‑008：来自同一 V 区的工程抗体构型；研究以交联质谱等实验定位其独特表位，提示靠近 C′D、D/E 链区且位于 PD‑L1 结合面相反侧；未发现可用于逐原子界面计算的 PD‑1 共晶。")
    doc.add_paragraph("本报告的“直接接触”是几何距离定义，不等同于文献逐项验证的氢键/盐桥，也不表示所有接触残基贡献相同结合能。晶体结构中的缺失环区、糖基化状态、构建体边界及晶格接触会影响可见界面；彩图为去除水分子和非蛋白配体后的蛋白表面示意。")

    doc.add_heading("结构与文献来源", level=1)
    for a in DATA["antibodies"]:
        p = doc.add_paragraph(style="List Bullet")
        if a.get("pdb"):
            hyperlink(p, f"{a['name']} · RCSB PDB {a['pdb']}", f"https://www.rcsb.org/structure/{a['pdb']}")
            citation = a.get("citation_title") or a.get("title") or "RCSB PDB entry"
            tail = f" — {citation}"
            if a.get("citation_year"):
                tail += f" ({a['citation_year']})"
            if a.get("doi"):
                tail += f" · DOI: {a['doi']}"
            p.add_run(tail)
        else:
            if a.get("source_url"):
                hyperlink(p, f"{a['name']} · {a['evidence']}", a["source_url"])
            else:
                p.add_run(f"{a['name']} · {a['evidence']}")
    p = doc.add_paragraph(style="List Bullet")
    hyperlink(p, "Structural characterization review (2023)", "https://pmc.ncbi.nlm.nih.gov/articles/PMC10399482/")
    p.add_run(" — PD‑1 抗体复合物结构和表位比较。")
    p = doc.add_paragraph(style="List Bullet")
    hyperlink(p, "PD‑1 抗体结构区域归类与 UDIZ 表位研究 (2025)", "https://pmc.ncbi.nlm.nih.gov/articles/PMC12667665/")
    p.add_run(" — 结构区域分类与 UDIZ‑007/008 实验映射。")
    p = doc.add_paragraph(style="List Bullet")
    hyperlink(p, "D12 发现与复合物结构 (2022)", "https://pmc.ncbi.nlm.nih.gov/articles/PMC9667898/")
    p.add_run("；")
    hyperlink(p, "Clone 19 结构 (PDB 9HK1)", "https://www.rcsb.org/structure/9HK1")
    p.add_run("。")

    docx_path = OUT / "PD1抗体名称来源与结合表位图谱.docx"
    doc.save(docx_path)
    return docx_path


def build_markdown():
    lines = [
        "# PD‑1 抗体名称、来源与结合表位图谱", "",
        "检索核对：2026‑09‑26。结构参考：人 PD‑1 UniProt Q15116；16 个 PD‑1–抗体实验复合物，另含 2 个非共晶证据条目。", "",
        "## 表位分区", "",
        "| 分区 | 区域/残基 | 抗体 |", "|---|---|---|",
    ]
    groups = {
        "N端环区 + FG环复合表位": ["Nivolumab", "GY-5", "NBO1a", "D12"],
        "前侧β片层（CC′/C′D）": ["Pembrolizumab", "Tislelizumab", "Sasanlimab"],
        "BC/C′D–FG环连续脊": ["Toripalimab", "Camrelizumab", "Cemiplimab", "Dostarlimab", "MW11-h317", "mAb059c", "Serplulimab", "GY-14"],
        "膜近端背面": ["Clone 19"],
        "侧向 C′D/D–E 面（非晶体映射）": ["UDIZ-007 / UDIZ-008"],
        "FG/配体面重叠模型（非晶体）": ["Prolgolimab"],
    }
    for name, members in groups.items():
        desc = next((f"{region} {pos}" for title, region, pos, _ in EPITOPE_ZONES if title in name or name.startswith(title)), "跨区/研究定义")
        lines.append(f"| {name} | {desc} | {', '.join(members)} |")
    lines += ["", "## 抗体、来源与界面", "", "接触定义：重原子最短距离 ≤4.5 Å（直接接触）；PD‑1 UniProt Q15116 前体编号。图中表位颜色：外围 4.5–8 Å 浅黄；近程红橙；≤3.5 Å 深红。", "",
              "| 抗体 | 来源/开发方 | 来源类型 | PDB/分辨率 | 直接接触 PD‑1 残基 | 证据与来源 |", "|---|---|---|---|---|---|"]
    for a in DATA["antibodies"]:
        s = a.get("structure")
        contact = fmt_residues(s["contacts_4_5A"]) if s else "—"
        if a.get("pdb"):
            res = f"{a['resolution_A']:.2f} Å" if a.get("resolution_A") else "—"
            pdb = f"[{a['pdb']}](https://www.rcsb.org/structure/{a['pdb']}) ({res})"
            source = f"{pdb}; [文献](https://doi.org/{a['doi']})" if a.get("doi") else pdb
            ev = "X-ray 共晶"
        else:
            pdb, ev = "无共晶", a.get("evidence", "非晶体映射")
            source = f"[{ev}]({a['source_url']})" if a.get("source_url") else ev
        lines.append(f"| {a['name']} ({a['zh']}) | {a['source']} | {a['origin']} | {pdb} | {contact} | {ev}; {source} |")
    lines += ["", "## 共晶结构图", ""]
    for a in DATA["antibodies"]:
        if not a.get("structure"):
            continue
        img = f"figures/{image_filename(a)}"
        lines += [f"### {a['name']} · PDB {a['pdb']}", "", f"![{a['name']} 与 PD‑1 共晶结构]({img})", "",
                  f"PD‑1链 {a['structure']['antigen_chain']}；Fab链 {','.join(a['structure']['binder_chains'])}；直接接触位点：{fmt_residues(a['structure']['contacts_4_5A'])}。", ""]
    lines += ["## 方法与证据限制", "", "复合物按 PD‑1 CA 原子序列匹配到 5WT9 链 G 并结构叠合；界面以任一抗体重链/轻链重原子到 PD‑1 残基的最短距离计算。区域名称按 PD‑1 结构拓扑归纳，多个抗体可横跨或部分共享区域。彩色表面图已去掉水及非蛋白配体。", "",
              "Prolgolimab 的界面来自计算模型；UDIZ‑007/008 来自交联质谱等映射，均不应解读成逐原子共晶界面。", "",
              "参考综述：[PD‑1/PD‑L1 抗体复合物结构综述](https://pmc.ncbi.nlm.nih.gov/articles/PMC10399482/)；[2025 年结构分区与 UDIZ 研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC12667665/)；[D12 结构研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC9667898/)。"]
    (OUT / "PD1抗体名称来源与结合表位图谱.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    expected = [FIGURES / image_filename(a) for a in DATA["antibodies"] if a.get("structure")]
    image_paths = expected if all(p.exists() and p.stat().st_size > 20_000 for p in expected) else render_all()
    docx = build_docx(image_paths)
    build_markdown()
    print(f"Rendered {len(image_paths)} complex figures; wrote {docx}")
