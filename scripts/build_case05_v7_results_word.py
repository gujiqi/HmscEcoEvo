"""Build a Chinese Case05 v7 results atlas from the completed output folder."""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.style import WD_STYLE_TYPE
from docx.enum.table import WD_ALIGN_VERTICAL
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor


OUTPUT = Path(sys.argv[1]).resolve()
REPORT_DIR = OUTPUT / "06_report"
MAP_DIR = OUTPUT / "05_shared_scale_maps"
AUDIT_DIR = OUTPUT / "01_input_audit"
QUALITY_DIR = OUTPUT / "06_quality"

TIME_POINTS = [325, 300, 250, 200, 150, 100, 65, 20, 0]
METRICS = [
    "posterior_lineage_location_density",
    "posterior_lineage_location_density_log10",
    "posterior_weighted_ancestral_HMSC_linear_predictor",
    "Arias_destination_landscape_weight",
    "terminal_calibrated_lineage_support_intensity",
    "support_mixture_shannon",
    "support_mixture_effective_lineages",
]

METRIC_INFO = {
    "posterior_lineage_location_density": (
        "树条件谱系位置密度",
        "在板块搬运和球面扩散模型、以及现代终点记录条件下，所选活跃谱系在每个 4 度状态格中的相对位置密度总和。",
        "可比较不同时间和地点的相对空间支持；在 0 Ma 应与现代记录格的空间分布一致。",
        "不是出现概率、占据概率、丰度或真实物种丰富度。",
    ),
    "posterior_lineage_location_density_log10": (
        "树条件谱系位置密度 对数尺度",
        "与位置密度相同，但以 log10 尺度显示，便于观察弱而连续的长距离扩散尾部。",
        "适合判断是否存在连续扩散通道和低密度远距离位置支持。",
        "颜色差异表示对数密度差异，不可直接当作概率差或扩散速率。",
    ),
    "posterior_weighted_ancestral_HMSC_linear_predictor": (
        "祖先 HMSC 环境线性响应诊断",
        "由 HMSC posterior Beta 沿定年树重建祖先响应后，与对应时间古环境相乘得到的线性预测量。",
        "正值表示在该响应函数下环境支持较高，负值表示较低；用于解释环境过滤的空间格局。",
        "它没有进入位置扩散似然，不能被称为最终古分布或历史占据。",
    ),
    "Arias_destination_landscape_weight": (
        "Arias 风格目的地地貌权重",
        "由真实 PALEOMAP 陆地比例和地形通透性形成的目的地权重；海洋保留低但非零通行权重。",
        "用于识别扩散转移中较容易被到达的陆地地貌状态。",
        "不是植物生境质量、适宜度或局地存续概率。",
    ),
    "terminal_calibrated_lineage_support_intensity": (
        "扩散条件谱系支持强度",
        "将每个谱系的扩散位置密度按其 0 Ma 实际记录格数校准后，在一个格网中求和得到的支持强度。",
        "值高表示更多扩散条件谱系支持场在该处重叠；可比较空间重叠热点。",
        "不是历史物种数、alpha diversity、PD、FD 或占据概率。",
    ),
    "support_mixture_shannon": (
        "扩散条件谱系支持混合 Shannon 多样性",
        "在支持强度内将各谱系支持场归一化为混合权重后计算的 Shannon 指数。",
        "值高表示局地支持不是由单一谱系主导，而是由较均匀的多个谱系共同构成。",
        "不能解释为调查群落 Shannon 指数或完整历史植物群落多样性。",
    ),
    "support_mixture_effective_lineages": (
        "扩散条件有效谱系数",
        "support mixture Shannon 指数的指数变换 exp(H)，即等权混合时对应的有效支持谱系数量。",
        "便于以接近谱系数的直观尺度比较位置场混合程度。",
        "仍然是扩散条件支持诊断，不是已重建的古代物种丰富度。",
    ),
}


def set_cell_shading(cell, fill: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:fill"), fill)
    tc_pr.append(shd)


def set_cell_margins(cell, top=90, start=110, bottom=90, end=110) -> None:
    tc = cell._tc
    tc_pr = tc.get_or_add_tcPr()
    mar = tc_pr.first_child_found_in("w:tcMar")
    if mar is None:
        mar = OxmlElement("w:tcMar")
        tc_pr.append(mar)
    for side, value in (("top", top), ("start", start), ("bottom", bottom), ("end", end)):
        node = mar.find(qn(f"w:{side}"))
        if node is None:
            node = OxmlElement(f"w:{side}")
            mar.append(node)
        node.set(qn("w:w"), str(value))
        node.set(qn("w:type"), "dxa")


def set_repeat_table_header(row) -> None:
    tr_pr = row._tr.get_or_add_trPr()
    tbl_header = OxmlElement("w:tblHeader")
    tbl_header.set(qn("w:val"), "true")
    tr_pr.append(tbl_header)


def style_document(doc: Document) -> None:
    normal = doc.styles["Normal"]
    normal.font.name = "Aptos"
    normal._element.rPr.rFonts.set(qn("w:eastAsia"), "Microsoft YaHei")
    normal.font.size = Pt(9.5)
    normal.paragraph_format.space_after = Pt(5)
    normal.paragraph_format.line_spacing = 1.15
    for style_name, size in (("Title", 22), ("Heading 1", 15), ("Heading 2", 12), ("Heading 3", 10.5)):
        style = doc.styles[style_name]
        style.font.name = "Aptos Display" if style_name == "Title" else "Aptos"
        style._element.rPr.rFonts.set(qn("w:eastAsia"), "Microsoft YaHei")
        style.font.size = Pt(size)
        style.font.color.rgb = RGBColor(0, 0, 0)
        style.font.bold = True
        style.paragraph_format.space_before = Pt(12 if style_name != "Title" else 0)
        style.paragraph_format.space_after = Pt(6)
    caption = doc.styles["Caption"]
    caption.font.name = "Aptos"
    caption._element.rPr.rFonts.set(qn("w:eastAsia"), "Microsoft YaHei")
    caption.font.size = Pt(8.3)
    caption.font.color.rgb = RGBColor(0, 0, 0)
    caption.font.italic = True
    caption.paragraph_format.space_after = Pt(7)


def add_footer(section) -> None:
    p = section.footer.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("HmscEcoEvo Case05 v7  扩散条件谱系支持多样性结果教程")
    r.font.size = Pt(8)
    r.font.color.rgb = RGBColor(90, 90, 90)


def add_heading(doc: Document, text: str, level: int = 1) -> None:
    doc.add_heading(text, level=level)


def add_text(doc: Document, text: str, bold_prefix: str | None = None) -> None:
    p = doc.add_paragraph()
    if bold_prefix and text.startswith(bold_prefix):
        p.add_run(bold_prefix).bold = True
        p.add_run(text[len(bold_prefix):])
    else:
        p.add_run(text)


def load_csv(path: Path):
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def add_table(doc: Document, headers, rows, widths=None) -> None:
    table = doc.add_table(rows=1, cols=len(headers))
    table.style = "Table Grid"
    table.autofit = False
    hdr = table.rows[0]
    set_repeat_table_header(hdr)
    for i, header in enumerate(headers):
        cell = hdr.cells[i]
        cell.text = str(header)
        cell.vertical_alignment = WD_ALIGN_VERTICAL.CENTER
        set_cell_shading(cell, "1F4E78")
        set_cell_margins(cell)
        for run in cell.paragraphs[0].runs:
            run.font.bold = True
            run.font.color.rgb = RGBColor(255, 255, 255)
            run.font.size = Pt(8.5)
        if widths:
            cell.width = Cm(widths[i])
    for row_i, row in enumerate(rows):
        cells = table.add_row().cells
        for i, value in enumerate(row):
            cells[i].text = str(value)
            cells[i].vertical_alignment = WD_ALIGN_VERTICAL.CENTER
            set_cell_margins(cells[i])
            if row_i % 2:
                set_cell_shading(cells[i], "F2F6FA")
            for paragraph in cells[i].paragraphs:
                for run in paragraph.runs:
                    run.font.size = Pt(8.2)
            if widths:
                cells[i].width = Cm(widths[i])
    doc.add_paragraph()


def read_map_index():
    rows = load_csv(MAP_DIR / "case05_v7_map_index.csv")
    output = {}
    for row in rows:
        output[(row["metric"], int(float(row["time_ma"])))] = row
    return output


def add_metric_figure(doc: Document, row: dict, time_ma: int) -> None:
    metric = row["metric"]
    title, definition, reading, limitation = METRIC_INFO[metric]
    image_path = Path(row["png"])
    if not image_path.exists():
        add_text(doc, f"图件缺失：{image_path}")
        return
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run().add_picture(str(image_path), width=Cm(15.8))
    cap = doc.add_paragraph(style="Caption")
    cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
    cap.add_run(f"图 {title}，{time_ma} Ma。该指标在全部 66 个时间片使用固定色标 [{row['scale_min']}, {row['scale_max']}]。")
    add_text(doc, f"表示什么：{definition}", "表示什么：")
    add_text(doc, f"如何读：{reading}", "如何读：")
    add_text(doc, f"解释边界：{limitation}", "解释边界：")


def main() -> None:
    required = [MAP_DIR / "case05_v7_map_index.csv", QUALITY_DIR / "case05_v7_quality_gates.csv",
                AUDIT_DIR / "terminal_range_support_calibration.csv", AUDIT_DIR / "terminal_conditioning_audit.csv"]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("Missing required Case05 output files: " + "; ".join(missing))
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    out_docx = REPORT_DIR / "Case05_Plant200_扩散条件谱系支持多样性结果教程_zh.docx"
    index = read_map_index()
    quality = load_csv(QUALITY_DIR / "case05_v7_quality_gates.csv")
    terminal = load_csv(AUDIT_DIR / "terminal_conditioning_audit.csv")
    calibration = load_csv(AUDIT_DIR / "terminal_range_support_calibration.csv")

    doc = Document()
    section = doc.sections[0]
    section.top_margin = Cm(1.8)
    section.bottom_margin = Cm(1.7)
    section.left_margin = Cm(2.0)
    section.right_margin = Cm(2.0)
    style_document(doc)
    add_footer(section)

    title = doc.add_paragraph(style="Title")
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    title.add_run("Case05 Plant200 扩散条件谱系支持多样性结果教程")
    subtitle = doc.add_paragraph()
    subtitle.alignment = WD_ALIGN_PARAGRAPH.CENTER
    subtitle.add_run("全球 4 度古地理格网  325 Ma 至 0 Ma  200 个现生植物谱系").italic = True
    doc.add_paragraph()
    add_text(doc, "本报告记录 Case05 v7 的正式全量运行结果。它采用真实 PALEOMAP 稳定 carrier track 的板块搬运、全局球面稀疏扩散转移和现代物种记录的终点条件，输出谱系位置场及扩散条件的谱系支持多样性。报告覆盖 200 个植物谱系、66 个显示时间片和 5 个 HMSC 响应抽样。")
    add_text(doc, "核心结论：扩散过程保留在模型中，但未经历史重复调查、化石过程率或独立人口统计校准的定殖、存续、局地灭绝和最终古占据层被明确关闭。多样性输出因此被命名为扩散条件谱系支持多样性，而不是历史物种丰富度。", "核心结论：")

    add_heading(doc, "一 结果使用说明")
    add_table(doc,
              ["层次", "本次运行实际计算", "不能替代"],
              [
                  ["板块和扩散", "真实 PALEOMAP carrier track 聚合后的板块搬运，加上 Arias 风格球面扩散转移", "自动跨大陆定殖率或真实迁移事件"],
                  ["环境响应", "HMSC posterior Beta 的树上祖先响应重建和古环境线性预测", "最终出现概率或历史占据"],
                  ["扩散多样性", "终点记录格数校准的谱系支持强度、Shannon 混合多样性和有效谱系数", "真实古代植物物种丰富度、PD、FD 或群落 alpha 多样性"],
                  ["未运行过程", "定殖、存续、局地灭绝、最终占据", "零过程或已被数据证明不存在"],
              ], widths=[3.1, 7.5, 5.2])
    add_text(doc, "时间以 Ma 表示，数值越大越古老。325 Ma 早于此批现生谱系树的根年龄之外的时期不作生物学解释。")

    add_heading(doc, "二 方法和结果的连接")
    add_text(doc, "每个谱系的位置场按时间递推：先将上一时间点的空间状态通过板块搬运映射到下一古地理位置，再施加球面扩散转移。现代 0 Ma 的 terminal distribution 被设置为该物种真实记录格的面积加权密度。该设置使深时位置场被现代分布条件化，而不会把现代 HMSC 预测图平均回祖先。")
    add_text(doc, "对位置密度 pi_lc 使用 a_lc = min(1, kappa_l pi_lc) 得到谱系支持场，其中 kappa_l 使 0 Ma 支持格数等于 99.5% 的实际记录格数。随后 I_c = sum_l a_lc，H_c = -sum_l w_lc log(w_lc)，effective_c = exp(H_c)。")
    add_text(doc, "祖先 HMSC 线性响应图与扩散图并列输出，而不相乘。它回答环境响应在空间上的方向和强度；扩散支持图回答在球面扩散和端点条件下哪些地方具有谱系位置支持。")

    add_heading(doc, "三 正式运行质量门")
    add_table(doc, ["检查", "结果", "解释"],
              [[q["check"], q["value"], "通过" if q["pass"].strip().lower() == "true" else "未通过"] for q in quality],
              widths=[6.5, 5.6, 3.7])
    masses = [float(x["terminal_density_mass"]) for x in terminal]
    errors = [abs(float(x["target_terminal_support_cells"]) - float(x["achieved_terminal_support_cells"])) for x in calibration]
    add_text(doc, f"终点条件审计：{len(terminal)} 个现代 tip 的位置密度质量均为 1；终点支持格数校准最大绝对误差为 {max(errors):.3e}。这保证 0 Ma 的支持多样性由实际记录格约束，而不是由短末端枝数值误差决定。")

    add_heading(doc, "四 地图指标索引")
    metric_rows = []
    for metric in METRICS:
        title, definition, reading, limitation = METRIC_INFO[metric]
        ref = index[(metric, 0)]
        metric_rows.append([title, definition, ref["unit"], f"{ref['scale_min']} 至 {ref['scale_max']}"])
    add_table(doc, ["图层", "核心含义", "单位", "全时段固定色标"], metric_rows,
              widths=[4.1, 7.4, 2.6, 2.0])
    add_text(doc, f"全部正式地图的文件级索引为：{MAP_DIR / 'case05_v7_map_index.csv'}。每一行给出指标、时间、统一色标范围、PNG 路径和 GIS 可读 GeoTIFF 路径。")

    add_heading(doc, "五 代表性时间图集")
    add_text(doc, "以下图集按 325、300、250、200、150、100、65、20 和 0 Ma 展示。每一指标在全部 66 个时间片采用固定色标，因而颜色可跨时间直接比较。灰色代表没有对应的有效图层支持或不可用地貌状态，不表示低丰富度。")
    for time_ma in TIME_POINTS:
        doc.add_page_break()
        add_heading(doc, f"时间切片 {time_ma} Ma", level=1)
        if time_ma == 0:
            add_text(doc, "0 Ma 是现代终点：位置与支持图受真实现代记录格条件化。因此它们应与记录分布相一致，但不必与纯环境支持图完全相同。")
        elif time_ma == 325:
            add_text(doc, "325 Ma 接近此抽样树的根时间。此处谱系数较少，低的有效谱系数主要反映活跃祖先 branch 较少，不可与现代 200 个 tip 的数量直接比较。")
        else:
            add_text(doc, "此时间片的空间格局由相邻时间段的板块搬运、球面扩散、地貌权重和树上活跃谱系共同决定。")
        for metric in METRICS:
            row = index.get((metric, time_ma))
            if row is None:
                add_text(doc, f"缺少 {metric} 在 {time_ma} Ma 的正式地图。")
            else:
                add_metric_figure(doc, row, time_ma)

    doc.add_page_break()
    add_heading(doc, "六 如何在论文或报告中表述")
    add_text(doc, "可以写：在一个以现代记录条件化、板块搬运和球面扩散为基础的谱系地理位置模型中，某些古地理区域显示出较高的扩散条件谱系支持强度和支持混合多样性。结果描述的是模型条件下的相对空间支持与谱系混合结构。")
    add_text(doc, "不能写：本模型重建了真实古代植物物种丰富度、真实局地占据率、真实定殖率、真实局地灭绝率或完整历史植物群落。当前现代记录同时服务于 HMSC 与终点条件，因此本运行是模块化 cut model，不是独立现代终点验证。")
    add_text(doc, "下一步若需要将支持多样性升级为古占据或丰富度，需要独立地校准定殖、局地存续和观测过程，例如重复监测、物种级或类群级化石、独立现代范围资料，以及可验证的人口统计或扩散先验。")

    add_heading(doc, "七 输出目录说明")
    add_table(doc, ["目录", "内容"], [
        ["00_contract", "科学合同，逐项说明哪些过程运行、哪些过程未运行。"],
        ["01_input_audit", "物种选择、现代终点记录、终点条件与支持尺度校准审计。"],
        ["02_ancestral_response", "树上重建的祖先 HMSC 环境响应 draw 与均值表。"],
        ["03_global_sphere_transitions", "球面扩散和板块搬运组合转移的归一化质量检查。"],
        ["04_phylogeographic_pruning", "各 HMSC 响应抽样的树条件位置推断对象与对数似然。"],
        ["05_shared_scale_maps", "66 个时间片的 PNG、GeoTIFF 与完整地图索引。"],
        ["06_quality", "正式质量门结果。"],
    ], widths=[4.7, 11.4])

    doc.save(out_docx)
    print(out_docx)


if __name__ == "__main__":
    main()
