#!/usr/bin/env python3
"""Build a Chinese Case05 tutorial from a Case05 output directory.

The report deliberately reads files rather than duplicating model calculations.
It can be rerun after a long Case05 execution to add its map gallery and current
validation status without altering any scientific output.
"""

from __future__ import annotations

import argparse
import csv
import os
import platform
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.style import WD_STYLE_TYPE
from docx.enum.table import WD_ALIGN_VERTICAL
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor
from PIL import Image, ImageDraw, ImageFont


REPORT_NAME = "Case05_Plant200_全过程结果教程_zh.docx"


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def as_number(value: str | None) -> float | None:
    try:
        return float(value) if value not in (None, "", "NA") else None
    except ValueError:
        return None


def set_cell_shading(cell, fill: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shade = tc_pr.find(qn("w:shd"))
    if shade is None:
        shade = OxmlElement("w:shd")
        tc_pr.append(shade)
    shade.set(qn("w:fill"), fill)


def set_cell_border(cell, color: str = "D9D9D9") -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    borders = tc_pr.first_child_found_in("w:tcBorders")
    if borders is None:
        borders = OxmlElement("w:tcBorders")
        tc_pr.append(borders)
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        tag = qn(f"w:{edge}")
        element = borders.find(tag)
        if element is None:
            element = OxmlElement(f"w:{edge}")
            borders.append(element)
        element.set(qn("w:val"), "single")
        element.set(qn("w:sz"), "4")
        element.set(qn("w:space"), "0")
        element.set(qn("w:color"), color)


def set_repeat_table_header(row) -> None:
    tr_pr = row._tr.get_or_add_trPr()
    element = OxmlElement("w:tblHeader")
    element.set(qn("w:val"), "true")
    tr_pr.append(element)


def set_cell_margins(cell, top=90, start=100, bottom=90, end=100) -> None:
    tc = cell._tc
    tc_pr = tc.get_or_add_tcPr()
    margins = tc_pr.first_child_found_in("w:tcMar")
    if margins is None:
        margins = OxmlElement("w:tcMar")
        tc_pr.append(margins)
    for side, value in (("top", top), ("start", start), ("bottom", bottom), ("end", end)):
        node = margins.find(qn(f"w:{side}"))
        if node is None:
            node = OxmlElement(f"w:{side}")
            margins.append(node)
        node.set(qn("w:w"), str(value))
        node.set(qn("w:type"), "dxa")


def add_field(paragraph, field: str) -> None:
    run = paragraph.add_run()
    begin = OxmlElement("w:fldChar")
    begin.set(qn("w:fldCharType"), "begin")
    instr = OxmlElement("w:instrText")
    instr.set(qn("xml:space"), "preserve")
    instr.text = field
    separate = OxmlElement("w:fldChar")
    separate.set(qn("w:fldCharType"), "separate")
    text = OxmlElement("w:t")
    text.text = "1"
    separate.append(text)
    end = OxmlElement("w:fldChar")
    end.set(qn("w:fldCharType"), "end")
    run._r.extend([begin, instr, separate, end])


def add_caption(document: Document, text: str) -> None:
    paragraph = document.add_paragraph(style="Caption")
    paragraph.alignment = WD_ALIGN_PARAGRAPH.LEFT
    run = paragraph.add_run(text)
    run.font.color.rgb = RGBColor(70, 70, 70)


def set_document_baseline(document: Document) -> None:
    section = document.sections[0]
    section.top_margin = Inches(0.72)
    section.bottom_margin = Inches(0.7)
    section.left_margin = Inches(0.78)
    section.right_margin = Inches(0.78)
    styles = document.styles
    normal = styles["Normal"]
    normal.font.name = "Microsoft YaHei"
    normal._element.rPr.rFonts.set(qn("w:eastAsia"), "Microsoft YaHei")
    normal.font.size = Pt(10.5)
    normal.paragraph_format.space_after = Pt(6)
    normal.paragraph_format.line_spacing = 1.22
    for name, size in (("Title", 22), ("Heading 1", 16), ("Heading 2", 13), ("Heading 3", 11)):
        style = styles[name]
        style.font.name = "Microsoft YaHei"
        style._element.rPr.rFonts.set(qn("w:eastAsia"), "Microsoft YaHei")
        style.font.size = Pt(size)
        style.font.bold = True
        style.font.color.rgb = RGBColor(0, 0, 0)
        style.paragraph_format.space_before = Pt(13 if name != "Title" else 0)
        style.paragraph_format.space_after = Pt(6)
    # LibreOffice renders Word's stock Title style with a blue bottom rule.
    # This manual is intentionally utilitarian, so remove that inherited border.
    title_pr = styles["Title"]._element.get_or_add_pPr()
    title_border = title_pr.find(qn("w:pBdr"))
    if title_border is not None:
        title_pr.remove(title_border)
    caption = styles["Caption"]
    caption.font.name = "Microsoft YaHei"
    caption._element.rPr.rFonts.set(qn("w:eastAsia"), "Microsoft YaHei")
    caption.font.size = Pt(8.5)
    caption.font.italic = True
    if "Table Text" not in styles:
        style = styles.add_style("Table Text", WD_STYLE_TYPE.PARAGRAPH)
        style.font.name = "Microsoft YaHei"
        style._element.rPr.rFonts.set(qn("w:eastAsia"), "Microsoft YaHei")
        style.font.size = Pt(8.5)
        style.paragraph_format.space_after = Pt(0)
        style.paragraph_format.line_spacing = 1.05
    footer = section.footer.paragraphs[0]
    footer.alignment = WD_ALIGN_PARAGRAPH.CENTER
    footer.add_run("Case05 Plant200 全过程结果教程  |  页 ")
    add_field(footer, "PAGE")


def add_table(document: Document, headers: list[str], rows: list[list[str]], widths=None) -> None:
    table = document.add_table(rows=1, cols=len(headers))
    table.style = "Table Grid"
    table.autofit = False
    head = table.rows[0]
    set_repeat_table_header(head)
    for index, text in enumerate(headers):
        cell = head.cells[index]
        cell.text = str(text)
        cell.vertical_alignment = WD_ALIGN_VERTICAL.CENTER
        set_cell_shading(cell, "1F4E78")
        set_cell_border(cell)
        for run in cell.paragraphs[0].runs:
            run.font.color.rgb = RGBColor(255, 255, 255)
            run.font.bold = True
            run.font.size = Pt(8.5)
        if widths:
            cell.width = Inches(widths[index])
        set_cell_margins(cell)
    for row_index, values in enumerate(rows):
        cells = table.add_row().cells
        for index, text in enumerate(values):
            cell = cells[index]
            cell.text = str(text)
            cell.vertical_alignment = WD_ALIGN_VERTICAL.CENTER
            if row_index % 2:
                set_cell_shading(cell, "F3F7FA")
            set_cell_border(cell)
            for paragraph in cell.paragraphs:
                paragraph.style = "Table Text"
            if widths:
                cell.width = Inches(widths[index])
            set_cell_margins(cell)
    document.add_paragraph()


def add_bullet(document: Document, text: str, level=0) -> None:
    style = "List Bullet" if level == 0 else "List Bullet 2"
    document.add_paragraph(text, style=style)


def add_step(document: Document, number: int, title: str, question: str,
             input_text: str, action: str, output_text: str,
             interpretation: str, boundary: str) -> None:
    document.add_heading(f"步骤 {number}  {title}", level=2)
    document.add_paragraph(question)
    add_table(document,
              ["输入", "实际计算", "主要输出", "如何读", "不能怎么解释"],
              [[input_text, action, output_text, interpretation, boundary]],
              [1.15, 1.9, 1.4, 1.9, 1.55])


def metric_explanations() -> dict[str, tuple[str, str, str, str]]:
    return {
        "expected_lineage_richness": (
            "期望采样存活谱系丰富度", "活跃谱系最终占据概率之和", "谱系越多、概率越高则值越大", "不是历史全球植物总物种丰富度"),
        "relative_environmental_support_sum": (
            "相对环境支持总量", "各活跃谱系无截距祖先环境响应的支持度之和", "只用于比较 P1 环境过滤在空间和时间上的相对格局", "不是出现概率、不是现代 HMSC nowcast，也不能叫物种丰富度"),
        "binary_lineage_richness": (
            "二值化谱系丰富度", "以验证阈值将占据概率转换为 0 或 1 后计数", "用于稳健性对照，不替代连续概率主图", "阈值化会丢失不确定性"),
        "mean_occupancy_probability": (
            "平均最终占据概率", "每格对当时活跃谱系的边际占据概率平均值", "综合 P1 P2 P4 P5 与局地状态更新", "不是观测到该物种的频率"),
        "mean_environmental_support": (
            "平均祖先环境支持度", "祖先响应乘以该时段古环境后的 probit 或 logit 预测", "这里环境对祖先谱系有多合适", "不是最终古分布概率"),
        "mean_arrival_hazard_per_myr": (
            "平均到达风险率 每百万年", "P2 移动核和来源占据状态形成的到达 hazard", "比较来源限制强弱时优先读此图", "不是 0 到 1 的出现概率"),
        "mean_arrival_probability_reporting_interval": (
            "报告窗口内平均到达概率", "由到达 hazard 按报告窗口转换", "与其他概率图同尺度比较", "窗口长度改变会改变数值"),
        "mean_colonisation_hazard_per_myr": (
            "平均定殖风险率 每百万年", "到达 hazard 乘以条件建立概率", "区分能到达但难建立的地点", "不是实测深时定殖率"),
        "mean_colonisation_probability_reporting_interval": (
            "报告窗口内平均定殖概率", "定殖 hazard 在固定窗口内的转换", "必须和 arrival 图一起看", "低值不能单独断言扩散受限"),
        "mean_persistence_probability": (
            "平均局地存续概率", "已占据 cell 在报告窗口后仍占据的 CTMC 概率", "高值表示模型情景下局地状态较稳定", "不是完整谱系未灭绝概率"),
        "mean_local_extinction_probability": (
            "平均局地消失概率", "局地存续概率的补集", "用于定位局地状态流失热点", "不能叫作全球物种灭绝率"),
        "biotic_competition_pressure": (
            "生物竞争压力", "实证核心模式固定为零", "该层证明 P3 未被伪造", "不能解释为 HMSC 残差竞争"),
        "biotic_facilitation_pressure": (
            "生物促进压力", "实证核心模式固定为零", "该层证明 P3 未被伪造", "不能解释为 HMSC 残差互利"),
        "speciation_inheritance_footprint": (
            "物种形成继承足迹", "定年树节点处 copy then diverge 的空间继承记录", "识别节点导致的谱系身份转换", "不是模型自动推断的真实成种地点"),
        "scenario_lineage_extinction_pressure": (
            "谱系灭绝情景压力", "实证核心不估计全局谱系灭绝", "应为关闭或边界说明层", "不可称作历史灭绝率"),
        "weighted_endemism": (
            "加权特有性", "按当时预测范围面积倒数加权的谱系贡献", "高值指向模型中范围较窄谱系聚集格网", "不等于实际保护优先级"),
        "shannon_diversity": (
            "占据加权 Shannon 熵", "在总占据质量大于容差时，将谱系占据概率归一化后计算 Shannon 熵", "同时反映动态占据的谱系数和均匀度；空 cell 固定为 0", "不是观察丰度 Shannon 指数，也不是完整历史植物群落 Shannon 指数"),
        "simpson_diversity": (
            "占据加权 Gini-Simpson", "在总占据质量大于容差时计算 1-sum(w^2)", "值越大表示动态占据权重更均匀；空 cell 固定为 0", "不是观察丰度 Simpson 指数"),
        "occupancy_weighted_effective_lineages_shannon": (
            "Shannon 有效谱系数", "exp(占据加权 Shannon 熵)", "以等效谱系数呈现权重均匀度", "仅是预测占据权重的描述性摘要"),
        "occupancy_weighted_effective_lineages_simpson": (
            "Simpson 有效谱系数", "1/sum(w^2)", "较强调占据权重占优的谱系结构", "仅在总占据质量大于容差时定义"),
        "occupancy_weighted_diversity_defined": (
            "多样性是否定义", "该 cell 的总动态占据质量是否超过容差", "0 表示没有可归一化的预测谱系 assemblage", "不是生境质量或采样缺失图层"),
        "hmsc_fixed_effect_nowcast_expected_richness": (
            "0 Ma HMSC nowcast 期望丰富度", "使用现代 tip 的原始 HMSC 截距和锁定环境配方，将物种出现概率相加", "应与现代记录和动态终点在大尺度格局上相符", "无 site random effect；必须结合 endpoint score 判断"),
        "dynamic_endpoint_expected_lineage_richness": (
            "0 Ma 动态终点期望谱系丰富度", "合并后验 draw 的 tip-level 动态占据概率之和", "与 0 Ma HMSC nowcast 并列，检查历史动态能否回到合理现代状态", "不是独立验证；当前训练资料仅作 diagnostic"),
        "niche_response_dispersion": (
            "环境响应离散度", "同一格活跃谱系祖先环境响应向量之间的离散程度", "比较生态策略差异与潜在功能分化", "不等于功能多样性或实际竞争强度"),
    }


def safe_file_stem(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in value)


def make_time_atlas(rows: list[dict[str, str]], panel_path: Path) -> Path | None:
    """Create a compact panel that retains every final PNG for one output type."""
    usable = [(row, Path(row.get("png", ""))) for row in rows]
    usable = [(row, png) for row, png in usable if png.exists()]
    if not usable:
        return None
    cols, tile_w, tile_h, label_h, margin = 4, 430, 250, 42, 18
    rows_n = (len(usable) + cols - 1) // cols
    canvas = Image.new("RGB", (cols * tile_w + (cols + 1) * margin,
                                rows_n * (tile_h + label_h) + (rows_n + 1) * margin), "white")
    draw = ImageDraw.Draw(canvas)
    try:
        font = ImageFont.truetype("C:/Windows/Fonts/arial.ttf", 22)
    except OSError:
        font = ImageFont.load_default()
    for index, (row, png) in enumerate(usable):
        image = Image.open(png).convert("RGB")
        image.thumbnail((tile_w, tile_h), Image.Resampling.LANCZOS)
        col, row_n = index % cols, index // cols
        x = margin + col * tile_w + col * margin
        y = margin + row_n * (tile_h + label_h) + row_n * margin
        canvas.paste(image, (x + (tile_w - image.width) // 2, y + (tile_h - image.height) // 2))
        draw.rectangle((x, y, x + tile_w, y + tile_h), outline="#B8B8B8", width=1)
        draw.text((x + 5, y + tile_h + 7), f"{row.get('time_ma', '')} Ma", fill="black", font=font)
    panel_path.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(panel_path, optimize=True)
    return panel_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, help="Case05 output directory")
    parser.add_argument("--docx", default=None, help="Optional output document path")
    args = parser.parse_args()
    output = Path(args.output).resolve()
    report_dir = output / "06_report"
    report_dir.mkdir(parents=True, exist_ok=True)
    docx_path = Path(args.docx).resolve() if args.docx else report_dir / REPORT_NAME
    config = {row.get("parameter", ""): row.get("value", "") for row in read_csv(output / "00_config" / "case05_parallel_posterior_execution.csv")}
    process_rows = read_csv(output / "00_config" / "case05_eight_process_contract.csv")
    validation_rows = read_csv(output / "05_summaries" / "case05_final_validation.csv")
    map_rows = read_csv(output / "05_summaries" / "case05_map_index.csv")
    folder_rows = read_csv(output / "00_config" / "case05_map_folder_codebook.csv")
    difference_rows = read_csv(output / "00_config" / "case05_difference_folder_codebook.csv")
    time_summary = read_csv(output / "05_summaries" / "case05_global_time_scheme_summary.csv")
    jobs = read_csv(output / "01_run_logs" / "case05_posterior_draw_jobs_final.csv")
    current_running = not (output / "05_summaries" / "case05_final_validation.csv").exists()
    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    document = Document()
    set_document_baseline(document)
    document.add_heading("Case05 Plant200 全过程结果教程", level=0)
    subtitle = document.add_paragraph()
    subtitle.alignment = WD_ALIGN_PARAGRAPH.LEFT
    run = subtitle.add_run("全球一度古陆地格网的祖先环境响应和动态占据分析")
    run.font.name = "Microsoft YaHei"
    run._element.rPr.rFonts.set(qn("w:eastAsia"), "Microsoft YaHei")
    run.font.size = Pt(13)
    document.add_paragraph(f"生成时间：{now}。输出目录：{output.as_posix()}")
    status = "正式计算仍在运行，地图和数值表将在任务完成后由本脚本刷新。" if current_running else "正式计算已完成，本教程读取当前输出清单、验证表和地图索引。"
    document.add_paragraph(status)
    document.add_paragraph(
        "本教程面向首次阅读 Case05 的生态学、宏演化、历史生物地理和古环境研究者。它解释每一步做什么、每类输出在哪里、地图颜色和零值意味着什么，以及哪些结果可作为过程约束情景而不能直接写成真实历史事实。"
    )

    document.add_heading("先读这三句话", level=1)
    add_bullet(document, "Case05 是以现代 HMSC 后验环境响应、定年系统树和动态古地球为条件的全球格网正向过程情景；它不是唯一真实历史的后验重建。")
    add_bullet(document, "三套分析只改变 P2 扩散核。P1 环境过滤、P4 演化、P5 定年树节点、根先验、定殖与存续参数、古环境与板块携带均保持一致。")
    add_bullet(document, "BioGeoBEARS、离散区域门控和 A_hist 不在 Case05 主模型中；所有扩散发生在完整古陆地 cell 图上。")
    add_bullet(document, "根年龄约 325 Ma；325 Ma 以前必须被视为 sampled tree 之外的时间域，而不是植物丰富度为零。")

    document.add_heading("一 运行配置和状态", level=1)
    runtime_rows = [
        ["现代数据", "Plant200：200 个现生植物物种；现代群落、性状、定年树和 HMSC posterior。"],
        ["时间域", "325–0 Ma；每个实际古环境层参与积分，代表性 20 个时间片输出地图。"],
        ["空间域", "每个时间片的完整有效一度古陆地格网；海洋和不可用生境均为 NA。"],
        ["后验抽样", f"{config.get('n_real_hmsc_posterior_draws', '25')} 个真实 HMSC response draw，随后等权汇总。"],
        ["并行和存储", f"{config.get('workers', '未记录')} 个工作进程；按扩散方案串行，按 posterior draw 并行，分片合并后清理。"],
        ["板块传输", "PALEOMAP H3 target centred nearest track carriage；它表示地理载运，不是主动扩散。"],
        ["生物相互作用", "实证核心模式关闭，因为没有独立历史相互作用数据。"],
        ["完整谱系灭绝", "现生树无法识别无现生成员的灭绝支系，因此不估计。"],
    ]
    add_table(document, ["主题", "本案例的实际设定"], runtime_rows, [1.35, 5.65])
    if jobs:
        statuses = defaultdict(int)
        for row in jobs:
            statuses[row.get("status", "unknown")] += 1
        document.add_paragraph("当前 posterior-draw 任务状态：" + "；".join(f"{key}={value}" for key, value in sorted(statuses.items())) + "。")

    document.add_heading("二 八过程的唯一职责", level=1)
    process_table = []
    for row in process_rows:
        process_table.append([
            row.get("process", ""), row.get("implementation", ""), row.get("empirical_status", ""),
            "该过程在本案例中只在此处计算；不能再作为额外乘数重复加入最终概率。"
        ])
    if not process_table:
        process_table = [["Environmental filtering", "HMSC posterior beta 与古环境", "data informed", "环境支持不是最终占据。"]]
    add_table(document, ["过程", "本案例的计算位置", "证据状态", "阅读边界"], process_table, [1.25, 2.7, 1.25, 1.8])

    document.add_heading("三 从输入到地图的总流程", level=1)
    document.add_paragraph("计算顺序不可颠倒。现代预测地图不会被向过去平均；被重建的是现代 HMSC 的环境响应参数 beta，再将祖先 response 放到对应年代的古环境上。")
    steps = [
        (1, "读取现代输入与名称审计", "哪些现生谱系和哪些现代环境响应进入案例？", "comm sites traits dated tree HMSC beta posterior", "核对 200 tips、现代训练格网、性状缺失标志与后验 draw 标识。", "01_input_audit 和 02_hmsc_posterior", "物种名称、树 tip 和 posterior draw 必须一一对应。", "GBIF 型现代记录的零不自动等于真实缺失。"),
        (2, "建立活跃谱系时间表", "在每个历史时间点，实际存在的是哪条祖先或后代 branch？", "定年系统树和古环境时间轴", "将树节点插入总时间轴，并仅让 active lineage 进入动态更新。", "03_lineage_time_tree", "分化前不应出现未来 daughter。", "选中的 200 tips 不代表所有历史植物谱系。"),
        (3, "重建祖先环境响应", "祖先适合什么环境？", "每个 HMSC posterior beta draw 与 dated tree", "按 draw 条件重建 branch 内 beta trajectory，再保留 posterior 不确定性。", "04_evolution 中 ancestral beta tables", "这是祖先环境 response surface 的来源。", "不能将现代预测地图平均为祖先地图，也不能把 HMSC 截距等同于祖先生态位。"),
        (4, "读取完整动态古地球", "当时哪些 cell 存在，环境和地形是什么？", "历史环境和一度古陆地切片", "逐时间层读取完整有效格网、环境变量、land habitat 状态和地形阻力。", "05_palaeo_environment 与 06_dynamic_earth_grid", "每个时间片的格网都可不同。", "无陆地 cell 是硬约束，不是低适合度。"),
        (5, "执行地理载运", "大陆移动后上一时间状态如何进入下一古地理位置？", "外部 PALEOMAP H3 transport index", "先做 target centred plate carriage，再构建主动扩散移动图。", "07_plate_transport", "这一步防止把大陆漂移误当主动扩散。", "nearest track 是近似，不是完整板块多边形 overlap。"),
        (6, "环境过滤", "古环境是否支持该祖先谱系？", "祖先 beta 与现代锁定 recipe 变换后的古环境", "计算 eta_env 与 s_env，并将 no analog 作为不确定性语境单列。", "09_environmental_filtering", "s_env 高表示环境支持高。", "s_env 不是最终历史占据概率。"),
        (7, "扩散与到达", "已有来源能否把繁殖体输送到目标 cell？", "上期 occupancy、完整陆地 graph、扩散性状代理和地形成本", "三套 P2 核分别计算 local movement、arrival hazard 及报告窗口概率。", "11_dispersal", "先读 arrival hazard，再读其窗口概率。", "零到达是模型情景结果，不是化石缺失证据。"),
        (8, "定殖", "到达以后是否能建立局地种群？", "arrival 与条件建立模型", "colonisation hazard 等于 arrival hazard 与 conditional establishment 的因果组合。", "12_colonisation", "arrival 低和 establishment 低是两种机制。", "不可称为由现代 HMSC 自动估计的真实深时定殖率。"),
        (9, "存续与局地消失", "已占据 cell 能否留存至下一子步？", "环境支持、救援输入与情景参数", "CTMC 计算 persistence 和 local loss；陆地消失强制清零。", "13_persistence 与 15_extinction", "高 local loss 指出过程情景下的局地不稳定。", "不能将 local loss 写成全球谱系灭绝。"),
        (10, "树条件物种形成", "节点发生后 parent 和 daughters 怎样交接？", "dated tree nodes", "主模式 copy then diverge；daughter 后续各自演化和扩散。", "14_speciation", "节点足迹解释谱系身份变更。", "不代表模型推断了真实成种地理机制。"),
        (11, "动态占据更新", "老种群留下来与新定殖如何合成为下一时刻状态？", "colonisation hazard persistence hazard land state active lineage", "在内部子步执行 CTMC 更新，随后应用 land habitat 与 lineage hard constraints。", "16_state_space_inference", "最终 occupancy 汇集多个过程。", "没有 endpoint fossil likelihood 时，它仍是 forward scenario。"),
        (12, "多样性与过程归因", "各过程累积后留下什么空间格局？", "lineage by cell by time occupancy 与 trait response summaries", "在线计算谱系丰富度、Shannon、Simpson、加权特有性和 response dispersion。", "18_diversity_maps 和 05_summaries", "必须明确为 sampled surviving lineage 指标。", "不代表全部史前植物群落。"),
    ]
    for step in steps:
        add_step(document, *step)

    document.add_heading("四 三种扩散方案如何比较", level=1)
    scheme_rows = [
        ["forward geographic", "陆地八邻域与大圆距离", "地理距离基线；无额外地形阻力", "不是无障碍全球扩散，仍受海陆与距离限制"],
        ["topographic limited", "同一局地图加连续 relief elevation step 阻力", "检验地形增加有效移动成本后的差异", "地形阻力不是观测到的历史迁移路线"],
        ["particle topographic", "同一地形核加有限繁殖体 Monte Carlo 路由", "检验有限繁殖体传播近似的影响", "不是 posterior particle filter，也不是基因流估计"],
    ]
    add_table(document, ["方案", "P2 内核", "用于回答", "不可写成"], scheme_rows, [1.35, 2.05, 1.8, 1.8])
    document.add_paragraph("严格控制：三套方案共享同一输入、同一 HMSC draw、同一祖先 response、同一树、同一根斑块、同一板块传输、同一定殖和存续参数、同一随机种子。只有 P2 movement kernel 改变。")

    document.add_heading("五 常规地图的颜色和读法", level=1)
    document.add_paragraph("每个指标在三个方案与全部实际时间层使用一个固定色标范围；不能把单张图的颜色深浅跨指标直接比较。所有常规 PNG 均有经纬度、Ma 时间、图例和输出属性；同名 GeoTIFF 是 GIS 可读栅格。")
    add_table(document, ["地图颜色", "含义", "应如何处理"], [
        ["深灰", "海洋、不存在的陆地或不可用生境；栅格值为 NA", "不要当作零或低适合度。"],
        ["浅黄", "有效古陆地格网的数值零", "这是模型结果，仍属于完整地图。"],
        ["蓝色梯度", "指标从低到高的共同尺度", "先核对该指标的单位、上限和时间窗口。"],
        ["图例上限", "该指标跨所有方案和时间的共享最大范围，或理论上限", "同一指标可横向比较；不同指标不应混比。"],
    ], [1.2, 3.2, 2.6])

    document.add_heading("六 每一个常规结果指标", level=1)
    explanations = metric_explanations()
    metric_rows = []
    for folder in folder_rows:
        metric = folder.get("metric", "")
        name, calculation, reading, boundary = explanations.get(metric, (metric, "见输出表", "见图例", "需阅读元数据"))
        metric_rows.append([metric, name, calculation, reading, boundary, folder.get("map_folder", "")])
    if not metric_rows:
        for metric, values in explanations.items():
            metric_rows.append([metric, *values, "地图生成后见 map index"])
    add_table(document,
              ["输出列", "中文名称", "它是什么", "图怎么读", "限制", "短目录"],
              metric_rows, [1.15, 1.15, 1.5, 1.35, 1.35, 0.75])
    document.add_paragraph("短目录只为避免 Windows 260 字符路径截断。完整英文指标名、标签、色标范围、PNG 与 GeoTIFF 完整路径均在 05_summaries/case05_map_index.csv 中。")

    document.add_heading("七 差值地图如何读", level=1)
    document.add_paragraph("差值固定为 numerator minus denominator。上位目录的比较方向由 case05_difference_folder_codebook.csv 定义，不能只凭颜色猜测。每张差值图配有 comparison state GeoTIFF，保证零值、相等值和缺失值可机器读取。")
    add_table(document, ["颜色或状态", "含义", "对 colonisation 的正确解释"], [
        ["红", "numerator 高于 denominator", "说明指定窗口内定殖概率在 numerator 情景更高；仍需核对 arrival 与 establishment。"],
        ["蓝", "numerator 低于 denominator", "说明指定窗口内定殖概率在 numerator 情景更低；不能直接等同为灭绝。"],
        ["白", "两方案相等且非零", "两种 P2 方案在该指标无数值差异。"],
        ["中灰", "两方案都为有效零", "是有效古陆地上的共同零，不是被切掉的格网。"],
        ["深灰", "海洋或不可用生境", "不可比较；没有过程数值。"],
    ], [1.2, 2.0, 3.8])
    if difference_rows:
        add_caption(document, "表  差值目录代码表")
        add_table(document, ["比较方向", "短目录", "文件代码"], [[r.get("comparison_id", ""), r.get("map_folder", ""), r.get("file_code", "")] for r in difference_rows], [2.5, 2.0, 2.5])

    document.add_heading("八 关键审核表和输出文件", level=1)
    output_rows = [
        ["00_config/case05_parallel_posterior_execution.csv", "运行规格", "确认 200 物种、25 draws、时间选择、内部步长、板块与根先验。"],
        ["00_config/case05_eight_process_contract.csv", "过程合同", "确认八过程的主计算位置及 data informed、scenario calibrated 或 off 状态。"],
        ["00_config/case05_shared_metric_limits.csv", "色标合同", "验证每个常规指标在全部方案和时间共享范围。"],
        ["02_scenario_runs/<scheme>/16_state_space_inference/cell_metrics_<time>Ma.csv", "每时间片完整格网表", "每个有效古陆地 cell 的过程与多样性数值；正式汇总结果。"],
        ["03_shared_scale_maps", "常规地图", "按扩散方案和指标保存 GeoTIFF 与 PNG。"],
        ["04_scheme_difference_maps", "差值地图", "用于隔离 P2 方案影响；含 comparison state GeoTIFF。"],
        ["05_summaries/case05_colonisation_zero_arrival_audit.csv", "到达定殖因果审核", "必须满足 zero arrival implies zero colonisation。"],
        ["05_summaries/case05_difference_state_audit.csv", "差值状态审核", "逐图统计 both zero、equal nonzero、different cell 数。"],
        ["05_summaries/case05_final_validation.csv", "最终质量门", "应优先查看是否全为 PASS；失败项先修复，不写结论。"],
        ["05_summaries/case05_global_time_scheme_summary.csv", "全球时间轨迹", "比较方案间总丰富度、平均过程量和时间趋势。"],
        ["05_summaries/case05_map_index.csv", "地图索引", "每一张图的完整路径、指标、时间、色标范围、CRS 和来源。"],
    ]
    add_table(document, ["位置", "内容", "使用方式"], output_rows, [3.15, 1.35, 2.5])

    document.add_heading("九 当前质量控制", level=1)
    if validation_rows:
        add_table(document, ["检查", "状态", "解释"], [[r.get("check", ""), r.get("status", ""), r.get("detail", "")] for r in validation_rows], [2.1, 0.75, 4.15])
    else:
        document.add_paragraph("正式任务仍在进行，尚未生成最终质量门。完成后重新运行本脚本，教程会自动插入 case05_final_validation.csv。")
    document.add_paragraph("最低验收条件：每方案时间轴一致；每个有效古陆地 cell 都进入输出；概率在 [0, 1]；新陆地不自动继承；未来 daughter 不在分化前出现；零 arrival 必导致零 colonisation；所有被索引地图真实存在；每个常规指标色标一致。")

    document.add_heading("十 论文级解释边界", level=1)
    boundaries = [
        ["可以写", "在给定古环境、根先验、板块携带、树条件谱系身份和过程参数的情景下，某方案产生了较高或较低的 expected sampled surviving lineage richness。"],
        ["可以写", "祖先环境 response posterior 与古环境结合后，在某些 cell 显示较高的潜在环境支持；该图不等同于实际古分布。"],
        ["可以写", "在相同输入下，加入连续地形阻力改变了到达、定殖、局地存续及下游多样性预测。"],
        ["不能写", "本案例恢复了全球全部植物的真实 325 Myr 多样性历史。"],
        ["不能写", "HMSC residual association 证明了古代竞争或促进；实证核心已关闭 P3。"],
        ["不能写", "局地消失概率就是完整谱系灭绝率；现生树不能恢复未留下后代的灭绝谱系。"],
        ["不能写", "PALEOMAP target centred nearest track 是精确板块多边形传输或真实迁徙轨迹。"],
    ]
    add_table(document, ["措辞", "推荐或禁止的表述"], boundaries, [1.0, 6.0])

    document.add_heading("十一 如何从图到结论", level=1)
    document.add_paragraph("建议按固定顺序读，而不是先看 richness。首先看 land habitat 图确认该时段舞台；其次看 mean environmental support 区分环境支持；再次看 arrival hazard 与 arrival probability 判断来源限制；随后看 colonisation 和 persistence 分开建立与局地留存；最后读取 occupancy、expected lineage richness、Shannon、Simpson、weighted endemism 和 response dispersion。比较方案时先看相同时间、相同指标的差值图，再回查输入过程图。")
    add_bullet(document, "某 cell 的环境支持高而到达低：称为潜在适生但来源受限，不称为历史缺失。")
    add_bullet(document, "到达非零但定殖接近零：称为在当前建立参数下可能 establishment limited，不称为绝对无法定殖。")
    add_bullet(document, "occupancy 低且 local loss 高：应检查环境、救援、地理硬约束和时间步敏感性。")
    add_bullet(document, "richness 的显著差异：必须回溯是 lineage identity、环境支持、arrival、colonisation、persistence 还是强制地理丢失导致。")

    document.add_heading("十二 地图图集和结果索引", level=1)
    if not map_rows:
        document.add_paragraph("正式地图尚未生成。正式计算完成后再次运行本脚本；它将从 case05_map_index.csv 自动插入每个指标的代表时间地图与对应解释。")
    else:
        panels = []
        atlas_dir = report_dir / "_tutorial_time_atlases"
        by_key = defaultdict(list)
        for row in map_rows:
            by_key[(row.get("scenario_id", ""), row.get("metric", ""))].append(row)
        for (scenario, metric), rows in sorted(by_key.items()):
            rows.sort(key=lambda r: as_number(r.get("time_ma")) or -999, reverse=True)
            panel = make_time_atlas(rows, atlas_dir / f"{safe_file_stem(scenario)}__{safe_file_stem(metric)}.png")
            if panel:
                panels.append((scenario, metric, rows, panel))
        if not panels:
            document.add_paragraph("地图索引已存在，但尚没有可嵌入的 PNG；请检查 map index 与路径。")
        else:
            document.add_paragraph("每幅面板按时间排列并包含该方案与指标的全部导出 PNG。原始单图仍保存在 03_shared_scale_maps 或 04_scheme_difference_maps；map index 提供其完整路径、CRS 和数值范围。")
            for scenario, metric, rows, panel in panels:
                base_metric = metric.removesuffix("_difference")
                title, calculation, reading, boundary = explanations.get(base_metric, (metric, "", "", ""))
                document.add_heading(f"{scenario}  {title}", level=3)
                document.add_picture(str(panel), width=Inches(6.55))
                lower = rows[0].get("lower", "")
                upper = rows[0].get("upper", "")
                add_caption(document, f"时间图集：{metric}；共 {len(rows)} 个时间片；统一色标 [{lower}, {upper}]。")
                document.add_paragraph(f"读法：{reading}。限制：{boundary}。")

    document.add_heading("十三 完成后刷新本教程", level=1)
    document.add_paragraph("正式运行完成后，在包根目录执行下列命令即可用最终 map index、质量门和 PNG 图集覆盖本文件：")
    command = document.add_paragraph()
    command.style = "Normal"
    run = command.add_run(f"python scripts/build_case05_tutorial_docx.py --output \"{output.as_posix()}\"")
    run.font.name = "Consolas"
    run.font.size = Pt(9)
    document.add_paragraph("该刷新操作只读取已有输出，不重新计算 HMSC、祖先 response、扩散、定殖或多样性。")

    document.save(docx_path)
    print(docx_path)


if __name__ == "__main__":
    main()
