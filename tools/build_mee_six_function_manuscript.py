"""Build the six-function MEE manuscript from the audited Plant200 draft.

The source DOCX is a scientific source and supplies four audited figures. The
new document is a substantive rewrite, not a conversion or a figure-only edit.
"""

from __future__ import annotations

import csv
from io import BytesIO
from pathlib import Path

from docx import Document
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(r"C:\Users\Google\Documents\HMSC-HIST")
SOURCE = Path(
    r"C:\Users\Google\Desktop\HMSCEE"
    r"\HmscEE_MEE_complete_manuscript_landmask_audited.docx"
)
OUTPUT = Path(
    r"C:\Users\Google\Desktop\HMSCEE"
    r"\HmscEE_MEE_six_function_manuscript.docx"
)
FIGURE_DIR = ROOT / "HmscEcoEvo" / "outputs" / "mee_six_function_manuscript"
MARINE_CSV = (
    ROOT / "HmscEcoEvo" / "docs" / "manuscript"
    / "MEE_full_manuscript_20260924" / "source_data"
    / "geographic_landmask_audit_all_times.csv"
)


def source_image(doc: Document, index: int) -> BytesIO:
    inline = doc.inline_shapes[index]._inline
    rel = inline.graphic.graphicData.pic.blipFill.blip.embed
    return BytesIO(doc.part.related_parts[rel].blob)


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    name = "arialbd.ttf" if bold else "arial.ttf"
    return ImageFont.truetype(str(Path(r"C:\Windows\Fonts") / name), size)


def centered(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int],
             title: str, detail: str, fill: str, stroke: str) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=19, fill=fill, outline=stroke, width=3)
    title_font = font(28, True)
    detail_font = font(21)
    tw = draw.textlength(title, font=title_font)
    dw = draw.textlength(detail, font=detail_font)
    draw.text(((x0 + x1 - tw) / 2, y0 + 28), title, font=title_font,
              fill="#203449")
    draw.text(((x0 + x1 - dw) / 2, y0 + 69), detail, font=detail_font,
              fill="#435466")


def arrow(draw: ImageDraw.ImageDraw, x: int, y0: int, y1: int) -> None:
    draw.line((x, y0, x, y1 - 12), fill="#648095", width=5)
    draw.polygon([(x - 10, y1 - 17), (x + 10, y1 - 17), (x, y1)],
                 fill="#648095")


def draw_framework(path: Path) -> None:
    im = Image.new("RGB", (1900, 870), "white")
    d = ImageDraw.Draw(im)
    title = font(30, True)
    d.text((62, 22), "Evidence and dynamic Earth inputs", fill="#203449", font=title)
    boxes = [
        ((65, 80, 615, 202), "Modern communities", "HMSC response draws"),
        ((675, 80, 1225, 202), "Dated tree and traits", "Active lineages and ancestry"),
        ((1285, 80, 1835, 202), "Palaeo Earth grid", "Climate, land and plate carriage"),
    ]
    for box, a, b in boxes:
        centered(d, box, a, b, "#F2F7FA", "#7EA4BE")
    arrow(d, 38, 207, 258)
    d.text((62, 273), "Six biological function families", fill="#203449", font=title)
    parts = [
        (65, 325, "Evolution", "Ancestral response"),
        (675, 325, "Speciation", "Tree node events"),
        (1285, 325, "Environmental filtering", "Environmental support"),
        (65, 473, "Dispersal", "Landscape transition"),
        (675, 473, "Biotic filtering", "External interactions only"),
        (1285, 473, "Extinction", "Loss and survival limits"),
    ]
    for x, y, a, b in parts:
        fill = "#F0F7F1" if a in {"Evolution", "Speciation", "Environmental filtering"} else "#F4F3F8"
        stroke = "#83AA91" if fill == "#F0F7F1" else "#A49BBC"
        centered(d, (x, y, x + 550, y + 122), a, b, fill, stroke)
    arrow(d, 38, 603, 655)
    d.text((62, 669), "Distinct outputs and evidence status", fill="#203449", font=title)
    for box, a, b in [
        ((65, 713, 615, 849), "Environmental support", "Calculated from response and climate"),
        ((675, 713, 1225, 849), "Geographic support", "Scenario and endpoint conditioned"),
        ((1285, 713, 1835, 849), "Occupancy and diversity", "Not estimated in this application"),
    ]:
        centered(d, box, a, b, "#FFF8EF", "#D8B47F")
    im.save(path, dpi=(250, 250))


def draw_marine_audit(path: Path) -> None:
    with MARINE_CSV.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 66
    assert sum(float(row["marine_location_fraction"]) > 0.5 for row in rows) == 59
    im = Image.new("RGB", (1800, 710), "white")
    d = ImageDraw.Draw(im)
    x0, x1, y0, y1 = 175, 1680, 90, 550

    def xy(age: float, fraction: float) -> tuple[float, float]:
        return (x0 + (325 - age) / 325 * (x1 - x0),
                y1 - fraction * (y1 - y0))

    for level in (0.0, 0.25, 0.5, 0.75, 1.0):
        y = int(xy(325, level)[1])
        color = "#637487" if level == 0.5 else "#DEE4E9"
        d.line((x0, y, x1, y), fill=color, width=3 if level == 0.5 else 2)
        d.text((65, y - 17), f"{int(level * 100)}%", font=font(27), fill="#33495A")
    for age in (325, 300, 250, 200, 150, 100, 50, 0):
        x = int(xy(age, 0)[0])
        d.line((x, y0, x, y1), fill="#EFF2F5", width=2)
        d.text((x - 19, 568), str(age), font=font(25), fill="#33495A")
    d.line((x0, y0, x0, y1), fill="#263847", width=3)
    d.line((x0, y1, x1, y1), fill="#263847", width=3)
    for field, color in (
        ("marine_location_fraction", "#AE5B31"),
        ("marine_support_fraction", "#2C6489"),
    ):
        points = [xy(float(row["time_ma"]), float(row[field])) for row in rows]
        d.line(points, fill=color, width=6, joint="curve")
        for x, y in points:
            d.ellipse((x - 3, y - 3, x + 3, y + 3), fill=color)
    d.text((x0, 25), "Marine-state mass across all 66 mapped times",
           font=font(36, True), fill="#203449")
    d.text((660, 624), "Time before present Ma", font=font(26), fill="#33495A")
    d.line((240, 671, 290, 671), fill="#AE5B31", width=6)
    d.text((302, 653), "Location mass", font=font(25), fill="#33495A")
    d.line((635, 671, 685, 671), fill="#2C6489", width=6)
    d.text((698, 653), "Calibrated support", font=font(25), fill="#33495A")
    d.line((1150, 671, 1200, 671), fill="#637487", width=3)
    d.text((1213, 653), "50% diagnostic threshold", font=font(25), fill="#33495A")
    im.save(path, dpi=(250, 250))


def set_cell_shading(cell, color: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:fill"), color)
    tc_pr.append(shd)


def set_cell_borders(cell) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    borders = OxmlElement("w:tcBorders")
    for edge in ("top", "left", "bottom", "right"):
        item = OxmlElement("w:" + edge)
        item.set(qn("w:val"), "single")
        item.set(qn("w:sz"), "5")
        item.set(qn("w:color"), "D9D9D9")
        borders.append(item)
    tc_pr.append(borders)


def set_cell_margin(cell, top: int = 80, start: int = 90,
                    bottom: int = 80, end: int = 90) -> None:
    tc = cell._tc
    tc_pr = tc.get_or_add_tcPr()
    mar = tc_pr.first_child_found_in("w:tcMar")
    if mar is None:
        mar = OxmlElement("w:tcMar")
        tc_pr.append(mar)
    for side, value in (("top", top), ("start", start),
                        ("bottom", bottom), ("end", end)):
        el = mar.find(qn("w:" + side))
        if el is None:
            el = OxmlElement("w:" + side)
            mar.append(el)
        el.set(qn("w:w"), str(value))
        el.set(qn("w:type"), "dxa")


def add_table(doc: Document, headers: list[str], rows: list[list[str]],
              widths: list[float]) -> None:
    tab = doc.add_table(rows=1, cols=len(headers))
    tab.autofit = False
    for i, value in enumerate(headers):
        tab.columns[i].width = Inches(widths[i])
        cell = tab.rows[0].cells[i]
        cell.width = Inches(widths[i])
        cell.text = value
        set_cell_shading(cell, "E8EDF1")
        for run in cell.paragraphs[0].runs:
            run.bold = True
            run.font.size = Pt(8.8)
    tr_pr = tab.rows[0]._tr.get_or_add_trPr()
    repeat = OxmlElement("w:tblHeader")
    repeat.set(qn("w:val"), "true")
    tr_pr.append(repeat)
    for row_i, values in enumerate(rows):
        cells = tab.add_row().cells
        for i, value in enumerate(values):
            cells[i].width = Inches(widths[i])
            cells[i].text = value
            if row_i % 2:
                set_cell_shading(cells[i], "F8FAFB")
    for row in tab.rows:
        for cell in row.cells:
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            set_cell_margin(cell)
            set_cell_borders(cell)
            for paragraph in cell.paragraphs:
                paragraph.paragraph_format.space_after = Pt(0)
                paragraph.paragraph_format.line_spacing = 1.08
                for run in paragraph.runs:
                    run.font.name = "Arial"
                    if run.font.size is None:
                        run.font.size = Pt(8.6)
    doc.add_paragraph()


def add_page_number(paragraph) -> None:
    paragraph.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    run = paragraph.add_run("HmscEE methods manuscript  |  ")
    run.font.size = Pt(8)
    field = OxmlElement("w:fldSimple")
    field.set(qn("w:instr"), "PAGE")
    paragraph._p.append(field)


def set_styles(doc: Document) -> None:
    sec = doc.sections[0]
    sec.page_width = Inches(8.5)
    sec.page_height = Inches(11)
    sec.top_margin = Inches(0.8)
    sec.bottom_margin = Inches(0.75)
    sec.left_margin = Inches(0.86)
    sec.right_margin = Inches(0.86)
    sec.header_distance = Inches(0.35)
    sec.footer_distance = Inches(0.4)
    normal = doc.styles["Normal"]
    normal.font.name = "Cambria"
    normal.font.size = Pt(10.5)
    normal.font.color.rgb = RGBColor(0, 0, 0)
    normal.paragraph_format.line_spacing = 1.17
    normal.paragraph_format.space_after = Pt(7)
    for style_name, size, before, after in [
        ("Title", 16.5, 0, 13),
        ("Heading 1", 13.5, 18, 7),
        ("Heading 2", 11.5, 12, 5),
    ]:
        style = doc.styles[style_name]
        style.font.name = "Cambria"
        style.font.size = Pt(size)
        style.font.bold = True
        style.font.color.rgb = RGBColor(0, 0, 0)
        style.paragraph_format.space_before = Pt(before)
        style.paragraph_format.space_after = Pt(after)
        style.paragraph_format.keep_with_next = True
    title_pr = doc.styles["Title"]._element.get_or_add_pPr()
    title_border = title_pr.find(qn("w:pBdr"))
    if title_border is not None:
        title_pr.remove(title_border)
    add_page_number(sec.footer.paragraphs[0])


def prose(doc: Document, text: str) -> None:
    p = doc.add_paragraph()
    labels = (
        "Environmental filtering.", "Dispersal.", "Biotic filtering.",
        "Evolution.", "Speciation.", "Extinction.",
    )
    prefix = next((label for label in labels if text.startswith(label)), None)
    if prefix:
        p.add_run(prefix).bold = True
        p.add_run(text[len(prefix):])
    else:
        p.add_run(text)
    p.paragraph_format.widow_control = True


def heading(doc: Document, text: str, level: int = 1) -> None:
    doc.add_heading(text, level)


def caption(doc: Document, text: str) -> None:
    p = doc.add_paragraph()
    p.style = doc.styles["Normal"]
    p.paragraph_format.space_before = Pt(4)
    p.paragraph_format.space_after = Pt(11)
    p.paragraph_format.keep_together = True
    run = p.add_run(text)
    run.font.size = Pt(9.2)


def figure(doc: Document, image: BytesIO | Path, width: float,
           text: str) -> None:
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before = Pt(6)
    p.paragraph_format.space_after = Pt(2)
    p.paragraph_format.keep_with_next = True
    descriptor = str(image) if isinstance(image, Path) else image
    p.add_run().add_picture(descriptor, width=Inches(width))
    caption(doc, text)


def build() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    framework_path = FIGURE_DIR / "Figure_1_six_function_framework.png"
    marine_path = FIGURE_DIR / "Figure_5_marine_state_audit.png"
    draw_framework(framework_path)
    draw_marine_audit(marine_path)
    source = Document(SOURCE)
    doc = Document()
    set_styles(doc)
    doc.core_properties.title = (
        "HmscEE A modular R framework for linking contemporary communities "
        "to dynamic Earth histories"
    )
    doc.core_properties.subject = "Methods in Ecology and Evolution manuscript draft"

    p = doc.add_paragraph(style="Title")
    p.add_run(
        "HmscEE A modular R framework for linking contemporary communities "
        "to dynamic Earth histories"
    )
    prose(doc, "Methods and Resources manuscript draft for Methods in Ecology and Evolution")

    heading(doc, "Abstract")
    prose(doc,
        "Contemporary community data record the joint consequences of environmental "
        "responses, movement and evolutionary history, but they do not identify those "
        "processes on their own. Existing workflows often combine modern species "
        "distribution models, phylogenies and palaeoenvironmental maps without "
        "distinguishing environmental support from historical location or occupancy. "
        "We introduce HmscEE, implemented in the R package HmscEcoEvo, as a modular "
        "workflow that keeps these inferential targets separate. Its present conceptual "
        "configuration has six biological function families: environmental filtering, "
        "dispersal, biotic filtering, evolution, speciation and extinction. Dynamic "
        "climate and geography provide external inputs rather than extra biological "
        "probability factors. We describe each function's data requirements, outputs "
        "and evidence status, then demonstrate the workflow with 200 extant plants "
        "from 68 families. Modern HMSC responses are projected through a dated tree "
        "onto 66 palaeoenvironmental slices, and a separate geographic scenario is "
        "audited on a four-degree grid. The held-out modern occurrence diagnostic gave "
        "a Brier score of 0.0181, but a terrestrial-state check found that 59 of 66 "
        "historical slices placed more than half of geographic location mass in marine "
        "transit states. These maps therefore demonstrate an auditable workflow, not "
        "a validated terrestrial occupancy history. By reporting fitted estimates, "
        "conditional reconstructions and uncalibrated scenarios separately, HmscEE "
        "supports more reproducible tests of how biodiversity may have assembled "
        "on a changing Earth."
    )
    prose(doc,
        "Keywords  ancestral environmental response; community ecology; dynamic "
        "palaeogeography; HMSC; phylogeny; spatial inference; uncertainty"
    )

    heading(doc, "Introduction")
    prose(doc,
        "Patterns of biodiversity reflect both contemporary conditions and histories "
        "that cannot be read directly from a modern map. Temperature and moisture may "
        "support a lineage in places it never reached; an occupied place may become "
        "unsuitable or disappear; and lineage identity changes at speciation nodes. "
        "Historical analysis therefore needs to distinguish ecological response, "
        "geographic movement and the evidence used to constrain them. Joint species "
        "distribution models estimate species-specific responses to measured "
        "environments, while phylogenies and palaeo-Earth data place those responses "
        "in temporal and spatial context (Ovaskainen et al., 2017; Tikhonov et al., "
        "2020; Guillory & Brown, 2021). None of these sources alone identifies a full "
        "history of local population occupancy."
    )
    prose(doc,
        "gen3sis demonstrates the value of a common engine with explicit landscapes, "
        "configurable biological functions and observable outputs (Hagen et al., "
        "2021). HmscEE adopts that clarity of organisation, but addresses a different "
        "inferential task. gen3sis primarily explores outcomes generated by specified "
        "mechanisms; HmscEE links a fitted modern community model to dated lineages "
        "and palaeoenvironmental inputs, while making clear which historical quantities "
        "are reconstructed, conditioned or only explored as scenarios. Its six "
        "biological function families are not gen3sis's six configuration functions "
        "and should not be treated as interchangeable."
    )
    prose(doc,
        "A particular risk in deep-time workflows is to multiply several measures "
        "that describe the same movement event. Another is to average modern "
        "prediction maps backwards in time, thereby carrying today's climate and "
        "geography into an ancestral response. HmscEE instead reconstructs ancestral "
        "environmental coefficients from HMSC posterior draws and the dated tree, "
        "then evaluates those coefficients against time-matched palaeoenvironmental "
        "grids. Geographic movement is represented separately. Geological carriage "
        "changes the location of a land unit; it is not itself organismal dispersal. "
        "Where a regional history is supplied, it conditions movement at that scale "
        "rather than adding an independent accessibility multiplier. The Plant200 "
        "illustration described here uses the region-free global-grid route and does "
        "not run BioGeoBEARS."
    )
    prose(doc,
        "We present the inputs and six functions as a reproducible configuration, "
        "describe the corresponding output classes, and use a worked plant analysis "
        "to show both what the current workflow produces and where it fails an "
        "essential terrestrial validation test. The aim is not to infer every "
        "historical plant population from one extant tree. It is to give ecological "
        "and historical evidence a defined role, so that more strongly identified "
        "analyses can be built and tested."
    )

    heading(doc, "Framework and implementation")
    heading(doc, "Inputs configuration and state", 2)
    prose(doc,
        "HmscEE accepts a modern community matrix, matched contemporary environmental "
        "covariates and coordinates, a dated phylogeny, optional measured traits, and "
        "a sequence of palaeo-Earth grids. Each historical grid records time, cell "
        "identity, coordinates, environmental covariates and habitat or land status. "
        "Movement analyses additionally need landscape costs or a sparse movement "
        "graph and, when continents move, a plate-carriage correspondence between "
        "successive grids. Fossils and independent modern observations may constrain "
        "or validate outputs at their actual taxonomic and temporal resolution."
    )
    prose(doc,
        "The configuration declares the target of each run before computation. "
        "Environmental support describes the response of an active lineage to a "
        "specified environment. Geographic-location support describes uncertainty "
        "about a lineage's position under a movement model and endpoint condition. "
        "Demographic occupancy describes whether a population occupies a cell. "
        "A source-normalised location density is not an occupancy probability, and "
        "its sum across lineages is not expected richness. This distinction governs "
        "which maps and diversity statistics can be labelled as empirical results."
    )
    prose(doc,
        "The present manuscript uses six top-level biological function families. "
        "Arrival, establishment and local retention may appear inside a future "
        "calibrated occupancy transition, but colonisation and persistence are not "
        "separate headline functions in this revised design. Climate, geology, "
        "geography, habitat, the dated tree and observation data are inputs or "
        "constraints, not additional biological functions (Figure 1; Table 1)."
    )
    figure(doc, framework_path, 6.7,
           "Figure 1. HmscEE evidence flow and six biological function families. "
           "The functions are interfaces with different evidence status, not six "
           "independently fitted probabilities. Plate carriage belongs to the Earth "
           "input layer. The present Plant200 analysis estimates modern responses "
           "and constructs conditional environmental and geographic support; it "
           "does not estimate demographic occupancy.")

    caption(doc,
        "Table 1. Public function examples, responsibilities and Plant200 status. "
        "A function's availability does not imply its parameters were estimated.")
    add_table(doc,
        ["Function family", "Representative public entry point", "Role and status in this analysis"],
        [
            ["Environmental filtering", "hee_environmental_filtering_ancestral_suitability()",
             "Projects reconstructed response coefficients onto historical environments; environmental support calculated."],
            ["Dispersal", "hee_dispersal_spherical_kernel()",
             "Builds a distance and cost dependent movement kernel; geographic movement is a fixed scenario, not a fitted rate."],
            ["Biotic filtering", "hee_biotic_filtering_effect()",
             "Accepts an externally justified interaction effect; set to zero without independent interaction evidence."],
            ["Evolution", "hee_evolution_ancestral_response_direct()",
             "Reconstructs ancestral environmental responses from tip-level HMSC draws and a dated tree; approximate in this run."],
            ["Speciation", "hee_speciation_events()",
             "Reads dated-tree splitting events and changes active lineage identity; no novel speciation rate is fitted."],
            ["Extinction", "hee_extinction_layers()",
             "Provides a loss-summary interface for compatible dynamic states; full lineage extinction not estimated here."],
        ],
        [1.25, 2.45, 2.95])

    heading(doc, "The six biological functions", 2)
    prose(doc,
        "Environmental filtering. The HMSC fit supplies posterior draws of the "
        "modern species' environmental coefficients. "
        "hee_environmental_filtering_ancestral_suitability() applies the modern "
        "predictor recipe to historical environmental cells, obtains the response "
        "of each lineage active at the requested time and returns a link-scale "
        "environmental score and its inverse-link support. The fitted intercept and "
        "sampling-effort term are not silently interpreted as ancestral physiological "
        "tolerance. Support is conditional on the specified response model and "
        "available environment; it is not a historical occurrence probability."
    )
    prose(doc,
        "Dispersal. hee_dispersal_spherical_kernel() represents organismal movement "
        "among cells through a declared distance and landscape cost. Related graph "
        "and transition functions can account for changing connectivity, while "
        "plate-carriage correspondences move geographic reference units between "
        "time slices. These operators must remain distinct: plate movement does "
        "not demonstrate seed movement, and a landscape kernel does not estimate "
        "establishment unless a separate model and data support that interpretation. "
        "The Plant200 calculation uses a fixed four-degree, eight-neighbour "
        "geographic movement scenario and a low marine-transit weight; neither "
        "parameter was calibrated from historical plant movement."
    )
    prose(doc,
        "Biotic filtering. hee_biotic_filtering_effect() accepts an independently "
        "supplied interaction effect when the relevant partners and their effects "
        "are observed or separately modelled. Its default value is zero when such "
        "evidence is absent. HMSC residual association is not automatically a "
        "competition or facilitation coefficient because shared unmeasured "
        "environments and observation processes can also create association "
        "(Blanchet et al., 2020). Plant200 therefore does not infer ancient biotic "
        "interactions."
    )
    prose(doc,
        "Evolution. hee_evolution_ancestral_response_direct() takes selected draws "
        "of modern HMSC coefficients as tip-level responses and reconstructs "
        "time-specific responses along the dated tree. This reconstructs response "
        "functions, not modern prediction maps. A trait-mediated route can relate "
        "ancestral traits to environmental responses, but requires trait measurements, "
        "missing-data uncertainty and a justified trait-response model. In this "
        "plant analysis, coefficients were treated separately and branch values "
        "received marginal Brownian-bridge perturbations around fitted node "
        "estimates. The result is a conditional approximation, not a joint "
        "multivariate evolutionary posterior."
    )
    prose(doc,
        "Speciation. hee_speciation_events() standardises the splitting times "
        "and parent-daughter identities in the supplied dated tree. A parent "
        "lineage stops at its node, and its daughters become active thereafter. "
        "This prevents a modern daughter species being projected into a time "
        "before its origin. The function does not infer a new macroevolutionary "
        "speciation rate or prove a particular geographic mode of speciation."
    )
    prose(doc,
        "Extinction. hee_extinction_layers() can summarise loss states when an "
        "appropriate dynamic state and extinction input are available. Habitat "
        "disappearance is a hard geographic constraint; it is not an independently "
        "estimated extinction event. More importantly, a tree containing only "
        "surviving tips cannot recover all lineages that left no extant descendants. "
        "Plant200 has no fitted global lineage-extinction process. Extinction "
        "maps or rates would require additional fossil, extinct-tip or demographic "
        "evidence, and the supported layer must be stated explicitly."
    )

    heading(doc, "Historical predictors and output recording", 2)
    prose(doc,
        "The retained hmscHist pathway is a supporting data workflow, not a "
        "seventh biological process. hee_history_data(), "
        "hee_history_build_tables() and related conversion functions align "
        "species, sites, traits and imported historical summaries for a modern "
        "HMSC design. Historical variables computed from the same response "
        "community must be labelled as response-derived descriptors; using them "
        "as if they were independent drivers can leak information into a fit or "
        "cross-validation fold. The small three-site, two-species example in the "
        "project documentation is synthetic and tests this interface only."
    )
    prose(doc,
        "An observer-style output layer records the active lineage set, environmental "
        "support, geographic-location support, diagnostic flags and scenario "
        "metadata for each mapped time. Outputs that require calibrated occupancy "
        "states, such as expected richness, time-sliced phylogenetic diversity, "
        "functional diversity and turnover, are not calculated from location "
        "density by relabelling or rescaling it. Nonlinear summaries must be "
        "computed within joint state draws before posterior aggregation."
    )

    heading(doc, "Worked Plant200 application")
    heading(doc, "Data and modern observation model", 2)
    prose(doc,
        "The worked analysis selects 200 extant terrestrial plant species from "
        "68 families. Its prepared occurrence matrix contains 2,168 occupied "
        "one-degree cells. After aggregation and environmental matching, the "
        "four-degree HMSC analysis used 554 cells, with 432 in training and 122 "
        "in a spatially held-out diagnostic. The 200-tip V.PhyloMaker2 S3 "
        "backbone has an original root age of about 325.05 Ma; the geographic "
        "runner made a small numerical rescaling to 325 Ma to match its oldest "
        "carrier. This is not a new fossil recalibration (Jin & Qian, 2022). "
        "Historical environments extend from 540 to 0 Ma, but lineage analyses "
        "before the sampled-tree root are not applicable rather than zero. "
        "Ten core traits have observed coverage from 59.0% to 97.5%, although "
        "traits did not enter the fitted HMSC model (Figure 2)."
    )
    figure(doc, source_image(source, 1), 6.15,
           "Figure 2. Plant200 input coverage and represented lineages. The "
           "modern records, trait availability and tree-time counts describe "
           "the selected surviving sample, not the complete historical plant "
           "biota.")
    prose(doc,
        "Occurrence zeros are target-group non-records, not standardised survey "
        "absences. The probit HMSC fit used eight preprocessed environmental "
        "axes: annual temperature and precipitation, elevation, temperature "
        "and precipitation seasonality, moisture availability, wetland potential "
        "and local relief. Standardised log target-group record count entered "
        "only as an observation-effort covariate. Hmsc::Hmsc() was run with "
        "XScale = FALSE because environmental columns were already transformed. "
        "The analysed fit did not include measured traits, phylogenetic "
        "covariance or a spatial random effect. Four chains retained 1,000 "
        "samples each; the archived coefficient check covers only 32 terms."
    )

    heading(doc, "Ancestral projection and geographic scenario", 2)
    prose(doc,
        "Five selected HMSC posterior draws were carried into the illustrated "
        "ancestral-response reconstruction. For each draw, coefficient-wise "
        "node estimates were fitted on the dated tree and branch-time values "
        "were approximated. The modern environmental transformation was held "
        "fixed when the reconstructed responses were evaluated against "
        "palaeoenvironmental cells. Time-specific environmental support was "
        "summarised across active lineages on a common zero-to-one scale. "
        "A modern prediction surface was never averaged backwards to form "
        "an ancestral map."
    )
    prose(doc,
        "The separate geographic calculation used 4,050 global four-degree "
        "states, 66 mapped times from 325 to 0 Ma and a 255-node master "
        "timeline that includes dated-tree events. One-degree PALEOMAP carrier "
        "tables supplied palaeocoordinates and land fractions, aggregated "
        "onto the four-degree lattice. The movement kernel was run on an "
        "eight-neighbour graph with a fixed low marine-transit weight. "
        "Recorded modern cells supplied terminal location fields; root "
        "weights and time-ordered transitions supplied ancestral location "
        "support. The HMSC response draws did not modify this geographic "
        "kernel or its likelihood. Accordingly, variation among the five "
        "response draws does not quantify geographic uncertainty."
    )
    prose(doc,
        "The geographic output is a source-normalised location field, "
        "conditioned on the supplied modern endpoint. A separate terminal "
        "calibration produces a lineage-support intensity useful for displaying "
        "the calculation, but it does not turn location density into the "
        "probability that a population occupied each cell. This run had no "
        "independently fitted arrival, establishment or local-retention model, "
        "no independent fossil likelihood, and no fitted biological interaction "
        "network. It therefore cannot produce empirical historical occupancy "
        "richness or historical extinction rates."
    )

    heading(doc, "Results")
    heading(doc, "Modern predictive performance", 2)
    prose(doc,
        "In the 122 held-out four-degree cells, observed recorded prevalence "
        "was 0.01951 and mean predicted prevalence was 0.02083. The aggregate "
        "Brier score was 0.01811 and mean log score was -0.07612. Observed "
        "recorded-species counts averaged 3.902 per cell, compared with "
        "4.167 summed predictions; their correlation was 0.8151 (Figure 3). "
        "Calibration intercept and slope were -0.750 and 0.743, so similarity "
        "of the two prevalence means did not establish perfect calibration. "
        "The held-out selection retained at least one occupied training cell "
        "per species, weakening strict spatial-block independence. The saved "
        "prediction diagnostic averaged 12 selected posterior draws, not the "
        "full retained posterior."
    )
    prose(doc,
        "For the 32 archived coefficient diagnostics, Rhat ranged from "
        "1.00017 to 1.03332 and effective sample size from 247.4 to "
        "2,767.3. Some checked values therefore exceed a 1.01 Rhat "
        "reference; the compact check cannot establish convergence of every "
        "parameter. These are diagnostics for the recorded occurrence process, "
        "not independent validation of ancient ranges."
    )
    figure(doc, source_image(source, 2), 6.15,
           "Figure 3. Spatial holdout predictions and compact MCMC diagnostics "
           "for the modern Plant200 HMSC fit. The diagnostic concerns recorded "
           "occurrence under the stated target-group observation design.")

    heading(doc, "Ancestral environmental responses", 2)
    prose(doc,
        "The active-lineage count was one at the root, two at 300 and 200 Ma, "
        "38 at 100 Ma and 200 at the present. The ancestral environmental "
        "support surfaces varied among historical land cells even when only "
        "one or two sampled branches were active (Figure 4). The pattern "
        "illustrates the different roles of the response and landscape "
        "inputs. It does not demonstrate that the lineage reached every "
        "environmentally suitable cell. The 0 Ma environmental-support panel "
        "is likewise a conditional response over land, not a map of GBIF "
        "records or detection-corrected occupancy."
    )
    figure(doc, source_image(source, 3), 6.15,
           "Figure 4. Ancestral environmental support at selected times. "
           "Time-matched palaeoenvironmental grids were evaluated using "
           "reconstructed active-lineage response coefficients. Each map "
           "shows environmental support, not realised historical occupancy.")

    heading(doc, "Geographic output and terrestrial state audit", 2)
    prose(doc,
        "All 462 indexed PNG-GeoTIFF pairs for seven mapped variables and 66 "
        "times were present, and the largest archived transition-column "
        "normalisation error was approximately 6.66e-16. These are "
        "software and numerical checks, not evidence of terrestrial "
        "predictive validity. Reaggregation also corrected an artificial "
        "root-time Shannon value: for one contributing lineage, the revised "
        "mixture Shannon is zero. That correction changes the displayed "
        "support diagnostic but not the stored location fields."
    )
    prose(doc,
        "The stronger audit concerns geography. Approximately 100.00% of "
        "normalised location mass lay in marine-classified states at 300 Ma, "
        "99.90% at 200 Ma, 92.94% at 100 Ma and 0.92% at 0 Ma. Fifty-nine "
        "of 66 mapped times had more than half of location mass in marine "
        "states (Figure 5). Three recorded terminal cells were also classified "
        "marine in the carrier alignment. For terrestrial plants, the "
        "deep-time result fails a basic state-domain check. Blanking oceans "
        "in a PNG makes a cleaner terrestrial figure but does not change "
        "the underlying transition calculation. The historical geographic "
        "maps therefore cannot be reported as inferred terrestrial "
        "palaeodistributions or biodiversity patterns."
    )
    figure(doc, marine_path, 6.15,
           "Figure 5. Marine-state audit of the Plant200 geographic scenario. "
           "The fraction of normalised location mass assigned to marine "
           "carrier states exceeds 50% for 59 of 66 mapped times. This is "
           "a model diagnostic failure for terrestrial palaeodistribution, "
           "not evidence of ancestral marine plants. The series are drawn "
           "from the archived 66-time land-mask audit.")

    caption(doc, "Table 2. Key data and diagnostic values in the audited application.")
    add_table(doc,
        ["Quantity", "Audited value", "Interpretation"],
        [
            ["Taxa and tree", "200 species; 68 families; root about 325 Ma", "Selected extant-surviving lineages only"],
            ["Modern HMSC design", "432 training; 122 held-out four-degree cells", "Occurrence non-records, not survey absences"],
            ["Historical grid", "4,050 four-degree states; 66 mapped times", "Includes marine transit states"],
            ["Brier / count correlation", "0.01811 / 0.8151", "Modern recorded-occurrence diagnostic"],
            ["Maps verified", "462 PNG and 462 GeoTIFF", "File existence and indexing, not model validity"],
            ["Marine-dominated times", "59 of 66", "Terrestrial historical inference fails"],
            ["Trait coverage", "59.0% to 97.5%", "Traits audited, not fitted in analysed HMSC"],
        ],
        [1.82, 2.02, 2.81])

    heading(doc, "Discussion")
    heading(doc, "A configurable framework is not a claim of complete inference", 2)
    prose(doc,
        "The useful lesson from gen3sis is methodological organisation: "
        "landscapes, initial conditions, biological functions and observable "
        "outputs should be specified separately and compared under known "
        "configurations (Hagen et al., 2021). HmscEE applies that organisation "
        "to data-linked historical analysis. Its six named functions state "
        "who is responsible for an effect, but availability of a function "
        "does not identify its parameters. Environmental filtering and "
        "tree-conditioned response evolution are exercised in Plant200; "
        "geographic movement is a fixed scenario; biotic filtering and "
        "global lineage extinction are not estimated. Dated-tree "
        "speciation controls identity but does not discover all historical "
        "species. The functions are connected through their inputs and "
        "state transitions rather than multiplied as unrelated scores."
    )
    prose(doc,
        "This distinction also helps interpret the apparently different "
        "map products. Ancestral environmental support can be high on land "
        "that the lineage never reached. Location support integrates a "
        "movement assumption and modern endpoint but remains normalised "
        "over space. A calibrated occupancy model would additionally need "
        "observations or defensible priors for establishment and loss, with "
        "a clear observation process. Without those pieces, support-mixture "
        "Shannon and effective support counts describe the location-field "
        "mixture, not ancient community Shannon diversity. The marine-state "
        "failure shows why this semantic separation is a necessary part of "
        "model checking, not merely careful wording."
    )

    heading(doc, "What the worked example does not establish", 2)
    prose(doc,
        "Several limits are substantive. The 200 species span 68 families "
        "but are not a dense, complete clade; extinct tips are absent. "
        "The modern response matrix derives from occurrence records with "
        "imperfect observation. The five selected response draws, a single "
        "tree and carrier reconstruction, fixed movement settings and "
        "coefficient-wise ancestral approximation understate historical "
        "uncertainty. Response draws were not coupled to the geographic "
        "transition model. The recent endpoint is imposed from records that "
        "also contributed to the HMSC fit; reproducing it is therefore "
        "conditioning, not an independent prediction. Family-level fossils "
        "could assess clade presence over age and location intervals, but "
        "they did not enter a fitted likelihood in this run and cannot "
        "validate a named modern species at deep time."
    )
    prose(doc,
        "A future terrestrial analysis needs the movement graph and "
        "plate-carriage mapping to preserve land-state support through "
        "time, including explicit checks for coastal and island cells. "
        "It must then pass grid and time-step refinement checks and be "
        "tested against simulations with known ancestral states. "
        "Independent modern holdouts or other endpoint data are needed "
        "before terminal fit can be presented as predictive validation. "
        "Ancient occupancy, phylogenetic diversity, functional diversity, "
        "refugia and turnover should be reported only after joint state "
        "draws are generated and externally checked."
    )

    heading(doc, "Conclusions")
    prose(doc,
        "HmscEE organises a deep-time community workflow around six "
        "biological function families while keeping dynamic Earth inputs, "
        "historical constraints and diagnostic outputs separate. Its "
        "current Plant200 application supports a modern HMSC response "
        "analysis, approximate ancestral environmental projection and "
        "an auditable geographic scenario. It does not yet support "
        "validated terrestrial occupancy histories or empirical "
        "ancient-diversity maps. The explicit six-function design makes "
        "that boundary visible and provides a tractable basis for "
        "subsequent, better-constrained tests."
    )

    heading(doc, "Data and code availability")
    prose(doc,
        "The analysed software is HmscEcoEvo version 1.0.3. Local "
        "reproducibility materials include source scripts, data inventories, "
        "figure-source tables, indexed raster paths and audit summaries. "
        "A permanent public repository and accession number have not yet "
        "been established. Original occurrence, trait, tree and "
        "palaeoenvironment providers' licences and identifiers must be "
        "retained when the reproducibility bundle is deposited."
    )
    prose(doc,
        "Release-alignment note. The archived 1.0.3 package catalog still "
        "reports colonisation and persistence as two legacy process labels. "
        "The six-function organisation in this manuscript is the revised "
        "conceptual configuration, and the public catalog and documentation "
        "need alignment before submission or a new package release. The "
        "legacy helper exports should not be mistaken for empirical "
        "modules activated in the Plant200 application."
    )

    heading(doc, "Author declarations")
    prose(doc,
        "Author names, affiliations, funding, contributions and competing "
        "interest statements require author confirmation. AI assistance "
        "was used for manuscript restructuring and language drafting; "
        "the authors remain responsible for data provenance, scientific "
        "interpretation, code and reference verification before submission."
    )

    heading(doc, "References")
    refs = [
        "Arias, J. S. (2024). Phylogenetic biogeography inference using dynamic paleogeography models and explicit geographic ranges. Systematic Biology, 73, 995-1014. https://doi.org/10.1093/sysbio/syae051",
        "Blanchet, F. G., Cazelles, K., & Gravel, D. (2020). Co-occurrence is not evidence of ecological interactions. Ecology Letters, 23, 1050-1063. https://doi.org/10.1111/ele.13525",
        "Flannery-Sutherland, J. T., Elsler, A., Farnsworth, A., Lunt, D. J., & Benton, M. J. (2025). Landscape-explicit phylogeography illuminates the ecographic radiation of early archosauromorph reptiles. Nature Ecology & Evolution, 9, 1138-1152. https://doi.org/10.1038/s41559-025-02739-y",
        "Guillory, W. X., & Brown, J. L. (2021). A new method for integrating ecological niche modeling with phylogenetics to estimate ancestral distributions. Systematic Biology, 70, 1033-1045. https://doi.org/10.1093/sysbio/syab016",
        "Hagen, O., Fl\u00fcck, B., Fopp, F., Cabral, J. S., Hartig, F., Pontarp, M., Rangel, T. F., & Pellissier, L. (2021). gen3sis: A general engine for eco-evolutionary simulations of the processes that shape Earth's biodiversity. PLOS Biology, 19, e3001340. https://doi.org/10.1371/journal.pbio.3001340",
        "Jin, Y., & Qian, H. (2022). V.PhyloMaker2: An updated and enlarged R package that can generate very large phylogenies for vascular plants. Plant Diversity, 44, 335-339. https://doi.org/10.1016/j.pld.2022.05.005",
        "Ovaskainen, O., Tikhonov, G., Norberg, A., Blanchet, F. G., Duan, L., Dunson, D., Roslin, T., & Abrego, N. (2017). How to make more out of community data? A conceptual framework and its implementation as models and software. Ecology Letters, 20, 561-576. https://doi.org/10.1111/ele.12757",
        "Tikhonov, G., Opedal, \u00d8. H., Abrego, N., Lehikoinen, A., de Jonge, M. M. J., Oksanen, J., & Ovaskainen, O. (2020). Joint species distribution modelling with the R-package Hmsc. Methods in Ecology and Evolution, 11, 442-447. https://doi.org/10.1111/2041-210X.13345",
    ]
    for ref in refs:
        p = doc.add_paragraph(ref)
        p.paragraph_format.left_indent = Inches(0.23)
        p.paragraph_format.first_line_indent = Inches(-0.23)
        p.paragraph_format.space_after = Pt(5)
        for run in p.runs:
            run.font.size = Pt(9.4)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    doc.save(OUTPUT)
    print(OUTPUT)
    print("Paragraphs:", len(doc.paragraphs), "Tables:", len(doc.tables),
          "Figures:", len(doc.inline_shapes))


if __name__ == "__main__":
    build()
