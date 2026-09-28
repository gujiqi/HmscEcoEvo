#' Write a complete end-to-end workflow Word report
#'
#' Builds a report from an existing Case 1 or Case 2 output directory by
#' scanning all generated figures and tables, assigning them to workflow
#' sections, and rendering a Word document. This avoids the fragile earlier
#' pattern where one representative figure was manually paired with each
#' chapter and some outputs could silently disappear from the report.
#'
#' @param output_root Output directory containing `figures/`, `tables/`, and
#'   `reports/`.
#' @param output_file Destination `.docx`. Defaults to
#'   `reports/hmscecoevo_full_540Ma_workflow.docx` under `output_root`.
#' @param include_all_figures Logical; if `TRUE`, every PNG in `figures/` is
#'   inserted into a module section or into the appendix.
#' @param max_table_rows Number of rows shown for each previewed CSV table.
#' @return Normalized path to the generated Word document.
#' @export
#'
#' @examples
#' \dontrun{
#' hee_write_full_workflow_report("outputs/case01_simulated_full_540Ma")
#' }
hee_write_full_workflow_report <- function(output_root,
                                           output_file = NULL,
                                           include_all_figures = TRUE,
                                           max_table_rows = 8) {
  .require_pkg("rmarkdown", "complete Word report rendering")
  output_root <- normalizePath(output_root, winslash = "/", mustWork = TRUE)
  fig_dir <- file.path(output_root, "figures")
  tab_dir <- file.path(output_root, "tables")
  rep_dir <- file.path(output_root, "reports")
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(output_file)) {
    output_file <- file.path(rep_dir, "hmscecoevo_full_540Ma_workflow.docx")
  }
  output_file <- normalizePath(output_file, winslash = "/", mustWork = FALSE)

  figs <- if (dir.exists(fig_dir)) {
    sort(list.files(fig_dir, pattern = "\\.png$", full.names = FALSE))
  } else character()
  tabs <- if (dir.exists(tab_dir)) {
    sort(list.files(tab_dir, pattern = "\\.csv$", full.names = FALSE))
  } else character()

  section_specs <- list(
    list(id = "1", title = "Executive Summary And Run Inventory",
         fig = c("^output_", "^model_agreement_map", "^uncertainty_budget"),
         tab = c("^case01_run_summary", "^case02_run_summary",
                 "^run_metadata", "^output_file_inventory"),
         note = "This section gives the run-level audit trail: output root, time-axis coverage, figure/table counts, MCMC mode and output inventory."),
    list(id = "2", title = "HMSC Four Input Matrices And Data Preparation",
         fig = c("^0[1-9]_", "^10_", "^11_", "^12_"),
         tab = c("^Y_matrix", "^XData_", "^TrData", "^phylogenetic_correlation",
                 "^tip_ranges", "^modern_sites", "^plate_points",
                 "^pseudo_fossils"),
         note = "These panels document Y, XData, TrData, the dated tree/C matrix, tip ranges, 0 Ma environments and modern site inputs."),
    list(id = "3", title = "HMSC MCMC, Posterior Summaries And Modern Validation",
         fig = c("^1[3-9]_", "^20_", "^21_", "^22"),
         tab = c("^hmsc_", "^cv_", "^hmsc_species_performance"),
         note = "These outputs are produced from a real HMSC model and MCMC run. Short local runs are workflow checks; inference needs full convergence."),
    list(id = "4", title = "BioGeoBEARS And Historical Accessibility",
         fig = c("^23", "^24", "^25", "^26", "^27", "^28"),
         tab = c("^bgb_"),
         note = "Case 1 BioGeoBEARS-like objects are labelled mock/pseudo workflow data. Case 2 must use externally supplied BioGeoBEARS outputs."),
    list(id = "5", title = "Phylogenetic Time Mask And Plate Correction",
         fig = c("^phylo_mask", "^valid_species", "^plate_"),
         tab = c("^phylo_time_mask", "^valid_species", "^plate_"),
         note = "Species are masked before their origin time. Plate correction uses supplied reconstructed coordinates or user tables; the package does not infer plate motion from present coordinates alone."),
    list(id = "6", title = "Deep-Time Projection And Representative Maps",
         fig = c("^deep_", "^projection_"),
         tab = c("^projection_", "^M5_dynamic_dispersal_representative_maps"),
         note = "These maps cover representative geological times. The script now keeps both geological event times and the older 60 Myr reference times."),
    list(id = "7", title = "M1-M5 Projection Models And Dispersal Contrasts",
         fig = c("^M1_M5", "^M3_M5", "^M5_", "^dispersal_",
                 "^model_agreement", "^projection_component"),
         tab = c("^M1_M5", "^M3_M5", "^M5_", "^projection_model",
                 "^projection_component", "^M5_minus", "^dynamic_diag"),
         note = "M1-M5 are projection scenarios, not BioGeoBEARS candidate models. M3-M5 isolate no/static/dynamic dispersal effects."),
    list(id = "8", title = "Geological Processes And Dynamic Earth-Biota Assembly",
         fig = c("^geoprocess_"),
         tab = c("^geoprocess_"),
         note = "These diagnostics cover arena change, opportunity, connectivity, source pressure, rescue, dynamic colonisation, extinction layers, species pools, M0-M8, H1-H8, limitation mechanisms and reliability."),
    list(id = "9", title = "Community Metrics, Turnover, Refugia, CWM, PD And Composition",
         fig = c("^turnover", "^refugia", "^hotspot", "^CWM",
                 "^phylogenetic_diversity", "^composition"),
         tab = c("^turnover", "^refugia", "^CWM", "^phylogenetic",
                 "^composition"),
         note = "Community maps are derived from model probabilities and should be read as expected summaries, not observed census counts."),
    list(id = "10", title = "Fossil, Pseudo-Hindcast, Simulation Recovery And Uncertainty",
         fig = c("^fossil", "^pseudo", "^simulation", "^coverage",
                 "^bias", "^uncertainty"),
         tab = c("^fossil", "^simulation", "^uncertainty"),
         note = "Case 1 validation is simulated or pseudo-validation. Real fossil validation requires independent fossil records and spatial/temporal matching."),
    list(id = "11", title = "Audit, Completeness And Reproducibility",
         fig = character(),
         tab = c("^output_file_inventory", "^case01_run_summary",
                 "^case02_run_summary", "^run_metadata"),
         note = "The completeness index below records which generated files were inserted or previewed in this report.")
  )

  matches_any <- function(x, patterns) {
    if (length(patterns) == 0L || length(x) == 0L) return(rep(FALSE, length(x)))
    Reduce(`|`, lapply(patterns, function(p) grepl(p, x)))
  }

  used_figs <- character()
  used_tabs <- character()
  section_assets <- lapply(section_specs, function(spec) {
    sf <- figs[matches_any(figs, spec$fig)]
    st <- tabs[matches_any(tabs, spec$tab)]
    used_figs <<- unique(c(used_figs, sf))
    used_tabs <<- unique(c(used_tabs, st))
    list(figs = sf, tabs = st)
  })
  remaining_figs <- setdiff(figs, used_figs)
  remaining_tabs <- setdiff(tabs, used_tabs)

  idx <- data.frame(
    file = c(file.path("figures", figs), file.path("tables", tabs)),
    type = c(rep("figure", length(figs)), rep("table", length(tabs))),
    included_in_report = c(figs %in% used_figs | figs %in% remaining_figs,
                           rep(TRUE, length(tabs))),
    section = NA_character_,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(section_specs)) {
    spec <- section_specs[[i]]
    assets <- section_assets[[i]]
    idx$section[idx$file %in% file.path("figures", assets$figs)] <-
      paste(spec$id, spec$title)
    idx$section[idx$file %in% file.path("tables", assets$tabs)] <-
      paste(spec$id, spec$title)
  }
  idx$section[idx$file %in% file.path("figures", remaining_figs)] <-
    "12 Complete Figure Appendix"
  idx$section[idx$file %in% file.path("tables", remaining_tabs)] <-
    "11 Audit, Completeness And Reproducibility"
  idx_path <- file.path(tab_dir, "full_workflow_report_completeness_index.csv")
  if (dir.exists(tab_dir)) {
    utils::write.csv(idx, idx_path, row.names = FALSE)
  }

  fig_rel <- function(name) file.path("..", "figures", name)
  tab_rel <- function(name) file.path("..", "tables", name)
  safe_chunk <- function(prefix, i, j = NULL) {
    paste0(prefix, "_", i, if (!is.null(j)) paste0("_", j) else "")
  }
  add_table_preview <- function(lines, tab, chunk_name) {
    c(lines,
      paste0("### Table Preview: `", tab, "`"),
      "",
      "```{r " %+% chunk_name %+% ", echo=FALSE}",
      "p <- '" %+% tab_rel(tab) %+% "'",
      "if (file.exists(p)) {",
      "  x <- utils::read.csv(p, check.names = FALSE)",
      "  knitr::kable(utils::head(x, " %+% max_table_rows %+% "))",
      "} else {",
      "  cat('Missing table:', p)",
      "}",
      "```",
      "")
  }

  `%+%` <- function(a, b) paste0(a, b)
  lines <- c(
    "---",
    "title: \"HmscEcoEvo Complete 540 Ma Workflow Report\"",
    "output:",
    "  word_document:",
    "    toc: true",
    "    toc_depth: 3",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
    "```",
    "",
    "# 0. How To Read This Report",
    "",
    "This report is generated from the output directory, not from a hand-picked figure list. It keeps the previous workflow outputs and the newer geological-process outputs together. Figures are grouped by file-name semantics; any unclassified figure is retained in the complete appendix so it is not silently lost.",
    "",
    paste0("- Output root: `", output_root, "`"),
    paste0("- Figures detected: ", length(figs)),
    paste0("- CSV tables detected: ", length(tabs)),
    paste0("- Completeness index: `tables/full_workflow_report_completeness_index.csv`"),
    ""
  )

  for (i in seq_along(section_specs)) {
    spec <- section_specs[[i]]
    assets <- section_assets[[i]]
    lines <- c(lines, paste0("# ", spec$id, ". ", spec$title), "",
               spec$note, "")
    if (length(assets$figs) > 0L) {
      lines <- c(lines, paste0("## ", spec$id, ".1 Figures"), "")
      for (j in seq_along(assets$figs)) {
        fig <- assets$figs[[j]]
        lines <- c(lines,
                   paste0("### ", spec$id, ".1.", j, " `", fig, "`"),
                   "",
                   paste0("![](", fig_rel(fig), "){width=100%}"),
                   "")
      }
    }
    if (length(assets$tabs) > 0L) {
      lines <- c(lines, paste0("## ", spec$id, ".2 Table Previews"), "")
      show_tabs <- assets$tabs[seq_len(min(length(assets$tabs), 12L))]
      for (j in seq_along(show_tabs)) {
        lines <- add_table_preview(lines, show_tabs[[j]],
                                   safe_chunk("tab", i, j))
      }
      if (length(assets$tabs) > length(show_tabs)) {
        lines <- c(lines, paste0("Additional tables in this section: ",
                                 paste(sprintf("`%s`", assets$tabs[-seq_along(show_tabs)]),
                                       collapse = ", ")), "")
      }
    }
  }

  if (include_all_figures && length(remaining_figs) > 0L) {
    lines <- c(lines, "# 12. Complete Figure Appendix", "",
               "The following figures were not matched by the primary section rules, so they are still included here to prevent silent omissions.", "")
    for (j in seq_along(remaining_figs)) {
      fig <- remaining_figs[[j]]
      lines <- c(lines,
                 paste0("## 12.", j, " `", fig, "`"),
                 "",
                 paste0("![](", fig_rel(fig), "){width=100%}"),
                 "")
    }
  }

  if (length(remaining_tabs) > 0L) {
    lines <- c(lines, "# 13. Complete Table Inventory", "",
               "These CSV tables are not previewed above but are present in the output directory and listed in the completeness index.", "")
    lines <- add_table_preview(lines,
                               "full_workflow_report_completeness_index.csv",
                               "complete_index")
  }

  lines <- c(lines,
             "# 14. Scientific Boundaries",
             "",
             "- Deep-time species-level maps older than 20 Ma should be read cautiously; older than 50 Ma they are primarily lineage/clade-level potential-suitability hypotheses.",
             "- Modern HMSC spatial random effects and residual associations are not projected into deep time by default.",
             "- BioGeoBEARS, plate reconstruction, fossil calibration and true speciation-rate estimation require external specialist data or models.",
             "- Geological-process opportunity, extinction-layer, cradle/museum/grave and limitation labels are transparent proxy diagnostics unless independently calibrated.",
             "- Colonisation, extinction, rescue, connectivity, speciation-opportunity, reliability, and M1-M5 dispersal contrasts are heuristic or proxy scenario diagnostics unless fitted to independent process data.",
             "- These outputs are non-causal summaries of the supplied model components; they should not be read as proof of realised dispersal routes, true palaeo-distributions, or actual speciation/extinction rates.",
             "- Prediction reliability is an epistemic diagnostic, not a posterior credible interval.",
             "",
             "# 15. Session Info",
             "",
             "```{r session_info}",
             "sessionInfo()",
             "```",
             "")

  rmd <- file.path(rep_dir, "hmscecoevo_full_540Ma_workflow.Rmd")
  writeLines(lines, rmd)
  hee_render_word_report(rmd, output_file)
}
