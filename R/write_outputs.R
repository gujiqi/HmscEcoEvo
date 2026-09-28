#' Build an ordered output catalog
#'
#' Create a human-readable catalog for workflow outputs. The catalog assigns
#' hierarchical numbers such as `1.1.1`, `2.3.2`, and `4.5.09` to figures,
#' tables, reports, task tables, metadata, and RDS files so outputs can be read
#' in workflow order rather than filesystem order.
#'
#' @param output_root Case output directory containing `figures`, `tables`,
#'   `reports`, `metadata`, `tasks`, and optionally `rds`.
#' @param include_rds Logical. Include RDS objects in the catalog.
#'
#' @return A data frame with original and ordered output paths.
#' @export
#'
#' @examples
#' \dontrun{
#' catalog <- hee_ordered_output_catalog("outputs/case01_simulated_full_540Ma")
#' }
hee_ordered_output_catalog <- function(output_root, include_rds = TRUE) {
  output_root <- normalizePath(output_root, winslash = "/", mustWork = TRUE)
  files <- list.files(output_root, recursive = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir == FALSE]
  rel <- normalizePath(files, winslash = "/", mustWork = TRUE)
  rel <- substring(rel, nchar(output_root) + 2L)
  rel <- rel[!grepl("^ordered_results/", rel)]
  rel <- rel[!grepl("^reports/(cn_explained|.*rendered_pages|lo_pdf_test|_cn_explanation_image_cache)", rel)]
  if (!include_rds) rel <- rel[!grepl("^rds/", rel)]
  if (length(rel) == 0L) {
    return(data.frame(
      order_number = character(), section = character(), section_title = character(),
      item_title = character(), file_type = character(), original_relpath = character(),
      ordered_relpath = character(), stringsAsFactors = FALSE
    ))
  }

  section_dirs <- c(
    "1" = "1_input_data_four_matrices",
    "2" = "2_hmsc_model_mcmc_validation",
    "3" = "3_biogeography_phylogeny_plate",
    "4" = "4_deep_time_projection",
    "5" = "5_M1_M5_model_comparison",
    "6" = "6_community_results",
    "7" = "7_validation_uncertainty",
    "8" = "8_reports_metadata_tasks",
    "9" = "9_rds_objects"
  )
  section_titles <- c(
    "1" = "Data inputs and HMSC four matrices",
    "2" = "HMSC model, posterior, MCMC and modern validation",
    "3" = "Historical biogeography, phylogenetic time mask and plate correction",
    "4" = "Deep-time projection maps from 540 Ma to 0 Ma",
    "5" = "Five projection model comparison",
    "6" = "Community-level derived results",
    "7" = "Validation and uncertainty",
    "8" = "Reports, metadata, task tables and inventories",
    "9" = "RDS objects for reproducibility"
  )

  exact <- data.frame(
    basename = c(
      "01_Y_heatmap.png", "02_species_prevalence.png", "03_site_richness.png",
      "03b_zero_proportion.png", "Y_matrix.csv",
      "04_X_correlation.png", "05_X_environment_maps_0Ma.png",
      "05b_X_histograms.png", "05c_0Ma_environment_layer_map.png",
      "XData_0Ma_from_env_cube_raw.csv", "XData_0Ma_from_env_cube.csv",
      "06_traits_heatmap.png", "07_traits_pca.png", "07b_trait_missingness.png",
      "TrData_traits.csv",
      "08_phylo_tree.png", "08_phylo_tree.pdf", "09_phylo_correlation_matrix.png",
      "09b_species_origin_time.png", "10_tip_ranges_heatmap.png",
      "11_region_map_0Ma.png", "12_landmask_map_0Ma.png",
      "phylogenetic_correlation_C.csv", "tip_ranges.csv", "modern_sites.csv",
      "plate_points.csv", "pseudo_fossils.csv",
      "13_hmsc_beta_heatmap.png", "hmsc_beta_posterior_mean.csv",
      "14_hmsc_gamma_heatmap.png", "hmsc_gamma_posterior_mean.csv",
      "15_hmsc_rho_posterior.png", "hmsc_rho_posterior.csv",
      "16_hmsc_omega_heatmap.png", "hmsc_omega_posterior_mean.csv",
      "17_hmsc_variance_partitioning.png", "hmsc_model_fit_metrics.csv",
      "18_hmsc_traceplots.png", "19_hmsc_ess_histogram.png", "hmsc_ess.csv",
      "20_hmsc_rhat_histogram.png", "hmsc_rhat.csv",
      "21_hmsc_observed_vs_predicted.png", "22_hmsc_cv_performance.png",
      "22b_spatial_block_CV_map.png", "22c_richness_observed_vs_predicted.png",
      "22d_calibration_plot.png", "22e_CV_performance_heatmap.png",
      "22f_species_level_performance.png", "22g_environment_block_CV_plot.png",
      "cv_folds.csv", "cv_performance.csv", "hmsc_species_performance.csv",
      "23_bgb_model_weights.png", "23b_bgb_parameter_plot.png",
      "24_bgb_lrt_table.png", "25_bgb_ancestral_range_heatmap.png",
      "26_bgb_bsm_events_time.png", "27_bgb_source_sink_heatmap.png",
      "28_bgb_dispersal_network.png", "bgb_model_compare.csv", "bgb_lrt.csv",
      "bgb_bsm_events.csv", "bgb_accessibility.csv",
      "phylo_mask_heatmap.png", "valid_species_through_time.png",
      "phylo_time_mask.csv", "valid_species_through_time.csv",
      "plate_tectonic_tracks_through_time.png",
      "plate_corrected_vs_uncorrected_environment_difference.png",
      "plate_corrected_vs_uncorrected.csv",
      "plate_corrected_uncorrected_environment_difference.csv",
      "projection_task_table.csv", "extrapolation_summary_by_time.csv",
      "projection_model_definitions.csv",
      "M1_M5_richness_through_time.png", "M1_M5_richness_through_time.csv",
      "M1_M5_area_through_time.png", "M1_M5_occupied_area_through_time.csv",
      "M1_M5_centroid_shift.png", "M1_M5_centroid_shift.csv",
      "M5_minus_M3_richness_difference.png",
      "M5_minus_M3_richness_difference_representative.csv",
      "dispersal_vs_no_dispersal_difference.png", "model_agreement_map.png",
      "turnover_vs_present_maps.png", "turnover_vs_present_representative.csv",
      "refugia_score_map.png", "refugia_score.csv", "hotspot_persistence_map.png",
      "CWM_trait_maps.png", "phylogenetic_diversity_maps.png", "composition_pca_maps.png",
      "fossil_validation_map.png", "fossil_validation.csv",
      "fossil_probability_histogram.png", "pseudo_hindcast_validation.png",
      "simulation_recovery_richness.png", "simulation_recovery_probability.png",
      "simulation_recovery.csv", "coverage_plot.png", "bias_plot.png",
      "uncertainty_budget.png", "uncertainty_budget.csv",
      "hmscecoevo_full_540Ma_workflow.docx", "hmscecoevo_full_540Ma_workflow.Rmd",
      paste0("hmscecoevo_full_540Ma_workflow_", "\u4e2d\u6587\u8be6\u89e3", ".docx"),
      "hmscecoevo_full_540Ma_workflow_cn_explained.docx",
      "output_file_inventory.csv", "run_metadata.csv", "case01_run_summary.csv",
      "case02_run_summary.csv", "ordered_output_index.csv"
    ),
    order_number = c(
      "1.1.1", "1.1.2", "1.1.3", "1.1.4", "1.1.5",
      "1.2.1", "1.2.2", "1.2.3", "1.2.4", "1.2.5", "1.2.6",
      "1.3.1", "1.3.2", "1.3.3", "1.3.4",
      "1.4.1", "1.4.1", "1.4.2", "1.4.3", "1.4.4",
      "1.4.5", "1.4.6", "1.4.7", "1.4.8", "1.4.9",
      "1.4.10", "1.4.11",
      "2.1.1", "2.1.1", "2.1.2", "2.1.2", "2.1.3", "2.1.3",
      "2.1.4", "2.1.4", "2.1.5", "2.1.6",
      "2.2.1", "2.2.2", "2.2.2", "2.2.3", "2.2.3",
      "2.3.1", "2.3.2", "2.3.3", "2.3.4", "2.3.5", "2.3.6",
      "2.3.7", "2.3.8", "2.3.9", "2.3.10", "2.3.11",
      "3.1.1", "3.1.2", "3.1.3", "3.1.4", "3.1.5", "3.1.6",
      "3.1.7", "3.1.8", "3.1.9", "3.1.10", "3.1.11",
      "3.2.1", "3.2.2", "3.2.1", "3.2.2",
      "3.3.1", "3.3.2", "3.3.3", "3.3.4",
      "4.8.1", "4.8.2", "5.0.1",
      "5.1.1", "5.1.1", "5.1.2", "5.1.2", "5.1.3", "5.1.3",
      "5.2.1", "5.2.1", "5.2.2", "5.2.3",
      "6.1.1", "6.1.1", "6.2.1", "6.2.1", "6.3.1",
      "6.4.1", "6.5.1", "6.6.1",
      "7.1.1", "7.1.1", "7.1.2", "7.2.1",
      "7.3.1", "7.3.2", "7.3.3", "7.4.1", "7.4.2",
      "7.5.1", "7.5.1",
      "8.1.1", "8.1.2", "8.1.3", "8.1.4",
      "8.2.1", "8.2.2", "8.2.3", "8.2.4", "8.2.5"
    ),
    item_title = c(
      "Y community matrix heatmap", "Species prevalence", "Site richness",
      "Zero proportion", "Y community matrix table",
      "XData correlation", "0 Ma environmental maps", "XData histograms",
      "0 Ma environment layer map", "Raw XData from 0 Ma cube",
      "Recipe-scaled XData for HMSC",
      "Traits heatmap", "Traits PCA", "Trait missingness", "Trait table",
      "Time-calibrated phylogenetic tree", "Time-calibrated phylogenetic tree",
      "Phylogenetic correlation matrix heatmap", "Species origin time",
      "Tip ranges heatmap", "Region map at 0 Ma", "Landmask map at 0 Ma",
      "Phylogenetic correlation matrix table", "Tip ranges table",
      "Modern sites table", "Plate points table", "Pseudo-fossils table",
      "Beta posterior mean", "Beta posterior mean", "Gamma posterior mean",
      "Gamma posterior mean", "Rho posterior", "Rho posterior",
      "Omega residual association", "Omega residual association",
      "Variance partitioning", "HMSC model fit metrics",
      "MCMC traceplots", "ESS histogram", "ESS table", "Rhat histogram",
      "Rhat table", "Observed vs predicted", "CV performance",
      "Spatial block CV map", "Richness observed vs predicted",
      "Calibration plot", "CV performance heatmap",
      "Species-level performance", "Environment block CV plot",
      "CV fold assignments", "CV performance table", "Species performance table",
      "BGB model weights", "BGB parameter plot", "BGB LRT table plot",
      "Ancestral range probability heatmap", "BSM events through time",
      "BGB source sink heatmap", "BGB dispersal network",
      "BGB model comparison table", "BGB LRT table", "BGB BSM events table",
      "BGB accessibility table", "Phylogenetic time mask heatmap",
      "Valid species through time", "Phylogenetic time mask table",
      "Valid species through time table", "Tectonic tracks through time",
      "Plate-corrected environment difference", "Plate corrected vs uncorrected table",
      "Plate environment difference table", "Projection task table",
      "Extrapolation summary by time", "Projection model definitions",
      "M1-M5 richness through time", "M1-M5 richness through time",
      "M1-M5 occupied area through time", "M1-M5 occupied area through time",
      "M1-M5 centroid shift", "M1-M5 centroid shift",
      "M5 minus M3 richness difference", "M5 minus M3 richness difference",
      "Dispersal vs no-dispersal difference", "Model agreement map",
      "Turnover vs present maps", "Turnover vs present table",
      "Refugia score map", "Refugia score table", "Hotspot persistence map",
      "CWM trait maps", "Phylogenetic diversity maps", "Composition PCA maps",
      "Pseudo-fossil validation map", "Fossil validation table",
      "Pseudo-fossil probability histogram", "Pseudo-hindcast validation",
      "Simulation recovery richness", "Simulation recovery probability",
      "Simulation recovery table", "Coverage plot", "Bias plot",
      "Uncertainty budget", "Uncertainty budget table",
      "Full workflow Word report", "Full workflow R Markdown report",
      "Chinese explanation Word report", "Chinese explanation Word report ASCII copy",
      "Output file inventory", "Run metadata", "Case 1 run summary",
      "Case 2 run summary", "Ordered output index"
    ),
    stringsAsFactors = FALSE
  )
  exact <- rbind(
    exact,
    data.frame(
      basename = c(
        "28b_bgb_bsm_source_sink_turnover.png",
        "28b_bgb_bsm_source_sink_turnover.pdf",
        "bgb_bsm_region_metrics.csv",
        "bgb_bsm_dispersal_rates.csv",
        "bgb_bsm_source_sink_turnover_long.csv",
        "projection_model_definitions.csv",
        "projection_numeric_example.csv",
        "projection_numeric_example.png",
        "projection_numeric_example.pdf",
        "projection_component_summary_M5.csv",
        "projection_component_summary_M5.png",
        "projection_component_summary_M5.pdf",
        "M1_M5_projection_component_diagnostics.csv",
        "M3_M5_richness_through_time.csv",
        "M3_M5_richness_through_time.png",
        "M3_M5_occupied_area_through_time.csv",
        "M3_M5_occupied_area_through_time.png",
        "M3_M5_component_diagnostics_through_time.csv",
        "M3_M5_component_diagnostics_through_time.png",
        "M3_M5_representative_richness_maps.csv",
        "M3_M5_representative_richness_maps.png",
        "M3_M5_representative_probability_maps.csv",
        "M3_M5_representative_probability_maps.png",
        "M5_minus_M4_dynamic_vs_static_difference_representative.csv",
        "M5_minus_M4_dynamic_vs_static_difference.png",
        "M5_dynamic_dispersal_diagnostics.csv",
        "M5_dynamic_dispersal_diagnostics.png",
        "M5_dynamic_dispersal_representative_maps.csv",
        "M5_dynamic_dispersal_maps.png",
        "extrapolation_q_extrap_policy.csv"
      ),
      order_number = c(
        "3.1.12", "3.1.12", "3.1.13", "3.1.14", "3.1.15",
        "5.0.1", "5.0.2", "5.0.2", "5.0.2",
        "5.0.3", "5.0.3", "5.0.3",
        "5.0.4",
        "5.3.1", "5.3.1", "5.3.2", "5.3.2",
        "5.3.3", "5.3.3",
        "5.4.1", "5.4.1", "5.4.2", "5.4.2",
        "5.4.3", "5.4.3",
        "5.5.1", "5.5.1", "5.5.2", "5.5.2",
        "7.5.2"
      ),
      item_title = c(
        "BSM source sink turnover through time",
        "BSM source sink turnover through time",
        "BSM region source sink turnover metrics",
        "BSM dispersal rates",
        "BSM source sink turnover long table",
        "Projection model definitions",
        "Layered projection numeric example",
        "Layered projection numeric example",
        "Layered projection numeric example",
        "M5 projection component summary",
        "M5 projection component summary",
        "M5 projection component summary",
        "M1-M5 projection component diagnostics table",
        "M3-M5 richness through deep time",
        "M3-M5 richness through deep time",
        "M3-M5 occupied area through deep time",
        "M3-M5 occupied area through deep time",
        "M3-M5 component diagnostics through deep time",
        "M3-M5 component diagnostics through deep time",
        "M3-M5 representative richness maps",
        "M3-M5 representative richness maps",
        "M3-M5 representative probability maps",
        "M3-M5 representative probability maps",
        "M5 minus M4 dynamic vs static richness difference",
        "M5 minus M4 dynamic vs static richness difference",
        "M5 dynamic dispersal diagnostics through deep time",
        "M5 dynamic dispersal diagnostics through deep time",
        "M5 dynamic dispersal representative maps",
        "M5 dynamic dispersal representative maps",
        "Q_extrap extrapolation policy"
      ),
      stringsAsFactors = FALSE
    )
  )

  classify_deep_map <- function(base) {
    kind_map <- c(
      "deep_suitability_map_" = "4.1",
      "deep_accessibility_map_" = "4.2",
      "deep_phylo_mask_summary_" = "4.3",
      "deep_final_probability_example_map_" = "4.4",
      "deep_richness_map_" = "4.5",
      "deep_uncertainty_map_" = "4.6",
      "deep_extrapolation_risk_map_" = "4.7"
    )
    hit <- kind_map[startsWith(base, names(kind_map))]
    if (length(hit) == 0L) return(NULL)
    tm <- sub(".*_([0-9.]+)Ma\\.[^.]+$", "\\1", base)
    tm_num <- suppressWarnings(as.numeric(tm))
    if (is.na(tm_num)) tm_num <- -Inf
    list(base_number = unname(hit[1]), time_ma = tm_num)
  }

  out <- vector("list", length(rel))
  for (i in seq_along(rel)) {
    base <- basename(rel[i])
    folder <- dirname(rel[i])
    ext <- tools::file_ext(base)
    file_type <- if (nzchar(ext)) tolower(ext) else "file"
    exact_i <- match(base, exact$basename)
    if (is.na(exact_i) && identical(file_type, "pdf")) {
      exact_i <- match(paste0(tools::file_path_sans_ext(base), ".png"), exact$basename)
    }
    if (!is.na(exact_i)) {
      num <- exact$order_number[exact_i]
      title <- exact$item_title[exact_i]
      dynamic_base <- NA_character_
      time_ma <- NA_real_
    } else {
      dm <- classify_deep_map(base)
      if (!is.null(dm)) {
        num <- NA_character_
        title <- sub("\\.[^.]+$", "", base)
        dynamic_base <- dm$base_number
        time_ma <- dm$time_ma
      } else if (grepl("^rds/", rel[i])) {
        num <- NA_character_
        title <- sub("\\.[^.]+$", "", base)
        dynamic_base <- "9.1"
        time_ma <- NA_real_
      } else if (grepl("^metadata/", rel[i])) {
        num <- NA_character_
        title <- sub("\\.[^.]+$", "", base)
        dynamic_base <- "8.3"
        time_ma <- NA_real_
      } else if (grepl("^tasks/", rel[i])) {
        num <- NA_character_
        title <- sub("\\.[^.]+$", "", base)
        dynamic_base <- "4.8"
        time_ma <- NA_real_
      } else if (grepl("^reports/", rel[i])) {
        num <- NA_character_
        title <- sub("\\.[^.]+$", "", base)
        dynamic_base <- "8.1"
        time_ma <- NA_real_
      } else {
        num <- NA_character_
        title <- sub("\\.[^.]+$", "", base)
        dynamic_base <- "8.9"
        time_ma <- NA_real_
      }
    }
    section <- if (!is.na(num)) strsplit(num, "\\.")[[1]][1] else strsplit(dynamic_base, "\\.")[[1]][1]
    out[[i]] <- data.frame(
      order_number = num,
      dynamic_base = dynamic_base,
      time_ma = time_ma,
      section = section,
      section_title = section_titles[section],
      item_title = title,
      file_type = file_type,
      original_relpath = rel[i],
      stringsAsFactors = FALSE
    )
  }
  catalog <- do.call(rbind, out)

  fill_dynamic <- is.na(catalog$order_number)
  if (any(fill_dynamic)) {
    bases <- unique(catalog$dynamic_base[fill_dynamic])
    for (b in bases) {
      idx <- which(fill_dynamic & catalog$dynamic_base == b)
      if (startsWith(b, "4.") && any(is.finite(catalog$time_ma[idx]))) {
        stems <- tools::file_path_sans_ext(basename(catalog$original_relpath[idx]))
        unique_stems <- unique(data.frame(
          stem = stems,
          time_ma = catalog$time_ma[idx],
          stringsAsFactors = FALSE
        ))
        unique_stems <- unique_stems[order(-unique_stems$time_ma, unique_stems$stem), ]
        for (j in seq_len(nrow(unique_stems))) {
          same <- idx[stems == unique_stems$stem[j]]
          catalog$order_number[same] <- paste0(b, ".", sprintf("%02d", j))
        }
        next
      } else {
        idx <- idx[order(catalog$original_relpath[idx])]
      }
      catalog$order_number[idx] <- paste0(b, ".", sprintf("%02d", seq_along(idx)))
    }
  }

  catalog$section_dir <- section_dirs[catalog$section]
  catalog$ordered_filename <- paste0(
    .hee_filename_order_number(catalog$order_number), "_",
    .hee_sanitize_ordered_name(catalog$item_title),
    ifelse(nzchar(catalog$file_type), paste0(".", catalog$file_type), "")
  )
  catalog$ordered_relpath <- file.path(catalog$section_dir, catalog$ordered_filename)
  catalog <- catalog[order(.hee_order_key(catalog$order_number), catalog$original_relpath), ]
  rownames(catalog) <- NULL
  catalog[, c("order_number", "section", "section_title", "item_title", "file_type",
              "original_relpath", "ordered_relpath")]
}

#' Write ordered output copies
#'
#' Copy or hard-link workflow outputs into an `ordered_results` folder with
#' hierarchical numeric prefixes. This leaves the original output tree intact,
#' so existing report links are not broken.
#'
#' @param output_root Case output directory.
#' @param ordered_dir Destination directory. Defaults to
#'   `file.path(output_root, "ordered_results")`.
#' @param include_rds Logical. Include RDS objects.
#' @param overwrite Logical. Remove an existing `ordered_dir` before writing.
#' @param copy_mode Either `"link_or_copy"` or `"copy"`. Hard links avoid
#'   duplicating large RDS files when supported by the filesystem.
#'
#' @return The ordered catalog data frame, invisibly.
#' @export
#'
#' @examples
#' \dontrun{
#' hee_write_ordered_outputs("outputs/case01_simulated_full_540Ma")
#' }
hee_write_ordered_outputs <- function(output_root,
                                      ordered_dir = file.path(output_root, "ordered_results"),
                                      include_rds = TRUE,
                                      overwrite = TRUE,
                                      copy_mode = c("link_or_copy", "copy")) {
  copy_mode <- match.arg(copy_mode)
  output_root <- normalizePath(output_root, winslash = "/", mustWork = TRUE)
  ordered_dir <- normalizePath(ordered_dir, winslash = "/", mustWork = FALSE)
  catalog <- hee_ordered_output_catalog(output_root, include_rds = include_rds)
  if (overwrite && dir.exists(ordered_dir)) unlink(ordered_dir, recursive = TRUE, force = TRUE)
  dir.create(ordered_dir, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_len(nrow(catalog))) {
    src <- file.path(output_root, catalog$original_relpath[i])
    dst <- file.path(ordered_dir, catalog$ordered_relpath[i])
    dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(dst)) unlink(dst, force = TRUE)
    linked <- FALSE
    if (copy_mode == "link_or_copy") {
      linked <- tryCatch(file.link(src, dst), error = function(e) FALSE)
    }
    if (!linked) {
      ok <- file.copy(src, dst, overwrite = TRUE, copy.date = TRUE)
      if (!ok) stop("Failed to copy ordered output: ", src, " -> ", dst, call. = FALSE)
    }
  }

  index_path <- file.path(ordered_dir, "00_ordered_output_index.csv")
  utils::write.csv(catalog, index_path, row.names = FALSE, fileEncoding = "UTF-8")
  tab_dir <- file.path(output_root, "tables")
  if (dir.exists(tab_dir)) {
    utils::write.csv(catalog, file.path(tab_dir, "ordered_output_index.csv"),
                     row.names = FALSE, fileEncoding = "UTF-8")
  }
  invisible(catalog)
}

.hee_sanitize_ordered_name <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT", sub = "")
  x <- gsub("^[0-9]+[A-Za-z]*_+", "", x)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x[!nzchar(x)] <- "output"
  x
}

.hee_order_key <- function(number) {
  vapply(strsplit(number, "\\."), function(parts) {
    nums <- suppressWarnings(as.numeric(parts))
    nums[is.na(nums)] <- 0
    sum(nums / (1000 ^ seq_along(nums)))
  }, numeric(1))
}

.hee_filename_order_number <- function(number) {
  vapply(strsplit(number, "\\."), function(parts) {
    if (length(parts) >= 3L) {
      last <- suppressWarnings(as.integer(parts[length(parts)]))
      if (!is.na(last)) parts[length(parts)] <- sprintf("%02d", last)
    }
    paste(parts, collapse = ".")
  }, character(1))
}
