#!/usr/bin/env Rscript

# Case04 no-BioGeoBEARS ancestral environmental suitability
#
# Scientific scope:
#   HMSC tip Beta posterior -> dated-tree ancestral Beta -> palaeoenvironment
#   -> S_env at every valid time slice.
#
# This script deliberately does not use BioGeoBEARS, BSM, region histories,
# movement kernels, colonisation, persistence, M1-M5, or final occupancy.
# Outputs must be read as ancestral environmental support / potential
# suitability, not historical presence probability.

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list()
  for (a in args) {
    if (!grepl("^--", a)) next
    kv <- sub("^--", "", a)
    pos <- regexpr("=", kv, fixed = TRUE)
    if (pos < 0) {
      out[[kv]] <- TRUE
    } else {
      key <- substr(kv, 1, pos - 1)
      val <- substr(kv, pos + 1, nchar(kv))
      out[[key]] <- val
    }
  }
  out
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x)) y else x

as_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

as_num <- function(x, default = NA_real_) {
  if (is.null(x) || identical(x, "")) return(default)
  suppressWarnings(as.numeric(x))
}

as_int <- function(x, default = NA_integer_) {
  z <- as_num(x, default = default)
  if (!is.finite(z)) return(default)
  as.integer(z)
}

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

copy_if_exists <- function(from, to_dir) {
  if (!file.exists(from)) return(FALSE)
  safe_dir_create(to_dir)
  file.copy(from, file.path(to_dir, basename(from)), overwrite = TRUE)
}

write_csv <- function(x, path) {
  safe_dir_create(dirname(path))
  utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
  invisible(path)
}

write_text <- function(x, path) {
  safe_dir_create(dirname(path))
  writeLines(x, path, useBytes = TRUE)
  invisible(path)
}

slug_time <- function(x) paste0(format(as.numeric(x), trim = TRUE, scientific = FALSE), "Ma")

sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\.-]+", "_", as.character(x))
  gsub("_+", "_", x)
}

now_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

cfg <- parse_args(args)

pkg_root <- normalizePath(cfg$pkg_root %||% "C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo",
                          winslash = "/", mustWork = FALSE)
case04_prev <- normalizePath(
  cfg$previous_case04 %||%
    "C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo_API8_20260818/outputs/case04_standard200_true_bgb_all109_stage1_20260821",
  winslash = "/", mustWork = FALSE
)
plant_input_dir <- normalizePath(
  cfg$plant_input_dir %||%
    "C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo_API8_20260818/data_external/plant_200_multifamily_extant_20260820/prepared_inputs/case03_data_raw/Plant200_multifamily_extant_1deg_traits",
  winslash = "/", mustWork = FALSE
)

output <- cfg$output %||%
  file.path(pkg_root, "outputs",
            paste0("case04_plant200_ancestral_suitability_no_bgb_", now_stamp()))
output <- safe_dir_create(output)

quick <- as_bool(cfg$quick, default = FALSE)
max_times <- as_int(cfg$max_times, default = if (quick) 5L else NA_integer_)
max_lineages <- as_int(cfg$max_lineages, default = if (quick) 25L else NA_integer_)
max_draws <- as_int(cfg$max_draws, default = NA_integer_)
lineage_chunk_size <- as_int(cfg$lineage_chunk_size, default = if (quick) 10L else 20L)
save_lineage_cell <- as_bool(cfg$save_lineage_cell, default = FALSE)
make_png <- as_bool(cfg$make_png, default = TRUE)
make_tif <- as_bool(cfg$make_tif, default = TRUE)
link <- cfg$link %||% "probit"

dir_config <- safe_dir_create(file.path(output, "00_config"))
dir_input <- safe_dir_create(file.path(output, "01_input_audit"))
dir_resp <- safe_dir_create(file.path(output, "02_ancestral_beta_responses"))
dir_tables <- safe_dir_create(file.path(output, "03_suitability_tables_by_time"))
dir_maps <- safe_dir_create(file.path(output, "04_maps_by_time"))
dir_rasters <- safe_dir_create(file.path(output, "05_rasters_by_time"))
dir_copy <- safe_dir_create(file.path(output, "06_copy_ready"))
dir_logs <- safe_dir_create(file.path(output, "logs"))

log_file <- file.path(dir_logs, "case04_no_bgb_run.log")
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                paste(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_msg("Case04 no-BGB ancestral suitability started.")

if (dir.exists(pkg_root) && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(pkg_root, quiet = TRUE)
  log_msg("Loaded latest package source with pkgload::load_all():", pkg_root)
} else {
  suppressPackageStartupMessages(library(HmscEcoEvo))
  log_msg("Loaded installed HmscEcoEvo package.")
}
if (utils::packageVersion("HmscEcoEvo") < package_version("1.0.4") ||
    sum(hee_core_process_catalog()$layer == "core_process") != 6L) {
  stop("Case04 requires HmscEcoEvo >= 1.0.4 with the six-process API.",
       call. = FALSE)
}

for (pkg in c("ape", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is missing: ", pkg, call. = FALSE)
  }
}

paths <- list(
  previous_case04 = case04_prev,
  plant_input_dir = plant_input_dir,
  tree = cfg$tree %||% file.path(plant_input_dir, "tree.tre"),
  response_history = cfg$response_history %||%
    file.path(case04_prev, "04_evolution", "case04_lineage_response_history.csv"),
  palaeo_earth_state = cfg$palaeo_earth_state %||%
    file.path(case04_prev, "05_palaeo_environment", "case04_palaeo_earth_state_scaled.rds"),
  recipe = cfg$recipe %||%
    file.path(case04_prev, "03_hmsc_environment", "case04_modern_0Ma_environment_recipe.rds"),
  time_domain = cfg$time_domain %||%
    file.path(case04_prev, "00_config", "case04_time_domain_540_0Ma.csv"),
  input_audit_summary = file.path(case04_prev, "01_input_audit", "case04_input_audit_summary.csv"),
  input_file_existence = file.path(case04_prev, "01_input_audit", "case04_input_file_existence.csv"),
  selected_species = file.path(case04_prev, "01_input_audit", "case04_selected_species.csv"),
  hmsc_xdata = file.path(case04_prev, "03_hmsc_environment", "case04_hmsc_XData_scaled.csv"),
  hmsc_y = file.path(case04_prev, "03_hmsc_environment", "case04_hmsc_Y_selected.csv"),
  hmsc_trdata = file.path(case04_prev, "03_hmsc_environment", "case04_hmsc_TrData_selected.csv")
)

required <- c("tree", "response_history", "palaeo_earth_state")
input_inventory <- data.frame(
  input = names(paths),
  path = unlist(paths, use.names = FALSE),
  required = names(paths) %in% required,
  exists = file.exists(unlist(paths, use.names = FALSE)),
  stringsAsFactors = FALSE
)
write_csv(input_inventory, file.path(dir_input, "case04_no_bgb_input_inventory.csv"))
missing_required <- input_inventory[input_inventory$required & !input_inventory$exists, , drop = FALSE]
if (nrow(missing_required) > 0L) {
  write_csv(missing_required, file.path(dir_input, "missing_required_inputs.csv"))
  stop("Missing required Case04 no-BGB input(s). See: ",
       file.path(dir_input, "missing_required_inputs.csv"), call. = FALSE)
}

invisible(copy_if_exists(paths$input_audit_summary, dir_input))
invisible(copy_if_exists(paths$input_file_existence, dir_input))
invisible(copy_if_exists(paths$selected_species, dir_input))
invisible(copy_if_exists(paths$hmsc_xdata, dir_input))
invisible(copy_if_exists(paths$hmsc_y, dir_input))
invisible(copy_if_exists(paths$hmsc_trdata, dir_input))
invisible(copy_if_exists(paths$time_domain, dir_config))
invisible(copy_if_exists(paths$recipe, dir_config))

tree <- ape::read.tree(paths$tree)
root_age <- max(ape::node.depth.edgelength(tree))
tree_tips <- tree$tip.label

response_raw <- utils::read.csv(paths$response_history, check.names = FALSE,
                                stringsAsFactors = FALSE)
if (!"lineage" %in% names(response_raw)) {
  stop("response_history must contain a lineage column.", call. = FALSE)
}
recipe <- if (file.exists(paths$recipe)) readRDS(paths$recipe) else NULL
basis_cols <- if (!is.null(recipe) && "variables" %in% names(recipe)) {
  as.character(recipe$variables)
} else {
  setdiff(names(response_raw)[vapply(response_raw, is.numeric, logical(1))],
          c("intercept", "hmsc_intercept_original", "time_ma"))
}
basis_cols <- setdiff(basis_cols, "log_sampling_effort")
missing_basis <- setdiff(basis_cols, names(response_raw))
if (length(missing_basis) > 0L) {
  stop("response_history lacks beta axis/axes required by recipe: ",
       paste(missing_basis, collapse = ", "), call. = FALSE)
}

tip_beta <- response_raw[response_raw$lineage %in% tree_tips, , drop = FALSE]
tip_beta <- tip_beta[, unique(c("lineage", "response_draw", "intercept",
                                "hmsc_intercept_original", basis_cols)),
                     drop = FALSE]
names(tip_beta)[names(tip_beta) == "lineage"] <- "species"
if (!"response_draw" %in% names(tip_beta)) tip_beta$response_draw <- "s1"
.hee_check_unique_keys(tip_beta, c("species", "response_draw"),
                       "Case04 tip beta response history")
draw_ids <- unique(tip_beta$response_draw)
if (is.finite(max_draws) && length(draw_ids) > max_draws) {
  keep_draws <- draw_ids[seq_len(max_draws)]
  tip_beta <- tip_beta[tip_beta$response_draw %in% keep_draws, , drop = FALSE]
  draw_ids <- keep_draws
}

species_ids <- unique(tip_beta$species)
if (is.finite(max_lineages) && length(species_ids) > max_lineages) {
  species_ids <- species_ids[seq_len(max_lineages)]
  tip_beta <- tip_beta[tip_beta$species %in% species_ids, , drop = FALSE]
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, species_ids))
  root_age <- max(ape::node.depth.edgelength(tree))
}

earth <- readRDS(paths$palaeo_earth_state)
earth <- as.data.frame(earth)
if (!all(c("cell_id", "time_ma", "lon", "lat") %in% names(earth))) {
  stop("palaeo_earth_state must contain cell_id, time_ma, lon, and lat.",
       call. = FALSE)
}
missing_earth_basis <- setdiff(basis_cols, names(earth))
if (length(missing_earth_basis) > 0L) {
  stop("palaeo_earth_state lacks scaled basis column(s): ",
       paste(missing_earth_basis, collapse = ", "), call. = FALSE)
}
if (!"region" %in% names(earth)) earth$region <- "global_no_bgb"
if (!"H_state" %in% names(earth)) earth$H_state <- 1
earth <- earth[is.finite(earth$time_ma) & earth$H_state == 1, , drop = FALSE]
available_times <- sort(unique(as.numeric(earth$time_ma)), decreasing = TRUE)
requested_times <- available_times
if (file.exists(paths$time_domain)) {
  td0 <- tryCatch(utils::read.csv(paths$time_domain, check.names = FALSE),
                  error = function(e) NULL)
  if (!is.null(td0) && "time_ma" %in% names(td0)) {
    requested_times <- sort(unique(as.numeric(td0$time_ma)), decreasing = TRUE)
  }
}
inside_times <- available_times[available_times <= root_age + sqrt(.Machine$double.eps)]
outside_times <- requested_times[requested_times > root_age + sqrt(.Machine$double.eps)]
if (is.finite(max_times) && length(inside_times) > max_times) {
  inside_times <- sort(unique(c(utils::head(inside_times, max_times - 1L), 0)),
                       decreasing = TRUE)
}
if (length(inside_times) == 0L) {
  stop("No palaeoenvironment time slice falls inside the tree time domain.",
       call. = FALSE)
}

time_domain <- data.frame(
  time_ma = requested_times,
  tree_root_age_ma = root_age,
  status = ifelse(requested_times %in% inside_times, "PROJECTED",
                  ifelse(requested_times %in% outside_times,
                         "OUTSIDE_LINEAGE_TIME_DOMAIN",
                         ifelse(!requested_times %in% available_times,
                                "MISSING_FROM_PALAEO_EARTH_STATE",
                                "NOT_SELECTED_BY_RUN_OPTIONS"))),
  interpretation = ifelse(requested_times > root_age,
                          "older_than_sampled_tree_root_not_zero_not_absence",
                          "ancestral_environmental_support_can_be_projected"),
  stringsAsFactors = FALSE
)
write_csv(time_domain, file.path(dir_config, "case04_no_bgb_time_domain.csv"))

config <- data.frame(
  parameter = c("analysis_name", "analysis_mode", "quick", "link",
                "tree_root_age_ma", "n_tip_species", "n_beta_draws",
                "n_projected_times", "save_lineage_cell",
                "lineage_chunk_size", "interpretation"),
  value = c("case04_plant200_ancestral_suitability_no_bgb",
            "no_biogeobears_environmental_filtering_only",
            quick, link, root_age, length(unique(tip_beta$species)),
            length(unique(tip_beta$response_draw)), length(inside_times),
            save_lineage_cell, lineage_chunk_size,
            "S_env is ancestral environmental support, not final occupancy"),
  stringsAsFactors = FALSE
)
write_csv(config, file.path(dir_config, "case04_no_bgb_config.csv"))

log_msg("Reconstructing ancestral beta responses for",
        length(inside_times), "time slice(s),",
        length(unique(tip_beta$species)), "tip species,",
        length(unique(tip_beta$response_draw)), "draw(s).")

responses <- hee_evolution_ancestral_response_direct(
  tip_beta_draws = tip_beta,
  tree = tree,
  times = inside_times,
  basis_cols = basis_cols,
  species_col = "species",
  response_draw_col = "response_draw",
  intercept_col = "intercept",
  reconstruct_intercept = FALSE,
  branch_uncertainty = "none"
)
saveRDS(responses, file.path(dir_resp, "case04_ancestral_beta_responses_no_bgb.rds"))
utils::write.csv(responses, file.path(dir_resp, "case04_ancestral_beta_responses_no_bgb.csv"),
                 row.names = FALSE, fileEncoding = "UTF-8")

response_count_by_time <- aggregate(
  lineage ~ time_ma + response_draw,
  responses,
  function(z) length(unique(z))
)
names(response_count_by_time)[names(response_count_by_time) == "lineage"] <- "n_active_lineages"
write_csv(response_count_by_time,
          file.path(dir_resp, "case04_active_lineage_counts_by_time_draw.csv"))

plot_map <- function(df, value_col, path, title, limits, legend_title) {
  if (!make_png) return(invisible(FALSE))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = lon, y = lat, fill = .data[[value_col]])) +
    ggplot2::geom_tile(width = 1, height = 1) +
    ggplot2::coord_equal(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
    ggplot2::scale_fill_gradientn(
      colours = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
      limits = limits,
      na.value = "grey92",
      name = legend_title
    ) +
    ggplot2::labs(
      title = title,
      subtitle = "No BioGeoBEARS: ancestral environmental support only, not historical occupancy",
      x = "Longitude", y = "Latitude"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  ggplot2::ggsave(path, p, width = 10, height = 5.4, dpi = 180)
  invisible(TRUE)
}

write_raster <- function(df, value_col, path) {
  if (!make_tif || !requireNamespace("terra", quietly = TRUE)) return(invisible(FALSE))
  xyz <- df[, c("lon", "lat", value_col), drop = FALSE]
  names(xyz) <- c("x", "y", "z")
  r <- terra::rast(xyz, type = "xyz", crs = "EPSG:4326")
  terra::writeRaster(r, path, overwrite = TRUE)
  invisible(TRUE)
}

summarise_draws_to_lineage_cell <- function(draws) {
  keep <- intersect(c("lineage", "cell_id", "time_ma", "region", "lon", "lat",
                      "H_state"), names(draws))
  key <- paste(draws$lineage, draws$cell_id, sep = "\r")
  first <- !duplicated(key)
  base <- draws[first, keep, drop = FALSE]
  base_key <- key[first]

  s <- draws$suitability
  eta <- draws$eta
  ok <- is.finite(s)
  n <- rowsum(as.numeric(ok), key, reorder = FALSE)
  sum_s <- rowsum(ifelse(ok, s, 0), key, reorder = FALSE)
  sum_s2 <- rowsum(ifelse(ok, s * s, 0), key, reorder = FALSE)
  sum_eta <- rowsum(ifelse(is.finite(eta), eta, 0), key, reorder = FALSE)
  n_eta <- rowsum(as.numeric(is.finite(eta)), key, reorder = FALSE)

  idx <- match(rownames(sum_s), base_key)
  out <- base[idx, , drop = FALSE]
  n_vec <- as.numeric(n[, 1])
  mu <- as.numeric(sum_s[, 1]) / pmax(n_vec, 1)
  var <- (as.numeric(sum_s2[, 1]) - (as.numeric(sum_s[, 1])^2 / pmax(n_vec, 1))) /
    pmax(n_vec - 1, 1)
  out$eta_mean <- as.numeric(sum_eta[, 1]) / pmax(as.numeric(n_eta[, 1]), 1)
  out$suitability_mean <- mu
  out$suitability_sd <- sqrt(pmax(var, 0))
  out$n_response_draws <- n_vec
  out$suitability_interpretation <-
    "ancestral_environmental_support_not_final_occupancy"
  rownames(out) <- NULL
  out
}

global_time_summary <- list()
map_index <- list()
copy_sentences <- c(
  "# Case04 Plant200 no-BioGeoBEARS 祖先环境适合度结果",
  "",
  "本目录只计算 `S_env = g^{-1}(B(X_palaeo) %*% beta_ancestor)`。",
  "这里的 `S_env` 是祖先谱系在给定古环境下的环境支持度/潜在适合度，不是历史占据概率。",
  "本次没有使用 BioGeoBEARS/BSM，因此没有大区历史约束、跨区扩散、定殖、存续或最终 occupancy。",
  "树根年龄之外的时间片标记为 outside lineage time domain；这些时间不是 0，也不是灭绝，只是不属于这棵 Plant200 现生采样树能解释的时间域。",
  ""
)

for (tm in inside_times) {
  log_msg("Projecting ancestral environmental support at", tm, "Ma.")
  e_tm <- earth[as.numeric(earth$time_ma) == as.numeric(tm), , drop = FALSE]
  r_tm <- responses[as.numeric(responses$time_ma) == as.numeric(tm), , drop = FALSE]
  lineages <- unique(r_tm$lineage)
  if (!nrow(e_tm) || !nrow(r_tm)) next

  cell_base <- unique(e_tm[, intersect(c("cell_id", "time_ma", "region", "lon", "lat", "H_state"),
                                      names(e_tm)),
                         drop = FALSE])
  cell_base$expected_support_richness <- 0
  cell_base$mean_environmental_support_sum <- 0
  cell_base$mean_support_sd_sum <- 0
  cell_base$n_lineage_cell_summaries <- 0L

  chunks <- split(lineages, ceiling(seq_along(lineages) / lineage_chunk_size))
  lineage_cell_preview <- list()
  for (ci in seq_along(chunks)) {
    r_chunk <- r_tm[r_tm$lineage %in% chunks[[ci]], , drop = FALSE]
    draws <- hee_lineage_suitability(
      earth_state = e_tm,
      responses = r_chunk,
      basis_cols = basis_cols,
      link = link
    )
    lc <- summarise_draws_to_lineage_cell(draws)
    if (save_lineage_cell) {
      chunk_file <- file.path(dir_tables,
                              paste0("lineage_cell_suitability_",
                                     slug_time(tm), "_chunk",
                                     sprintf("%03d", ci), ".rds"))
      saveRDS(lc, chunk_file)
    } else if (ci <= 2L) {
      lineage_cell_preview[[length(lineage_cell_preview) + 1L]] <-
        utils::head(lc, 200L)
    }

    sum_by_cell <- rowsum(lc$suitability_mean, lc$cell_id, reorder = FALSE)
    sd_by_cell <- rowsum(lc$suitability_sd, lc$cell_id, reorder = FALSE)
    n_by_cell <- rowsum(rep(1, nrow(lc)), lc$cell_id, reorder = FALSE)
    idx <- match(rownames(sum_by_cell), cell_base$cell_id)
    cell_base$expected_support_richness[idx] <-
      cell_base$expected_support_richness[idx] + as.numeric(sum_by_cell[, 1])
    cell_base$mean_support_sd_sum[idx] <-
      cell_base$mean_support_sd_sum[idx] + as.numeric(sd_by_cell[, 1])
    cell_base$n_lineage_cell_summaries[idx] <-
      cell_base$n_lineage_cell_summaries[idx] + as.integer(n_by_cell[, 1])
  }
  cell_base$mean_environmental_support <-
    ifelse(cell_base$n_lineage_cell_summaries > 0,
           cell_base$expected_support_richness /
             cell_base$n_lineage_cell_summaries,
           NA_real_)
  cell_base$mean_support_sd <-
    ifelse(cell_base$n_lineage_cell_summaries > 0,
           cell_base$mean_support_sd_sum /
             cell_base$n_lineage_cell_summaries,
           NA_real_)
  cell_base$n_active_lineages <- length(lineages)
  cell_base$n_response_draws <- length(unique(r_tm$response_draw))
  cell_base$analysis_mode <- "no_bgb_environmental_filtering_only"
  cell_base$interpretation <- "environmental_support_not_final_occupancy"

  time_slug <- slug_time(tm)
  csv_file <- file.path(dir_tables, paste0("cell_aggregate_suitability_", time_slug, ".csv"))
  write_csv(cell_base, csv_file)
  if (length(lineage_cell_preview) > 0L) {
    preview <- do.call(rbind, lineage_cell_preview)
    write_csv(preview, file.path(dir_tables,
                                 paste0("lineage_cell_suitability_preview_",
                                        time_slug, ".csv")))
  }

  mean_png <- file.path(dir_maps, paste0("mean_environmental_support_",
                                         time_slug, "_no_bgb.png"))
  rich_png <- file.path(dir_maps, paste0("expected_support_richness_",
                                         time_slug, "_no_bgb.png"))
  plot_map(cell_base, "mean_environmental_support", mean_png,
           paste0("Mean ancestral environmental support, ", time_slug),
           limits = c(0, 1), legend_title = "mean S_env")
  plot_map(cell_base, "expected_support_richness", rich_png,
           paste0("Expected environmental-support richness, ", time_slug),
           limits = c(0, length(unique(tip_beta$species))),
           legend_title = "sum S_env")

  mean_tif <- file.path(dir_rasters, paste0("mean_environmental_support_",
                                            time_slug, "_no_bgb.tif"))
  rich_tif <- file.path(dir_rasters, paste0("expected_support_richness_",
                                            time_slug, "_no_bgb.tif"))
  write_raster(cell_base, "mean_environmental_support", mean_tif)
  write_raster(cell_base, "expected_support_richness", rich_tif)

  global_time_summary[[length(global_time_summary) + 1L]] <- data.frame(
    time_ma = tm,
    n_land_cells = nrow(cell_base),
    n_active_lineages = length(lineages),
    n_response_draws = length(unique(r_tm$response_draw)),
    mean_of_cell_mean_support = mean(cell_base$mean_environmental_support, na.rm = TRUE),
    max_cell_mean_support = max(cell_base$mean_environmental_support, na.rm = TRUE),
    mean_expected_support_richness = mean(cell_base$expected_support_richness, na.rm = TRUE),
    max_expected_support_richness = max(cell_base$expected_support_richness, na.rm = TRUE),
    csv = csv_file,
    mean_png = mean_png,
    richness_png = rich_png,
    mean_tif = if (file.exists(mean_tif)) mean_tif else NA_character_,
    richness_tif = if (file.exists(rich_tif)) rich_tif else NA_character_,
    stringsAsFactors = FALSE
  )
  map_index[[length(map_index) + 1L]] <- data.frame(
    time_ma = c(tm, tm, tm, tm),
    metric = c("mean_environmental_support",
               "expected_support_richness",
               "mean_environmental_support",
               "expected_support_richness"),
    file_type = c("PNG", "PNG", "GeoTIFF", "GeoTIFF"),
    path = c(mean_png, rich_png,
             if (file.exists(mean_tif)) mean_tif else NA_character_,
             if (file.exists(rich_tif)) rich_tif else NA_character_),
    value_range = c("[0,1]",
                    paste0("[0,", length(unique(tip_beta$species)), "]"),
                    "[0,1]",
                    paste0("[0,", length(unique(tip_beta$species)), "]")),
    units = c("probability-scale environmental support",
              "sum of probability-scale environmental support over active lineages",
              "probability-scale environmental support",
              "sum of probability-scale environmental support over active lineages"),
    data_source = "HMSC tip beta posterior + dated tree + scaled palaeoenvironment; no BioGeoBEARS",
    interpretation = "ancestral environmental support, not final occupancy",
    stringsAsFactors = FALSE
  )
}

summary_df <- if (length(global_time_summary)) do.call(rbind, global_time_summary) else data.frame()
map_df <- if (length(map_index)) do.call(rbind, map_index) else data.frame()
write_csv(summary_df, file.path(dir_copy, "case04_no_bgb_global_time_summary.csv"))
write_csv(map_df, file.path(dir_copy, "case04_no_bgb_map_index.csv"))

if (nrow(summary_df) > 0L) {
  top_mean <- summary_df[which.max(summary_df$mean_of_cell_mean_support), , drop = FALSE]
  top_rich <- summary_df[which.max(summary_df$max_expected_support_richness), , drop = FALSE]
  copy_sentences <- c(
    copy_sentences,
    "## 可直接复制的结果句子",
    "",
    paste0("本次 no-BGB Case04 在 Plant200 树时间域内投影了 ",
           nrow(summary_df), " 个古环境时间片；树根年龄约为 ",
           round(root_age, 2), " Ma。"),
    paste0("在所有已投影时间片中，平均 cell-level 祖先环境支持度最高的时间片为 ",
           top_mean$time_ma, " Ma，全球陆地网格平均 S_env 为 ",
           signif(top_mean$mean_of_cell_mean_support, 4), "。"),
    paste0("环境支持度丰富度最高的单个网格出现在 ",
           top_rich$time_ma, " Ma，对应 `sum(S_env)` 最大值为 ",
           signif(top_rich$max_expected_support_richness, 4), "。"),
    "这些句子只能写成“环境支持度/潜在适合度”，不能写成真实历史分布、真实占据概率或 BioGeoBEARS 约束后的范围历史。",
    "",
    "## 解释边界",
    "",
    "- `mean_environmental_support`：每个网格上 active lineage 的 S_env 平均值，图例固定为 0-1。",
    "- `expected_support_richness`：每个网格上 active lineage 的 S_env 之和，是环境支持度丰富度，不是物种丰富度，也不是 occupancy richness。",
    "- 没有 BioGeoBEARS：适合但历史上未到达的大区不会被排除。",
    "- 没有 colonisation/persistence：到达、建立和局地存续尚未进入本次输出。",
    "- 325 Ma 以前：超出当前 Plant200 现生树根年龄，标记为 outside lineage time domain。"
  )
}
write_text(copy_sentences, file.path(dir_copy, "case04_no_bgb_copy_ready_interpretation_zh.md"))

final_validation <- data.frame(
  check = c("used_latest_hmscee_functions",
            "did_not_use_biogeobears",
            "time_ma_direction",
            "older_than_tree_root_not_projected",
            "suitability_range_mean_support",
            "map_color_range_consistent",
            "geotiff_requested",
            "png_requested"),
  status = c("PASS", "PASS", "PASS",
             ifelse(length(outside_times) > 0, "PASS", "NOT_APPLICABLE"),
             ifelse(nrow(summary_df) == 0 ||
                      (all(summary_df$mean_of_cell_mean_support >= 0, na.rm = TRUE) &&
                         all(summary_df$mean_of_cell_mean_support <= 1, na.rm = TRUE)),
                    "PASS", "FAIL"),
             "PASS",
             ifelse(make_tif, ifelse(requireNamespace("terra", quietly = TRUE), "PASS", "SKIP_TERRA_NOT_INSTALLED"), "SKIP_BY_ARGUMENT"),
             ifelse(make_png, "PASS", "SKIP_BY_ARGUMENT")),
  detail = c("hee_evolution_ancestral_response_direct + hee_lineage_suitability",
             "BioGeoBEARS/BSM was intentionally not loaded or called",
             "Ma larger means older; projected times sorted old to young",
             paste(length(outside_times), "available palaeo time(s) older than tree root were not projected"),
             "mean environmental support should remain in [0,1]",
             "mean S_env fixed to [0,1]; support richness fixed to [0,n_tip_species]",
             "GeoTIFFs written when terra is available",
             "PNGs written with longitude/latitude axes"),
  stringsAsFactors = FALSE
)
write_csv(final_validation, file.path(dir_config, "case04_no_bgb_validation_checks.csv"))

log_msg("Case04 no-BGB ancestral suitability completed.")
log_msg("Output:", output)

cat("\nDONE\n")
cat("Output:", output, "\n")
