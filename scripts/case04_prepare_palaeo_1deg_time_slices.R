#!/usr/bin/env Rscript

# Prepare compact, complete 1-degree palaeo-grid time slices for Case04/05.
# The source RDS is read once here; dynamic workers subsequently load one full
# time slice at a time rather than the complete 540--0 Ma cube.

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x)) y else x
parse_args <- function(x) {
  out <- list()
  for (arg in x) {
    if (!grepl("^--", arg)) next
    kv <- sub("^--", "", arg)
    p <- regexpr("=", kv, fixed = TRUE)
    if (p < 0L) out[[kv]] <- TRUE else {
      out[[substr(kv, 1L, p - 1L)]] <- substr(kv, p + 1L, nchar(kv))
    }
  }
  out
}
safe_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}
write_csv <- function(x, path) {
  safe_dir(dirname(path))
  if (requireNamespace("data.table", quietly = TRUE)) data.table::fwrite(x, path) else {
    utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
  }
}

cfg <- parse_args(args)
pkg_root <- normalizePath(cfg$pkg_root %||%
                            "C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo",
                          winslash = "/", mustWork = FALSE)
previous_case04 <- normalizePath(
  cfg$previous_case04 %||%
    "C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo_API8_20260818/outputs/case04_standard200_true_bgb_all109_stage1_20260821",
  winslash = "/", mustWork = FALSE
)
output <- safe_dir(cfg$output %||%
                     file.path(pkg_root, "derived_inputs",
                               "case04_plant200_palaeo_1deg_slices_20260919"))
earth_path <- normalizePath(
  cfg$palaeo_earth_state %||%
    file.path(previous_case04, "05_palaeo_environment", "case04_palaeo_earth_state_scaled.rds"),
  winslash = "/", mustWork = FALSE
)
topography_path <- normalizePath(
  cfg$topographic_resistance %||%
    "C:/Users/Google/Documents/HMSC-HIST/derived_inputs/Phanerozoic_DispersalResistance_540_0Ma_5Myr_1deg_topography_v1/topographic_dispersal_resistance_540_0Ma_5Myr_1deg_v1.rds",
  winslash = "/", mustWork = FALSE
)
recipe_path <- normalizePath(
  cfg$recipe %||%
    file.path(previous_case04, "03_hmsc_environment", "case04_modern_0Ma_environment_recipe.rds"),
  winslash = "/", mustWork = FALSE
)
overwrite <- tolower(cfg$overwrite %||% "false") %in% c("true", "t", "1", "yes")

for (p in c(earth_path, topography_path, recipe_path)) {
  if (!file.exists(p)) stop("Missing required input: ", p, call. = FALSE)
}
if (!requireNamespace("HmscEcoEvo", quietly = TRUE)) {
  stop("Install the current HmscEcoEvo package before preparing time slices.", call. = FALSE)
}

recipe <- readRDS(recipe_path)
basis_cols <- setdiff(as.character(recipe$variables), "log_sampling_effort")
earth <- as.data.frame(readRDS(earth_path), stringsAsFactors = FALSE)
earth_required <- c("cell_id", "time_ma", "lon", "lat", "H_state", "land_area_km2", basis_cols)
if (!all(earth_required %in% names(earth))) {
  stop("Palaeo Earth state is missing: ",
       paste(setdiff(earth_required, names(earth)), collapse = ", "), call. = FALSE)
}
earth$H_state <- as.integer(suppressWarnings(as.numeric(earth$H_state)) > 0)
earth <- earth[is.finite(earth$time_ma), earth_required, drop = FALSE]

topography_raw <- readRDS(topography_path)
topography <- if (is.list(topography_raw) && "cells" %in% names(topography_raw)) {
  as.data.frame(topography_raw$cells, stringsAsFactors = FALSE)
} else as.data.frame(topography_raw, stringsAsFactors = FALSE)
topography_required <- c("cell_id", "time_ma", "topographic_resistance",
                         "topographic_permeability", "movement_allowed")
if (!all(topography_required %in% names(topography))) {
  stop("Topographic resistance is missing: ",
       paste(setdiff(topography_required, names(topography)), collapse = ", "), call. = FALSE)
}
topography <- topography[is.finite(topography$time_ma), topography_required,
                         drop = FALSE]

earth_index <- HmscEcoEvo::hee_write_timecube_slices(
  earth, file.path(output, "earth"), keep_cols = earth_required,
  key_cols = "cell_id", prefix = "case04_palaeo_earth_1deg",
  compression = "gzip", overwrite = overwrite
)
topography_index <- HmscEcoEvo::hee_write_timecube_slices(
  topography, file.path(output, "topography"), keep_cols = topography_required,
  key_cols = "cell_id", prefix = "case04_topographic_resistance_1deg",
  compression = "gzip", overwrite = overwrite
)
write_csv(earth_index, file.path(output, "case04_palaeo_earth_slice_index.csv"))
write_csv(topography_index, file.path(output, "case04_topography_slice_index.csv"))

time_check <- merge(
  earth_index[, c("time_ma", "n_rows")],
  topography_index[, c("time_ma", "n_rows")], by = "time_ma", all = TRUE,
  suffixes = c("_earth", "_topography"), sort = FALSE
)
time_check$earth_time_slice <- is.finite(time_check$n_rows_earth)
time_check$topography_time_slice <- is.finite(time_check$n_rows_topography)
time_check$matching_time_slice <- !time_check$earth_time_slice |
  time_check$topography_time_slice
time_check$topography_only_outside_earth_time_domain <-
  !time_check$earth_time_slice & time_check$topography_time_slice
write_csv(time_check, file.path(output, "case04_time_slice_alignment_audit.csv"))
if (!all(time_check$matching_time_slice)) {
  stop("At least one palaeoenvironment time slice lacks matching topographic resistance.",
       call. = FALSE)
}

metadata <- data.frame(
  field = c("earth_source", "earth_md5", "topography_source", "topography_md5",
            "recipe_source", "recipe_md5", "n_time_slices", "n_basis_columns",
            "cell_sampling", "scientific_role"),
  value = c(
    earth_path, unname(tools::md5sum(earth_path)),
    topography_path, unname(tools::md5sum(topography_path)),
    recipe_path, unname(tools::md5sum(recipe_path)),
    nrow(earth_index), length(basis_cols),
    "none_each_output_slice_contains_the_complete_source_1deg_cell_time_grid",
    "streaming_IO_only_no_environmental_rescaling_or_geographic_modification"
  ),
  stringsAsFactors = FALSE
)
write_csv(metadata, file.path(output, "case04_palaeo_slice_metadata.csv"))
writeLines(c(
  "# Case04 1-degree palaeo-grid slices",
  "",
  "Each RDS contains every source cell for one time_ma, retaining the locked HMSC basis variables, coordinates, land/habitat hard state, area, and the corresponding topographic movement fields.",
  "No cells were sampled, rescaled, reprojected or assigned a new region. Dynamic engines may select 40 QA points from a full slice, but occupancy, diversity and GIS maps always use the complete slice.",
  "The slice index is the authoritative time-to-file mapping for streaming Case04/Case05 runs.",
  "Topographic time slices may extend older than the Plant200 tree domain; they are retained for provenance but are not interpreted as biological zeros."
), file.path(output, "README.md"), useBytes = TRUE)

message("Prepared complete 1-degree time slices in: ", output)
