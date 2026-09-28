#!/usr/bin/env Rscript

# Derive the reference-response calibration used by Case05 v3.
#
# This script intentionally does not put the modern HMSC intercept into an
# ancestral niche. It learns a reference on the X beta scale from observed
# modern occupied environments, reconstructs that reference along the dated
# tree, and also reconstructs a positive response-contrast scale. The latter
# keeps the P1 covariate comparable across time: re-standardising X beta within
# each palaeo slice would erase genuine global deterioration or improvement in
# environmental support. This script records the unavoidable data-reuse
# boundary: the current Plant200 HMSC posterior and occurrence matrix originate
# from the same survey source.

or_null <- function(x, y) {
  if (is.null(x) || !length(x) || identical(x, "") || is.na(x)) y else x
}
parse_args <- function(args) {
  out <- list()
  for (arg in args) if (grepl("^--", arg)) {
    z <- sub("^--", "", arg)
    at <- regexpr("=", z, fixed = TRUE)
    if (at < 0L) out[[z]] <- TRUE else {
      out[[substr(z, 1L, at - 1L)]] <- substr(z, at + 1L, nchar(z))
    }
  }
  out
}
safe_dir <- function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE)
  normalizePath(x, winslash = "/", mustWork = FALSE)
}
write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("data.table", quietly = TRUE)) data.table::fwrite(x, path) else utils::write.csv(x, path, row.names = FALSE)
  invisible(path)
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
file_arg <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
if (!length(script)) script <- "scripts/case05_v3_prepare_reference_support.R"
pkg_root <- normalizePath(file.path(dirname(normalizePath(script[[1L]], winslash = "/", mustWork = FALSE)), ".."),
                          winslash = "/", mustWork = TRUE)
workspace <- normalizePath(file.path(pkg_root, ".."), winslash = "/", mustWork = TRUE)
beta_path <- normalizePath(or_null(cfg$beta_path, file.path(
  pkg_root, "derived_inputs", "case04_plant200_hmsc_beta_posterior_25draws_20260919",
  "case04_hmsc_beta_posterior_25draws.csv"
)), winslash = "/", mustWork = TRUE)
earth_index_path <- normalizePath(or_null(cfg$earth_index, file.path(
  pkg_root, "derived_inputs", "case04_plant200_palaeo_1deg_slices_20260919",
  "case04_palaeo_earth_slice_index.csv"
)), winslash = "/", mustWork = TRUE)
plant_input_dir <- normalizePath(or_null(cfg$plant_input_dir, file.path(
  workspace, "HmscEcoEvo_API8_20260818", "data_external",
  "plant_200_multifamily_extant_20260820", "prepared_inputs",
  "case03_data_raw", "Plant200_multifamily_extant_1deg_traits"
)), winslash = "/", mustWork = TRUE)
tree_path <- normalizePath(or_null(cfg$tree, file.path(plant_input_dir, "tree.tre")),
                           winslash = "/", mustWork = TRUE)
output <- safe_dir(or_null(cfg$output, file.path(
  pkg_root, "derived_inputs", "case05_v3_reference_environmental_support"
)))
presence_quantile <- suppressWarnings(as.numeric(or_null(cfg$presence_quantile, 0.1)))
if (!is.finite(presence_quantile) || presence_quantile < 0 || presence_quantile > 1) {
  stop("--presence_quantile must be in [0, 1].", call. = FALSE)
}
threshold_method <- match.arg(
  tolower(or_null(cfg$threshold_method, "prevalence_matched")),
  c("presence_quantile", "prevalence_matched")
)
if (!requireNamespace("pkgload", quietly = TRUE) || !requireNamespace("ape", quietly = TRUE)) {
  stop("pkgload and ape are required.", call. = FALSE)
}
pkgload::load_all(pkg_root, quiet = TRUE)

beta <- utils::read.csv(beta_path, check.names = FALSE, stringsAsFactors = FALSE)
need_beta <- c("lineage", "response_draw", "intercept", "hmsc_intercept_original")
if (!all(need_beta %in% names(beta))) stop("Beta posterior table is missing required columns.", call. = FALSE)
axes <- setdiff(names(beta)[vapply(beta, is.numeric, logical(1))],
                c("chain", "iteration", "response_weight", "intercept",
                  "hmsc_intercept_original"))
if (!length(axes)) stop("No HMSC environmental response axes were found.", call. = FALSE)

y <- as.matrix(utils::read.csv(file.path(plant_input_dir, "comm.csv"),
                               row.names = 1, check.names = FALSE))
storage.mode(y) <- "double"
if (!setequal(colnames(y), unique(beta$lineage))) {
  stop("Plant200 comm.csv species do not match HMSC posterior lineages.", call. = FALSE)
}
y <- y[, sort(colnames(y)), drop = FALSE]
sites <- utils::read.csv(file.path(plant_input_dir, "sites.csv"), stringsAsFactors = FALSE)
if (!all(c("site_id", "lon", "lat") %in% names(sites))) {
  stop("sites.csv lacks site_id, lon, or lat.", call. = FALSE)
}
if (!setequal(rownames(y), sites$site_id)) stop("comm.csv row names must match sites.csv site_id.", call. = FALSE)
sites <- sites[match(rownames(y), sites$site_id), , drop = FALSE]
earth_index <- utils::read.csv(earth_index_path, stringsAsFactors = FALSE)
earth <- as.data.frame(readRDS(earth_index$file[[match(0, earth_index$time_ma)]]))
if (!all(axes %in% names(earth))) {
  stop("0 Ma Earth grid lacks HMSC basis axes: ", paste(setdiff(axes, names(earth)), collapse = ", "), call. = FALSE)
}
earth_key <- paste(formatC(earth$lon, format = "f", digits = 6),
                   formatC(earth$lat, format = "f", digits = 6), sep = "|")
site_key <- paste(formatC(sites$lon, format = "f", digits = 6),
                  formatC(sites$lat, format = "f", digits = 6), sep = "|")
hit <- match(site_key, earth_key)
if (anyNA(hit)) {
  stop("Some HMSC training coordinates were not found in the complete 0 Ma Earth grid.", call. = FALSE)
}
X <- as.matrix(earth[hit, axes, drop = FALSE])
rownames(X) <- rownames(y)

draws <- unique(as.character(beta$response_draw))
threshold_rows <- vector("list", length(draws))
for (ii in seq_along(draws)) {
  draw <- draws[[ii]]
  b <- beta[beta$response_draw == draw, c("lineage", axes), drop = FALSE]
  b <- b[match(colnames(y), b$lineage), , drop = FALSE]
  eta <- X %*% t(as.matrix(b[, axes, drop = FALSE]))
  colnames(eta) <- b$lineage
  ref <- HmscEcoEvo::hee_environmental_filtering_reference_threshold(
    y[, b$lineage, drop = FALSE], eta,
    presence_quantile = presence_quantile,
    method = threshold_method
  )
  # A robust modern response width puts the reference contrast on a common,
  # lineage-specific scale. It is a P1 calibration quantity, not an HMSC
  # intercept and not a historical occurrence probability.
  eta_scale <- vapply(seq_len(ncol(eta)), function(j) {
    z <- eta[, j]
    z <- z[is.finite(z)]
    value <- if (length(z) > 1L) stats::IQR(z) / 1.349 else NA_real_
    if (!is.finite(value) || value <= 1e-8) {
      value <- if (length(z) > 1L) stats::sd(z) else NA_real_
    }
    if (!is.finite(value) || value <= 1e-8) value <- 1
    value
  }, numeric(1))
  ref$eta_scale <- unname(eta_scale[match(ref$lineage, colnames(eta))])
  ref$response_draw <- draw
  threshold_rows[[ii]] <- ref
}
tip_reference <- do.call(rbind, threshold_rows)
if (any(!is.finite(tip_reference$eta_reference))) {
  stop("At least one lineage has no observed modern occurrence after site matching.", call. = FALSE)
}

tree <- ape::read.tree(tree_path)
earth_times <- sort(unique(suppressWarnings(as.numeric(earth_index$time_ma))), decreasing = TRUE)
root_age <- max(ape::node.depth.edgelength(tree)[seq_along(tree$tip.label)])
node_age <- root_age - ape::node.depth.edgelength(tree)
node_times <- node_age[(length(tree$tip.label) + 1L):(length(tree$tip.label) + tree$Nnode)]
node_times <- node_times[node_times <= max(earth_times) + 1e-8 & node_times >= min(earth_times) - 1e-8]
# Internal-node responses are conventionally reported on the incoming branch.
# Add a tiny younger-side query so daughter thresholds are available immediately
# after a dated speciation event without moving the Earth slice.
eps <- min(1e-4, min(tree$edge.length[tree$edge.length > 0]) / 20)
post_node_times <- node_times - eps
post_node_times <- post_node_times[post_node_times >= min(earth_times) - 1e-8]
times <- sort(unique(c(earth_times, node_times, post_node_times)), decreasing = TRUE)
ref_input <- data.frame(
  species = tip_reference$lineage,
  response_draw = tip_reference$response_draw,
  intercept = 0,
  eta_reference = tip_reference$eta_reference,
  log_eta_scale = log(tip_reference$eta_scale),
  stringsAsFactors = FALSE
)
ancestral_reference <- HmscEcoEvo::hee_evolution_ancestral_response_direct(
  tip_beta_draws = ref_input, tree = tree, times = times,
  basis_cols = c("eta_reference", "log_eta_scale"), reconstruct_intercept = FALSE,
  branch_uncertainty = "none"
)
ancestral_reference$eta_scale <- exp(ancestral_reference$log_eta_scale)
if (any(!is.finite(ancestral_reference$eta_scale) |
        ancestral_reference$eta_scale <= 1e-8)) {
  stop("Ancestral response-scale reconstruction produced an invalid eta_scale.",
       call. = FALSE)
}
write_csv(tip_reference, file.path(output, "tip_response_reference_thresholds.csv"))
write_csv(ancestral_reference, file.path(output, "ancestral_response_reference_thresholds.csv"))
write_csv(data.frame(
  field = c("case", "presence_quantile", "n_hmsc_response_draws",
            "threshold_method",
            "n_training_sites", "n_tip_lineages", "n_tree_time_nodes",
            "post_speciation_response_epsilon_myr", "modern_intercept_role",
            "response_scale_role",
            "scientific_boundary"),
  value = c(
    "case05_v3_reference_environmental_support", presence_quantile,
    length(draws), threshold_method, nrow(y), ncol(y), length(node_times), eps,
    "excluded_from_ancestral_environmental_response; only X beta response shape is used",
    "eta_scale = robust IQR(X beta)/1.349 on the modern training grid; log scale is reconstructed along the dated tree and prevents time-slice re-centering of the P1 process covariate",
    "Reference thresholds reuse Plant200 modern occurrence data used to fit the current HMSC posterior. They calibrate relative environmental support, not an independent endpoint likelihood or historical occupancy probability."
  ),
  stringsAsFactors = FALSE
), file.path(output, "reference_support_metadata.csv"))
write_csv(data.frame(
  site_id = rownames(y), earth_cell_id = earth$cell_id[hit],
  lon = earth$lon[hit], lat = earth$lat[hit],
  match_source = "sites_csv_exact_lon_lat_to_complete_0Ma_earth_grid",
  stringsAsFactors = FALSE
), file.path(output, "modern_training_site_to_earth_grid.csv"))
message("Case05 v3 response-reference preparation complete: ", output)
