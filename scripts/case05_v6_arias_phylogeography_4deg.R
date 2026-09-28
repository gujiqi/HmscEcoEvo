#!/usr/bin/env Rscript

# Case05 v6: explicit-grid phylogeography, with 4-degree map aggregation ----
#
# This runner deliberately replaces the v5 demographic scenario layer.  It
# does NOT calculate colonisation, establishment, persistence, local
# extinction, occupancy probability, or species richness.  Instead it follows
# the landscape-explicit phylogeographic design of Arias (2024): a stable
# plate-carried grid defines location states, a source-normalised landscape
# diffusion process moves location probability, and Felsenstein pruning
# conditions ancestral locations on present-day occurrence densities and a
# dated tree.  The Flannery-Sutherland et al. (2025) landscape graph remains a
# diagnostic concept; it is not a second probability multiplier.
#
# The 1-degree carrier graph is retained internally because its stable
# `track_id` separates plate carriage from biological diffusion.  Results are
# aggregated to 4-degree cells only for fast, readable maps.  Collapsing the
# carrier process before propagation would incorrectly turn plate motion into
# apparent dispersal.

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x)[1L]) y else x

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) next
    kv <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    out[[kv[[1L]]]] <- paste(kv[-1L], collapse = "=")
  }
  out
}

as_integer <- function(x, default) {
  value <- suppressWarnings(as.integer(x %||% default))
  if (!is.finite(value) || value < 1L) stop("Expected a positive integer.", call. = FALSE)
  value
}

safe_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
  invisible(path)
}

normalise <- function(x, label) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x) | x < 0] <- 0
  total <- sum(x)
  if (!is.finite(total) || total <= 0) stop(label, " has no positive mass.", call. = FALSE)
  x / total
}

child_key <- function(tree, child) {
  if (child <= length(tree$tip.label)) tree$tip.label[[child]] else paste0("node_", child)
}

tree_geometry <- function(tree) {
  depth <- ape::node.depth.edgelength(tree)
  root <- setdiff(tree$edge[, 1L], tree$edge[, 2L])
  if (length(root) != 1L) stop("Tree must have exactly one root.", call. = FALSE)
  root_age <- max(depth[seq_along(tree$tip.label)])
  list(root = root, root_age = root_age, node_age = root_age - depth)
}

bin4 <- function(lon, lat) {
  ix <- pmin(89L, pmax(0L, floor((as.numeric(lon) + 180) / 4)))
  iy <- pmin(44L, pmax(0L, floor((as.numeric(lat) + 90) / 4)))
  as.integer(iy * 90L + ix + 1L)
}

make_map_grid <- function() {
  ix <- rep(0:89, times = 45)
  iy <- rep(0:44, each = 90)
  data.frame(
    bin = seq_len(90L * 45L),
    lon = -178 + 4 * ix,
    lat = -88 + 4 * iy,
    stringsAsFactors = FALSE
  )
}

named_logical <- function(carrier, states, field = "active_land") {
  result <- stats::setNames(rep(FALSE, length(states)), states)
  keep <- match(carrier$track_id, states)
  value <- if (field %in% names(carrier)) carrier[[field]] else FALSE
  result[keep[!is.na(keep)]] <- isTRUE(value) | (!is.na(value) & value > 0)
  result
}

empty_map <- function() rep(0, 90L * 45L)

map_from_state <- function(value, slice, states) {
  result <- empty_map()
  active <- slice$active_by_state
  bin <- slice$bin_by_state
  keep <- active & is.finite(bin) & is.finite(value) & value != 0
  if (any(keep)) {
    summed <- rowsum(value[keep], group = bin[keep], reorder = FALSE)
    result[as.integer(rownames(summed))] <- summed[, 1L]
  }
  result
}

map_land <- function(slice) {
  result <- logical(90L * 45L)
  keep <- slice$active_by_state & is.finite(slice$bin_by_state)
  if (any(keep)) result[unique(slice$bin_by_state[keep])] <- TRUE
  result
}

write_map_artifacts <- function(map_grid, value, land, time_ma, metric, limits,
                                dirs, unit, source, subtitle) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("terra", quietly = TRUE)) {
    stop("Case05 v6 map export requires ggplot2 and terra.", call. = FALSE)
  }
  tab <- map_grid
  tab$value <- as.numeric(value)
  tab$land_available <- as.logical(land)
  tab$value[!tab$land_available] <- NA_real_
  png_path <- file.path(dirs$maps, metric,
                        sprintf("%s_%03dMa.png", metric, round(time_ma)))
  tif_path <- file.path(dirs$maps, metric,
                        sprintf("%s_%03dMa.tif", metric, round(time_ma)))
  safe_dir(dirname(png_path))
  colour_scale <- if (limits[[1L]] < 0) {
    ggplot2::scale_fill_gradient2(
      low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0,
      limits = limits, oob = scales::squish, na.value = "grey80", name = unit
    )
  } else {
    ggplot2::scale_fill_viridis_c(
      option = "C", limits = limits, oob = scales::squish,
      na.value = "grey80", name = unit
    )
  }
  plot <- ggplot2::ggplot(tab, ggplot2::aes(x = lon, y = lat)) +
    ggplot2::geom_tile(data = tab[!tab$land_available, , drop = FALSE],
                       fill = "grey80", colour = NA, show.legend = FALSE) +
    ggplot2::geom_tile(data = tab[tab$land_available, , drop = FALSE],
                       ggplot2::aes(fill = value), colour = NA) +
    colour_scale +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::labs(
      title = paste0("Case05 v6: ", metric, "; ", round(time_ma), " Ma"),
      subtitle = subtitle,
      x = "Palaeographic longitude (degrees)",
      y = "Palaeographic latitude (degrees)",
      caption = paste0("Grey = ocean or unavailable palaeogeography. ", source)
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_line(colour = "grey88"))
  ggplot2::ggsave(png_path, plot, width = 12, height = 6.4, dpi = 150)

  raster <- terra::rast(ncols = 90, nrows = 45, xmin = -180, xmax = 180,
                         ymin = -90, ymax = 90, crs = "EPSG:4326")
  cells <- terra::cellFromXY(raster, as.matrix(tab[, c("lon", "lat")]))
  values <- rep(NA_real_, terra::ncell(raster))
  values[cells] <- tab$value
  terra::values(raster) <- values
  names(raster) <- metric
  terra::writeRaster(raster, tif_path, overwrite = TRUE,
                     gdal = c("COMPRESS=DEFLATE"))
  data.frame(
    metric = metric,
    time_ma = time_ma,
    png = normalizePath(png_path, winslash = "/", mustWork = TRUE),
    geotiff = normalizePath(tif_path, winslash = "/", mustWork = TRUE),
    unit = unit,
    scale_min = limits[[1L]],
    scale_max = limits[[2L]],
    source = source,
    stringsAsFactors = FALSE
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
pkg_root <- normalizePath(args$pkg_root %||%
  "C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo", winslash = "/", mustWork = TRUE)
output <- args$output %||% file.path(
  "C:/Users/Google/Documents/HMSC-HIST/outputs/HmscEcoEvo",
  "case05_plant200_v6_arias_phygeography_4deg_preflight_20260923"
)
output <- safe_dir(output)
n_tips <- as_integer(args$n_tips, 20L)
n_draws <- as_integer(args$n_response_draws, 5L)
seed <- as_integer(args$seed, 20260923L)

if (!requireNamespace("HmscEcoEvo", quietly = TRUE)) {
  stop("Install the current HmscEcoEvo package before running Case05 v6.", call. = FALSE)
}
for (package in c("ape", "Matrix", "ggplot2", "terra")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Case05 v6 requires package `", package, "`.", call. = FALSE)
  }
}

dirs <- list(
  config = safe_dir(file.path(output, "00_contract")),
  audit = safe_dir(file.path(output, "01_input_audit")),
  response = safe_dir(file.path(output, "02_ancestral_response")),
  pruning = safe_dir(file.path(output, "03_phylogeographic_pruning")),
  maps = safe_dir(file.path(output, "04_shared_scale_maps")),
  quality = safe_dir(file.path(output, "05_quality"))
)

contract <- data.frame(
  item = c(
    "analysis", "internal_geographic_state", "map_aggregation",
    "environmental_response", "diffusion", "terminal_constraint",
    "colonisation", "persistence", "local_extinction", "occupancy_probability",
    "richness", "HMSC_intercept", "time_ma"
  ),
  status = c(
    "landscape_explicit_phylogeographic_reconstruction",
    "stable_1_degree_plate_carrier_track_id",
    "4_degree_display_and_summary_only",
    "separate_ancestral_HMSC_linear_predictor_diagnostic",
    "Arias_source_normalised_location_transition",
    "Felsenstein_pruning_conditioned_on_modern_occurrence_density",
    "NOT_RUN_NOT_IDENTIFIABLE_FROM_CURRENT_DATA",
    "NOT_RUN_NOT_IDENTIFIABLE_FROM_CURRENT_DATA",
    "NOT_RUN_NOT_INFERRED",
    "NOT_REPORTED",
    "NOT_REPORTED",
    "excluded_to_avoid_modern_prevalence_as_ancestral_niche",
    "larger_Ma_values_are_older"
  ),
  explanation = c(
    "Arias-style dynamic-grid location reconstruction, not a demographic occupancy model.",
    "Plate carriage follows a stable carrier identity; active diffusion follows the landscape kernel.",
    "Native carrier transitions remain at 1 degree; maps are binned at 4 degrees after inference.",
    "HMSC posterior beta is reconstructed on the dated tree and evaluated against paleoenvironment; it is not multiplied into diffusion.",
    "Location movement and retention only; no establishment or survival coefficient.",
    "Modern occurrence densities provide a modular terminal likelihood for tree pruning.",
    "Removed from the empirical Case05 core; old scenario outputs remain legacy-only.",
    "Removed from the empirical Case05 core; old scenario outputs remain legacy-only.",
    "Geographic sink tracks unavailable land, not biological extinction.",
    "No cell occupancy probability is claimed by this runner.",
    "No richness, Shannon, Simpson, PD, or functional-diversity map is claimed without a calibrated range/occupancy layer.",
    "Modern HMSC intercept reflects prevalence/sampling and is not reconstructed as ancestral niche.",
    "The time axis runs from older (larger Ma) toward 0 Ma (present)."
  ),
  stringsAsFactors = FALSE
)
write_csv(contract, file.path(dirs$config, "case05_v6_scientific_contract.csv"))

input_root <- normalizePath(file.path(
  dirname(pkg_root), "HmscEcoEvo_API8_20260818", "data_external",
  "plant_200_multifamily_extant_20260820", "prepared_inputs", "case03_data_raw",
  "Plant200_multifamily_extant_1deg_traits"
), winslash = "/", mustWork = TRUE)
legacy_output <- normalizePath(file.path(
  dirname(pkg_root), "outputs", "HmscEcoEvo",
  "case05_plant200_4deg_effortaware_formal25_all66_20260923"
), winslash = "/", mustWork = TRUE)
carrier_root <- normalizePath(file.path(
  pkg_root, "derived_inputs", "case05_v2_paleomap_h3_carrier_325_0Ma_1deg",
  "movement_cache_arias_topographic_rate0p4"
), winslash = "/", mustWork = TRUE)

comm <- utils::read.csv(file.path(input_root, "comm.csv"), check.names = FALSE,
                        stringsAsFactors = FALSE)
sites <- utils::read.csv(file.path(input_root, "sites.csv"), check.names = FALSE,
                         stringsAsFactors = FALSE)
tree <- ape::read.tree(file.path(input_root, "tree.tre"))
beta <- readRDS(file.path(legacy_output, "02_hmsc_4deg",
                          "hmsc_environmental_beta_posterior_draws.rds"))
carrier_index <- utils::read.csv(file.path(carrier_root,
  "case05_v2_movement_cache_index.csv"), stringsAsFactors = FALSE)

if (!identical(comm$site_id, sites$site_id)) {
  stop("comm.csv and sites.csv site_id order differs; Case05 v6 refuses positional matching.", call. = FALSE)
}
species <- setdiff(names(comm), "site_id")
missing <- setdiff(species, tree$tip.label)
if (length(missing)) stop("Tree lacks comm species: ", paste(missing, collapse = ", "), call. = FALSE)
if (!all(c("lineage", "response_draw") %in% names(beta))) {
  stop("Stored HMSC beta posterior has an unexpected schema.", call. = FALSE)
}
axes <- intersect(
  c("MAT_pohl_C", "MAP_pohl_mm_yr", "elevation_m", "T_seasonality_pohl_sd",
    "P_seasonality_pohl_sd", "moisture_availability_index_z",
    "wetland_potential_index_z", "relief_3x3_m"), names(beta)
)
if (length(axes) < 2L) stop("Stored HMSC beta posterior lacks expected environmental axes.", call. = FALSE)

occurrences <- colSums(as.matrix(comm[, species, drop = FALSE]) > 0)
selected_tips <- names(sort(occurrences, decreasing = TRUE))[seq_len(min(n_tips, length(species)))]
tree <- ape::keep.tip(tree, selected_tips)
tree <- ape::reorder.phylo(tree, "cladewise")
geom <- tree_geometry(tree)
draw_ids <- unique(beta$response_draw)
draw_ids <- draw_ids[seq_len(min(n_draws, length(draw_ids)))]
tip_beta <- beta[beta$lineage %in% tree$tip.label & beta$response_draw %in% draw_ids,
                 c("lineage", "response_draw", axes), drop = FALSE]
names(tip_beta)[names(tip_beta) == "lineage"] <- "species"
if (nrow(tip_beta) != length(tree$tip.label) * length(draw_ids)) {
  stop("HMSC beta posterior does not cover every selected tip and draw.", call. = FALSE)
}

input_audit <- data.frame(
  item = c("HmscEcoEvo_version", "comm_sites", "tree_tips", "selected_tips",
           "HMSC_response_draws", "carrier_time_slices", "environment_axes",
           "tree_root_age_ma", "output_map_grid"),
  value = c(
    as.character(utils::packageVersion("HmscEcoEvo")), nrow(comm), length(tree$tip.label),
    length(selected_tips), length(draw_ids), nrow(carrier_index),
    paste(axes, collapse = ";"), format(geom$root_age, digits = 7), "4 degrees"
  ),
  stringsAsFactors = FALSE
)
write_csv(input_audit, file.path(dirs$audit, "input_audit.csv"))
write_csv(data.frame(species = names(occurrences), occurrences = as.integer(occurrences),
                     selected_for_preflight = names(occurrences) %in% selected_tips,
                     stringsAsFactors = FALSE),
          file.path(dirs$audit, "species_selection.csv"))

# Build the stable carrier state space first.  It is deliberately not binned
# before inference: track identity carries a population with its plate.
carrier_index <- carrier_index[order(carrier_index$time_ma, decreasing = TRUE), , drop = FALSE]
all_tracks <- character()
for (path in carrier_index$file) {
  all_tracks <- union(all_tracks, as.character(readRDS(path)$carrier$track_id))
}
states <- sort(unique(all_tracks))
sink_state <- "__outside_available_geography__"
if (sink_state %in% states) stop("Reserved sink state collides with track IDs.", call. = FALSE)

slices <- list()
for (i in seq_len(nrow(carrier_index))) {
  object <- readRDS(carrier_index$file[[i]])
  carrier <- object$carrier
  rownames(carrier) <- carrier$track_id
  carrier <- carrier[match(states, carrier$track_id), , drop = FALSE]
  active <- !is.na(carrier$track_id) & isTRUE(carrier$active_land) |
    (!is.na(carrier$active_land) & carrier$active_land > 0)
  bins <- rep(NA_integer_, length(states))
  valid_coord <- !is.na(carrier$paleo_lon) & !is.na(carrier$paleo_lat)
  bins[valid_coord] <- bin4(carrier$paleo_lon[valid_coord], carrier$paleo_lat[valid_coord])
  slices[[as.character(carrier_index$time_ma[[i]])]] <- list(
    time_ma = carrier_index$time_ma[[i]], carrier = carrier,
    active_by_state = active, bin_by_state = bins
  )
}

# Each older-to-younger stage gets one Arias-style probability transition.
# A carrier that becomes geographically unavailable exports through the old
# landscape before residual mass enters the explicit geographic sink.
transition_cache <- new.env(parent = emptyenv())
make_transition <- function(interval_index, delta_t) {
  key <- paste(interval_index, format(delta_t, scientific = FALSE, trim = TRUE), sep = "__")
  if (exists(key, envir = transition_cache, inherits = FALSE)) {
    return(get(key, envir = transition_cache, inherits = FALSE))
  }
  older <- slices[[as.character(carrier_index$time_ma[[interval_index]])]]
  younger <- slices[[as.character(carrier_index$time_ma[[interval_index + 1L]])]]
  movement <- readRDS(carrier_index$file[[interval_index]])$movement
  transition <- HmscEcoEvo::hee_dispersal_location_transition(
    movement_kernel = movement, cell_order = states, delta_t_myr = delta_t,
    source_available_cells = stats::setNames(older$active_by_state, states),
    target_available_cells = stats::setNames(younger$active_by_state, states),
    sink_state = sink_state, sparse = TRUE
  )
  record <- list(transition = transition$transition,
                 time_ma = younger$time_ma,
                 delta_t_myr = delta_t,
                 diagnostics = transition$diagnostics)
  assign(key, record, envir = transition_cache)
  record
}

# Only complete palaeoenvironment intervals are used within a branch.  The
# remaining sub-interval at an exact dated node is identity, rather than
# pretending that a 5-Myr geology layer resolves a shorter interval.
identity_transition <- function(time_ma) {
  n <- length(states) + 1L
  P <- Matrix::Diagonal(n = n, x = 1)
  dimnames(P) <- list(c(states, sink_state), c(states, sink_state))
  list(transition = P, time_ma = time_ma, delta_t_myr = 0,
       diagnostics = data.frame(identity_stage = TRUE, stringsAsFactors = FALSE))
}

interval_old <- carrier_index$time_ma[-nrow(carrier_index)]
interval_young <- carrier_index$time_ma[-1L]
branch_steps <- list()
for (edge_index in seq_len(nrow(tree$edge))) {
  parent <- tree$edge[edge_index, 1L]
  child <- tree$edge[edge_index, 2L]
  parent_age <- geom$node_age[[parent]]
  child_age <- geom$node_age[[child]]
  included <- which(interval_old <= parent_age + 1e-8 & interval_young >= child_age - 1e-8)
  steps <- lapply(included, function(i) make_transition(i, interval_old[[i]] - interval_young[[i]]))
  if (!length(steps) || abs(tail(vapply(steps, `[[`, numeric(1), "time_ma"), 1L) - child_age) > 1e-8) {
    steps[[length(steps) + 1L]] <- identity_transition(child_age)
  }
  branch_steps[[child_key(tree, child)]] <- steps
}

transition_qc <- do.call(rbind, lapply(ls(transition_cache), function(key) {
  z <- get(key, envir = transition_cache)
  cbind(cache_key = key, z$diagnostics, stringsAsFactors = FALSE)
}))
write_csv(as.data.frame(transition_qc), file.path(dirs$pruning, "transition_quality.csv"))

# Construct a terminal occurrence-density likelihood.  These are presences
# only, not absence likelihoods.  They use the same modern observations as the
# upstream HMSC and are therefore explicitly labelled as a modular/cut model.
terminal_slice <- slices[["0"]]
site_bin <- bin4(sites$lon, sites$lat)
tip_likelihood <- list()
terminal_audit <- list()
for (tip in tree$tip.label) {
  record_bins <- unique(site_bin[as.matrix(comm[, tip, drop = FALSE])[, 1L] > 0])
  hit <- terminal_slice$active_by_state & terminal_slice$bin_by_state %in% record_bins
  likelihood <- stats::setNames(rep(0, length(states) + 1L), c(states, sink_state))
  if (any(hit)) {
    area <- terminal_slice$carrier$land_area_km2[hit]
    likelihood[states[hit]] <- normalise(area, paste0("terminal likelihood for ", tip))
  } else {
    stop("No active 0 Ma carrier state matches the 4-degree occurrence bins for ", tip,
         ".", call. = FALSE)
  }
  tip_likelihood[[tip]] <- likelihood
  terminal_audit[[tip]] <- data.frame(
    species = tip, n_occurrence_sites = sum(comm[[tip]] > 0),
    n_occurrence_bins_4deg = length(record_bins), n_terminal_carrier_states = sum(hit),
    stringsAsFactors = FALSE
  )
}
terminal_audit <- do.call(rbind, terminal_audit)
write_csv(terminal_audit, file.path(dirs$audit, "terminal_occurrence_density_audit.csv"))

root_slice_time <- carrier_index$time_ma[[which.min(abs(carrier_index$time_ma - geom$root_age))]]
root_slice <- slices[[as.character(root_slice_time)]]
root_weight <- root_slice$carrier$land_area_km2 * root_slice$carrier$topographic_permeability
root_weight[!root_slice$active_by_state] <- 0
root_prior <- stats::setNames(c(normalise(root_weight, "root spatial prior"), 0),
                              c(states, sink_state))

analysis_times <- carrier_index$time_ma[carrier_index$time_ma <= geom$root_age + 1e-8]
responses <- HmscEcoEvo::hee_evolution_ancestral_response_direct(
  tip_beta_draws = tip_beta, tree = tree, times = analysis_times,
  basis_cols = axes, reconstruct_intercept = FALSE,
  branch_uncertainty = "brownian_bridge", seed = seed
)
saveRDS(responses, file.path(dirs$response, "ancestral_beta_response_draws.rds"), compress = "gzip")
response_summary <- aggregate(responses[, axes, drop = FALSE],
                              responses[, c("lineage", "time_ma"), drop = FALSE], mean)
write_csv(response_summary, file.path(dirs$response, "ancestral_beta_response_mean.csv"))

map_grid <- make_map_grid()
times_for_maps <- as.character(sort(unique(c(analysis_times, root_slice_time)), decreasing = TRUE))
density_sum <- stats::setNames(lapply(times_for_maps, function(x) empty_map()), times_for_maps)
eta_sum <- stats::setNames(lapply(times_for_maps, function(x) empty_map()), times_for_maps)
eta_weight <- stats::setNames(lapply(times_for_maps, function(x) empty_map()), times_for_maps)
callback_count <- stats::setNames(rep(0L, length(times_for_maps)), times_for_maps)
per_draw_loglik <- numeric(length(draw_ids)); names(per_draw_loglik) <- draw_ids

response_at <- function(draw_id, lineage, time_ma) {
  row <- responses[responses$response_draw == draw_id & responses$lineage == lineage &
                     abs(responses$time_ma - time_ma) < 1e-8, , drop = FALSE]
  if (!nrow(row)) return(NULL)
  as.numeric(row[1L, axes, drop = TRUE])
}

for (draw_id in draw_ids) {
  callback <- function(branch_key, child_node, time_ma, stage_index, posterior) {
    time_key <- as.character(time_ma)
    if (!time_key %in% names(density_sum)) return(invisible(NULL))
    slice <- slices[[time_key]]
    state_mass <- as.numeric(posterior[states])
    density_sum[[time_key]] <<- density_sum[[time_key]] + map_from_state(state_mass, slice, states)
    callback_count[[time_key]] <<- callback_count[[time_key]] + 1L
    beta_row <- response_at(draw_id, branch_key, time_ma)
    if (!is.null(beta_row)) {
      environment <- as.matrix(slice$carrier[, axes, drop = FALSE])
      eta <- as.vector(environment %*% beta_row)
      valid <- slice$active_by_state & is.finite(eta)
      eta_sum[[time_key]] <<- eta_sum[[time_key]] + map_from_state(state_mass * eta, slice, states)
      eta_weight[[time_key]] <<- eta_weight[[time_key]] + map_from_state(state_mass * valid, slice, states)
    }
    invisible(NULL)
  }
  result <- HmscEcoEvo::hee_dispersal_time_ordered_pruning(
    tree = tree, tip_likelihood = tip_likelihood, branch_steps = branch_steps,
    root_prior = root_prior, state_callback = callback
  )
  per_draw_loglik[[draw_id]] <- result$log_likelihood
  if (!is.finite(result$log_likelihood)) {
    stop("Phylogeographic pruning returned non-finite likelihood for draw ", draw_id,
         ".", call. = FALSE)
  }
  root_density <- map_from_state(as.numeric(result$root_posterior[states]), root_slice, states)
  density_sum[[as.character(root_slice_time)]] <- density_sum[[as.character(root_slice_time)]] + root_density
  callback_count[[as.character(root_slice_time)]] <- callback_count[[as.character(root_slice_time)]] + 1L
  saveRDS(result, file.path(dirs$pruning, paste0("pruning_", draw_id, ".rds")), compress = "gzip")
}
write_csv(data.frame(response_draw = names(per_draw_loglik), log_likelihood = per_draw_loglik,
                     stringsAsFactors = FALSE),
          file.path(dirs$pruning, "pruning_log_likelihoods.csv"))

maps <- list()
mm <- 0L
mean_density <- lapply(density_sum, function(x) x / length(draw_ids))
density_limit <- c(0, max(unlist(mean_density), na.rm = TRUE))
eta_mean <- stats::setNames(lapply(names(eta_sum), function(key) {
  answer <- eta_sum[[key]] / pmax(eta_weight[[key]], .Machine$double.eps)
  answer[eta_weight[[key]] <= 0] <- NA_real_
  answer
}), names(eta_sum))
eta_limit_value <- max(abs(unlist(eta_mean)), na.rm = TRUE)
eta_limit <- c(-eta_limit_value, eta_limit_value)
for (time_key in names(mean_density)) {
  slice <- slices[[time_key]]
  if (is.null(slice)) next
  land <- map_land(slice)
  mm <- mm + 1L
  maps[[mm]] <- write_map_artifacts(
    map_grid, mean_density[[time_key]], land, as.numeric(time_key),
    metric = "posterior_lineage_location_density", limits = density_limit, dirs = dirs,
    unit = "expected selected-lineage location density",
    source = "Arias-style diffusion + terminal occurrence pruning; not occupancy or richness.",
    subtitle = "4-degree aggregation after native 1-degree plate-carried diffusion; common scale across time"
  )
  mm <- mm + 1L
  maps[[mm]] <- write_map_artifacts(
    map_grid, eta_mean[[time_key]], land, as.numeric(time_key),
    metric = "posterior_weighted_ancestral_HMSC_linear_predictor", limits = eta_limit, dirs = dirs,
    unit = "posterior-location-weighted environmental linear predictor",
    source = "Separate ancestral HMSC beta diagnostic; never multiplied into diffusion.",
    subtitle = "Intercept excluded; larger absolute values indicate stronger modelled environmental association"
  )
}
map_index <- do.call(rbind, maps)
write_csv(map_index, file.path(dirs$maps, "case05_v6_map_index.csv"))

transition_error <- max(transition_qc$max_source_normalisation_error, na.rm = TRUE)
quality <- data.frame(
  check = c("colonisation_not_run", "persistence_not_run", "no_occupancy_or_richness_claim",
            "all_pruning_draws_finite", "source_normalisation_error", "native_grid", "map_grid",
            "n_selected_tips", "n_response_draws", "n_png_maps", "n_geotiff_maps"),
  value = c(TRUE, TRUE, TRUE, all(is.finite(per_draw_loglik)), transition_error,
            "1_degree_plate_carrier", "4_degree_aggregation", length(tree$tip.label),
            length(draw_ids), sum(file.exists(map_index$png)), sum(file.exists(map_index$geotiff))),
  pass = c(TRUE, TRUE, TRUE, all(is.finite(per_draw_loglik)), transition_error <= 1e-10,
           TRUE, TRUE, TRUE, TRUE, all(file.exists(map_index$png)), all(file.exists(map_index$geotiff))),
  stringsAsFactors = FALSE
)
write_csv(quality, file.path(dirs$quality, "case05_v6_quality_gates.csv"))
saveRDS(list(
  density = mean_density, environmental_linear_predictor = eta_mean,
  active_branch_callbacks = callback_count, map_times_ma = as.numeric(names(mean_density)),
  interpretation = "Posterior lineage location density and separate ancestral environmental predictor only."
), file.path(dirs$maps, "case05_v6_map_values.rds"), compress = "gzip")

message("Case05 v6 complete: ", output)
message("No colonisation, persistence, occupancy, extinction, richness, Shannon, Simpson, PD, or FD maps were produced.")
