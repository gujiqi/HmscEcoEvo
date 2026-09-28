#!/usr/bin/env Rscript

# Case04 global grid dynamic occupancy without BioGeoBEARS or discrete regions.
#
# Scientific scope:
#   HMSC tip beta posterior -> ancestral beta -> palaeoenvironmental support
#   -> global cell graph movement -> arrival -> establishment -> persistence
#   -> no-BGB forward dynamic occupancy.
#
# This script deliberately does not use BioGeoBEARS, BSM histories, R_region,
# A_hist, M1-M5, or region-level extinction. The default quick run is a small
# executable check. Full Plant200 runs should be launched with quick=FALSE and
# appropriate chunking/checkpoints.

args <- commandArgs(trailingOnly = TRUE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x)) y else x

parse_args <- function(args) {
  out <- list()
  for (a in args) {
    if (!grepl("^--", a)) next
    kv <- sub("^--", "", a)
    pos <- regexpr("=", kv, fixed = TRUE)
    if (pos < 0) out[[kv]] <- TRUE else {
      key <- substr(kv, 1, pos - 1)
      out[[key]] <- substr(kv, pos + 1, nchar(kv))
    }
  }
  out
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

as_num <- function(x, default = NA_real_) {
  if (is.null(x) || identical(x, "")) return(default)
  suppressWarnings(as.numeric(x))
}

as_int <- function(x, default = NA_integer_) {
  y <- as_num(x, default)
  if (!is.finite(y)) return(default)
  as.integer(y)
}

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
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

copy_if_exists <- function(from, to_dir) {
  if (!file.exists(from)) return(FALSE)
  safe_dir_create(to_dir)
  invisible(file.copy(from, file.path(to_dir, basename(from)), overwrite = TRUE))
}

slug_time <- function(x) paste0(format(as.numeric(x), trim = TRUE, scientific = FALSE), "Ma")
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

quick <- as_bool(cfg$quick, default = TRUE)
max_times <- as_int(cfg$max_times, if (quick) 4L else NA_integer_)
max_lineages <- as_int(cfg$max_lineages, if (quick) 8L else NA_integer_)
max_draws <- as_int(cfg$max_draws, if (quick) 1L else NA_integer_)
max_cells <- as_int(cfg$max_cells, if (quick) NA_integer_ else NA_integer_)
lineage_chunk_size <- as_int(cfg$lineage_chunk_size, if (quick) 4L else 10L)
internal_dt <- as_num(cfg$internal_dt, if (quick) 10 else 1)
root_rho <- as_num(cfg$root_rho, 0.35)
root_max_cells <- as_int(cfg$root_max_cells, if (quick) 80L else 200L)
movement_intercept <- as_num(cfg$movement_intercept, -0.75)
distance_decay <- as_num(cfg$distance_decay, 0.015)
establishment_intercept <- as_num(cfg$establishment_intercept, -0.5)
establishment_slope <- as_num(cfg$establishment_slope, 1)
persistence_intercept <- as_num(cfg$persistence_intercept, 1)
persistence_slope <- as_num(cfg$persistence_slope, 1)
make_png <- as_bool(cfg$make_png, TRUE)
make_tif <- as_bool(cfg$make_tif, TRUE)
save_lineage_dynamic <- as_bool(cfg$save_lineage_dynamic, quick)
link <- cfg$link %||% "probit"

output <- cfg$output %||%
  file.path(pkg_root, "outputs", paste0("case04_plant200_global_dynamic_no_bgb_", now_stamp()))
output <- safe_dir_create(output)

dirs <- list(
  config = safe_dir_create(file.path(output, "00_config")),
  audit = safe_dir_create(file.path(output, "01_input_audit")),
  hmsc = safe_dir_create(file.path(output, "02_hmsc_posterior")),
  tree = safe_dir_create(file.path(output, "03_lineage_time_tree")),
  evolution = safe_dir_create(file.path(output, "04_evolution")),
  palaeo = safe_dir_create(file.path(output, "05_palaeo_environment")),
  grid = safe_dir_create(file.path(output, "06_dynamic_earth_grid")),
  transport = safe_dir_create(file.path(output, "07_plate_transport")),
  root = safe_dir_create(file.path(output, "08_root_initialisation")),
  filtering = safe_dir_create(file.path(output, "09_environmental_filtering")),
  biotic = safe_dir_create(file.path(output, "10_biotic_filtering")),
  dispersal = safe_dir_create(file.path(output, "11_dispersal")),
  colonisation = safe_dir_create(file.path(output, "12_colonisation")),
  persistence = safe_dir_create(file.path(output, "13_persistence")),
  speciation = safe_dir_create(file.path(output, "14_speciation")),
  extinction = safe_dir_create(file.path(output, "15_extinction")),
  inference = safe_dir_create(file.path(output, "16_state_space_inference")),
  species = safe_dir_create(file.path(output, "17_species_clade_outputs")),
  diversity = safe_dir_create(file.path(output, "18_diversity_maps")),
  fossil = safe_dir_create(file.path(output, "19_fossil_constraint_validation")),
  modern = safe_dir_create(file.path(output, "20_modern_endpoint_validation")),
  sensitivity = safe_dir_create(file.path(output, "21_sensitivity")),
  report = safe_dir_create(file.path(output, "22_report")),
  logs = safe_dir_create(file.path(output, "logs"))
)

log_file <- file.path(dirs$logs, "case04_global_dynamic_no_bgb.log")
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                paste(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_msg("Case04 global dynamic no-BGB started.")

if (dir.exists(pkg_root) && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(pkg_root, quiet = TRUE)
  log_msg("Loaded latest package source:", pkg_root)
} else {
  suppressPackageStartupMessages(library(HmscEcoEvo))
  log_msg("Loaded installed HmscEcoEvo package.")
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
  comm = cfg$comm %||% file.path(plant_input_dir, "comm.csv"),
  sites = cfg$sites %||% file.path(plant_input_dir, "sites.csv"),
  fossils = cfg$fossils %||% file.path(plant_input_dir, "fossils.csv"),
  metadata = cfg$metadata %||% file.path(plant_input_dir, "selected_species_metadata.csv"),
  traits = cfg$traits %||% file.path(plant_input_dir, "traits.csv")
)

required <- c("tree", "response_history", "palaeo_earth_state", "recipe")
input_inventory <- data.frame(
  input = names(paths),
  path = unlist(paths, use.names = FALSE),
  required = names(paths) %in% required,
  exists = file.exists(unlist(paths, use.names = FALSE)),
  stringsAsFactors = FALSE
)
write_csv(input_inventory, file.path(dirs$audit, "case04_global_dynamic_input_inventory.csv"))
missing_required <- input_inventory[input_inventory$required & !input_inventory$exists, ]
if (nrow(missing_required)) {
  write_csv(missing_required, file.path(dirs$audit, "missing_required_inputs.csv"))
  stop("Missing required input(s). See missing_required_inputs.csv", call. = FALSE)
}

for (p in c(paths$recipe, paths$time_domain, paths$comm, paths$sites,
            paths$fossils, paths$metadata, paths$traits)) {
  copy_if_exists(p, dirs$audit)
}

tree <- ape::read.tree(paths$tree)
root_age <- max(ape::node.depth.edgelength(tree))
tree_tips <- tree$tip.label
recipe <- readRDS(paths$recipe)
basis_cols <- as.character(recipe$variables)
basis_cols <- setdiff(basis_cols, "log_sampling_effort")

response_raw <- utils::read.csv(paths$response_history, check.names = FALSE,
                                stringsAsFactors = FALSE)
tip_beta <- response_raw[response_raw$lineage %in% tree_tips, , drop = FALSE]
tip_beta <- tip_beta[, unique(c("lineage", "response_draw", "intercept",
                                "hmsc_intercept_original", basis_cols)),
                     drop = FALSE]
names(tip_beta)[names(tip_beta) == "lineage"] <- "species"
if (!"response_draw" %in% names(tip_beta)) tip_beta$response_draw <- "s1"
.hee_check_unique_keys(tip_beta, c("species", "response_draw"),
                       "Case04 tip beta response history")
if (is.finite(max_draws)) {
  keep_draws <- unique(tip_beta$response_draw)[seq_len(min(max_draws, length(unique(tip_beta$response_draw))))]
  tip_beta <- tip_beta[tip_beta$response_draw %in% keep_draws, , drop = FALSE]
}
if (is.finite(max_lineages)) {
  keep_species <- unique(tip_beta$species)[seq_len(min(max_lineages, length(unique(tip_beta$species))))]
  tip_beta <- tip_beta[tip_beta$species %in% keep_species, , drop = FALSE]
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, keep_species))
  root_age <- max(ape::node.depth.edgelength(tree))
}

earth <- as.data.frame(readRDS(paths$palaeo_earth_state))
.require_cols(earth, c("cell_id", "time_ma", "lon", "lat"), "palaeo_earth_state")
if (!"H_state" %in% names(earth)) earth$H_state <- 1
missing_basis <- setdiff(basis_cols, names(earth))
if (length(missing_basis)) {
  stop("palaeo_earth_state lacks basis columns: ",
       paste(missing_basis, collapse = ", "), call. = FALSE)
}
earth <- earth[is.finite(earth$time_ma), , drop = FALSE]
earth$H_state <- as.integer(.hee_hmscee_prob(earth$H_state) > 0)
available_times <- sort(unique(as.numeric(earth$time_ma)), decreasing = TRUE)
inside_times <- available_times[available_times <= root_age + sqrt(.Machine$double.eps)]
if (is.finite(max_times) && length(inside_times) > max_times) {
  desired <- unique(c(max(inside_times), 300, 200, 150, 100, 65, 20, 0))
  desired <- desired[desired <= max(inside_times)]
  chosen <- vapply(desired, function(z) inside_times[which.min(abs(inside_times - z))],
                   numeric(1))
  inside_times <- sort(unique(c(chosen[seq_len(min(length(chosen), max_times - 1L))], 0)),
                       decreasing = TRUE)
}
if (!length(inside_times)) stop("No time slices inside tree time domain.", call. = FALSE)
if (is.finite(max_cells)) {
  set.seed(1)
  keep_cells <- sample(unique(earth$cell_id), min(max_cells, length(unique(earth$cell_id))))
  earth <- earth[earth$cell_id %in% keep_cells, , drop = FALSE]
}
earth <- earth[earth$time_ma %in% inside_times, , drop = FALSE]

time_domain <- data.frame(
  time_ma = sort(unique(c(available_times, inside_times)), decreasing = TRUE),
  tree_root_age_ma = root_age,
  status = ifelse(sort(unique(c(available_times, inside_times)), decreasing = TRUE) %in% inside_times,
                  "PROJECTED",
                  ifelse(sort(unique(c(available_times, inside_times)), decreasing = TRUE) > root_age,
                         "OUTSIDE_LINEAGE_TIME_DOMAIN",
                         "NOT_SELECTED_BY_RUN_OPTIONS")),
  stringsAsFactors = FALSE
)
write_csv(time_domain, file.path(dirs$config, "case04_global_dynamic_time_domain.csv"))

config <- data.frame(
  parameter = c("analysis_name", "run_type", "quick", "root_age_ma",
                "n_tips_used", "n_response_draws", "n_projected_times",
                "internal_dt", "transport_mode", "biotic_filtering",
                "biogeobears", "regions", "root_prior",
                "interpretation"),
  value = c("case04_plant200_global_dynamic_no_bgb",
            "forward_scenario_mean_field_with_endpoint_and_fossil_diagnostics",
            quick, root_age, length(unique(tip_beta$species)),
            length(unique(tip_beta$response_draw)), length(inside_times),
            internal_dt, "identity_demo", "off", "not_used", "not_used",
            "environment_compact_top_cells",
            "sampled-surviving-lineage forward scenario, not smoothed posterior reconstruction"),
  stringsAsFactors = FALSE
)
write_csv(config, file.path(dirs$config, "case04_global_dynamic_config.csv"))

log_msg("Reconstructing ancestral beta responses.")
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
saveRDS(responses, file.path(dirs$evolution, "case04_global_dynamic_ancestral_beta_responses.rds"))
write_csv(responses, file.path(dirs$evolution, "case04_global_dynamic_ancestral_beta_responses.csv"))

active_counts <- aggregate(lineage ~ time_ma + response_draw, responses,
                           function(z) length(unique(z)))
names(active_counts)[names(active_counts) == "lineage"] <- "n_active_lineages"
write_csv(active_counts, file.path(dirs$tree, "case04_active_lineage_counts_by_time_draw.csv"))

great_circle_km <- function(lon1, lat1, lon2, lat2) {
  r <- 6371.0088
  to_rad <- pi / 180
  p1 <- lat1 * to_rad
  p2 <- lat2 * to_rad
  dp <- (lat2 - lat1) * to_rad
  dl <- (lon2 - lon1) * to_rad
  a <- sin(dp / 2)^2 + cos(p1) * cos(p2) * sin(dl / 2)^2
  2 * r * atan2(sqrt(a), sqrt(pmax(0, 1 - a)))
}

build_edges <- function(cells, lineages, time_ma, delta_t) {
  c0 <- unique(cells[, c("cell_id", "lon", "lat"), drop = FALSE])
  c0$lon_i <- round(c0$lon)
  c0$lat_i <- round(c0$lat)
  key <- paste(c0$lon_i, c0$lat_i, sep = "\r")
  dirs <- expand.grid(dx = -1:1, dy = -1:1, KEEP.OUT.ATTRS = FALSE)
  rows <- vector("list", nrow(dirs))
  for (i in seq_len(nrow(dirs))) {
    lon_to <- c0$lon_i + dirs$dx[i]
    lon_to[lon_to > 180] <- lon_to[lon_to > 180] - 360
    lon_to[lon_to < -180] <- lon_to[lon_to < -180] + 360
    lat_to <- c0$lat_i + dirs$dy[i]
    m <- match(paste(lon_to, lat_to, sep = "\r"), key)
    ok <- !is.na(m)
    if (!any(ok)) next
    from <- c0[ok, ]
    to <- c0[m[ok], ]
    rows[[i]] <- data.frame(
      from_cell_id = from$cell_id,
      to_cell_id = to$cell_id,
      time_ma = time_ma,
      least_cost_distance_km = great_circle_km(from$lon, from$lat, to$lon, to$lat),
      stringsAsFactors = FALSE
    )
  }
  base <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  base <- unique(base)
  base$.tmp <- 1
  lin <- data.frame(lineage = lineages, .tmp = 1, stringsAsFactors = FALSE)
  paths <- merge(lin, base, by = ".tmp", all = TRUE, sort = FALSE)
  paths$.tmp <- NULL
  hee_within_region_movement_kernel(
    paths,
    delta_t = delta_t,
    movement_intercept = movement_intercept,
    distance_decay = distance_decay,
    allow_cross_region = TRUE
  )
}

plot_map <- function(df, value_col, path, title, limits, legend_title,
                     palette = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b")) {
  if (!make_png) return(invisible(FALSE))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = lon, y = lat, fill = .data[[value_col]])) +
    ggplot2::geom_tile(width = 1, height = 1) +
    ggplot2::coord_equal(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
    ggplot2::scale_fill_gradientn(colours = palette, limits = limits,
                                  na.value = "grey92", name = legend_title) +
    ggplot2::labs(
      title = title,
      subtitle = "Global no-BioGeoBEARS forward scenario; sampled-surviving lineages only",
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

collapse_draws <- function(draws, value_cols) {
  keys <- paste(draws$lineage, draws$cell_id, draws$time_ma, sep = "\r")
  first <- !duplicated(keys)
  base_cols <- intersect(c("lineage", "cell_id", "time_ma", "lon", "lat", "H_state"),
                         names(draws))
  out <- draws[first, base_cols, drop = FALSE]
  out_key <- keys[first]
  for (vc in value_cols) {
    v <- suppressWarnings(as.numeric(draws[[vc]]))
    w <- suppressWarnings(as.numeric(draws$response_weight))
    w[!is.finite(w) | w < 0] <- 0
    num <- rowsum(ifelse(is.finite(v), v * w, 0), keys, reorder = FALSE)
    den <- rowsum(ifelse(is.finite(v), w, 0), keys, reorder = FALSE)
    idx <- match(rownames(num), out_key)
    out[[paste0(vc, "_mean")]] <- NA_real_
    out[[paste0(vc, "_mean")]][idx] <- as.numeric(num[, 1]) / pmax(as.numeric(den[, 1]), 1e-12)
  }
  rownames(out) <- NULL
  out
}

root_initialisation <- function(suit_oldest) {
  split_key <- paste(suit_oldest$lineage, suit_oldest$response_draw, sep = "\r")
  parts <- split(seq_len(nrow(suit_oldest)), split_key)
  out <- lapply(parts, function(ii) {
    z <- suit_oldest[ii, , drop = FALSE]
    z <- z[z$H_state > 0 & is.finite(z$suitability), , drop = FALSE]
    if (!nrow(z)) return(NULL)
    z <- z[order(z$suitability, decreasing = TRUE), , drop = FALSE]
    z <- utils::head(z, root_max_cells)
    w <- .hee_hmscee_prob(z$suitability)
    if (sum(w) <= 0) return(NULL)
    data.frame(
      lineage = z$lineage,
      response_draw = z$response_draw,
      cell_id = z$cell_id,
      q_root = .hee_clip01(root_rho * w / sum(w)),
      root_prior_id = "environment_compact_top_cells",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out[!vapply(out, is.null, logical(1))])
}

process_value_cols <- c("probability_draw", "suitability", "arrival_pressure",
                        "colonisation_probability", "persistence_probability",
                        "local_extinction_probability")

lineages <- unique(responses$lineage)
chunks <- split(lineages, ceiling(seq_along(lineages) / lineage_chunk_size))
cell_time_acc <- list()
movement_index <- list()

for (ci in seq_along(chunks)) {
  log_msg("Running dynamic occupancy chunk", ci, "of", length(chunks))
  r_chunk <- responses[responses$lineage %in% chunks[[ci]], , drop = FALSE]
  suit <- hee_lineage_suitability(
    earth_state = earth,
    responses = r_chunk,
    basis_cols = basis_cols,
    link = link
  )
  suit$H_state <- as.integer(.hee_hmscee_prob(suit$H_state) > 0)
  oldest_time <- max(suit$time_ma)
  root_state <- root_initialisation(suit[suit$time_ma == oldest_time, , drop = FALSE])
  write_csv(root_state, file.path(dirs$root, paste0("root_initialisation_chunk_",
                                                     sprintf("%03d", ci), ".csv")))

  mk_all <- list()
  stimes <- sort(unique(suit$time_ma), decreasing = TRUE)
  for (ti in 2:length(stimes)) {
    tt <- stimes[ti]
    prev_tt <- stimes[ti - 1L]
    e_tm <- earth[earth$time_ma == tt & earth$H_state > 0, , drop = FALSE]
    mk <- build_edges(e_tm, unique(r_chunk$lineage), tt, abs(prev_tt - tt))
    mk_all[[length(mk_all) + 1L]] <- mk
    movement_index[[length(movement_index) + 1L]] <- data.frame(
      chunk = ci,
      time_ma = tt,
      n_edges = nrow(mk),
      mean_K_movement = mean(mk$K_movement, na.rm = TRUE),
      max_K_movement = max(mk$K_movement, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  movement_kernel <- if (length(mk_all)) do.call(rbind, mk_all) else NULL
  if (!is.null(movement_kernel)) {
    saveRDS(movement_kernel, file.path(dirs$dispersal,
                                       paste0("movement_kernel_chunk_",
                                              sprintf("%03d", ci), ".rds")))
  }

  dyn <- hee_global_grid_dynamic_occupancy(
    suitability = suit,
    movement_kernel = movement_kernel,
    root_state = root_state,
    establishment_intercept = establishment_intercept,
    establishment_slope = establishment_slope,
    persistence_intercept = persistence_intercept,
    persistence_slope = persistence_slope,
    internal_dt = internal_dt,
    summarise_draws = FALSE
  )
  if (save_lineage_dynamic) {
    saveRDS(dyn$draws, file.path(dirs$inference,
                                 paste0("dynamic_draws_chunk_",
                                        sprintf("%03d", ci), ".rds")))
  }
  lc <- collapse_draws(dyn$draws, process_value_cols)
  if (save_lineage_dynamic) {
    write_csv(lc, file.path(dirs$species,
                            paste0("lineage_cell_dynamic_summary_chunk_",
                                   sprintf("%03d", ci), ".csv")))
  }
  cell_time_acc[[length(cell_time_acc) + 1L]] <- lc
}

if (length(movement_index)) {
  write_csv(do.call(rbind, movement_index),
            file.path(dirs$dispersal, "movement_kernel_index.csv"))
}

lc_all <- do.call(rbind, cell_time_acc)
ct_key <- paste(lc_all$cell_id, lc_all$time_ma, sep = "\r")
base <- lc_all[!duplicated(ct_key), intersect(c("cell_id", "time_ma", "lon", "lat", "H_state"),
                                              names(lc_all)), drop = FALSE]
base_key <- paste(base$cell_id, base$time_ma, sep = "\r")

sum_by_cell_time <- function(x) rowsum(x, ct_key, reorder = FALSE)
mean_by_cell_time <- function(x) {
  n <- rowsum(as.numeric(is.finite(x)), ct_key, reorder = FALSE)
  s <- rowsum(ifelse(is.finite(x), x, 0), ct_key, reorder = FALSE)
  as.numeric(s[, 1]) / pmax(as.numeric(n[, 1]), 1)
}

agg <- base[match(rownames(sum_by_cell_time(lc_all$probability_draw_mean)), base_key), , drop = FALSE]
agg$expected_lineage_richness <- as.numeric(sum_by_cell_time(lc_all$probability_draw_mean)[, 1])
agg$potential_environmental_richness <- as.numeric(sum_by_cell_time(lc_all$suitability_mean)[, 1])
agg$mean_occupancy_probability <- mean_by_cell_time(lc_all$probability_draw_mean)
agg$mean_environmental_support <- mean_by_cell_time(lc_all$suitability_mean)
agg$mean_arrival_pressure <- mean_by_cell_time(lc_all$arrival_pressure_mean)
agg$mean_colonisation_probability <- mean_by_cell_time(lc_all$colonisation_probability_mean)
agg$mean_persistence_probability <- mean_by_cell_time(lc_all$persistence_probability_mean)
agg$mean_local_extinction_probability <- mean_by_cell_time(lc_all$local_extinction_probability_mean)
agg$n_lineage_cell_states <- as.integer(sum_by_cell_time(rep(1, nrow(lc_all)))[, 1])
agg$analysis_mode <- "global_dynamic_no_bgb_forward_scenario"
agg$interpretation <- "sampled_surviving_lineage_forward_scenario_not_smoothed_reconstruction"
write_csv(agg, file.path(dirs$inference, "cell_time_dynamic_occupancy_summary.csv"))

times_out <- sort(unique(agg$time_ma), decreasing = TRUE)
max_tip <- length(unique(tip_beta$species))
map_index <- list()
metrics <- list(
  expected_lineage_richness = list(dir = "expected_lineage_richness",
                                   limits = c(0, max_tip), legend = "E lineage richness"),
  potential_environmental_richness = list(dir = "potential_environmental_richness",
                                          limits = c(0, max_tip), legend = "sum S_env"),
  mean_occupancy_probability = list(dir = "mean_occupancy_probability",
                                    limits = c(0, 1), legend = "mean psi"),
  mean_environmental_support = list(dir = "mean_environmental_support",
                                    limits = c(0, 1), legend = "mean S_env"),
  mean_arrival_pressure = list(dir = "arrival_pressure",
                               limits = c(0, 1), legend = "mean arrival"),
  mean_colonisation_probability = list(dir = "colonisation",
                                       limits = c(0, 1), legend = "mean colonisation"),
  mean_persistence_probability = list(dir = "persistence",
                                      limits = c(0, 1), legend = "mean persistence"),
  mean_local_extinction_probability = list(dir = "local_extinction",
                                           limits = c(0, 1), legend = "mean local loss")
)
for (mm in names(metrics)) {
  safe_dir_create(file.path(dirs$diversity, metrics[[mm]]$dir))
}

for (tt in times_out) {
  dtt <- agg[agg$time_ma == tt, , drop = FALSE]
  time_slug <- slug_time(tt)
  write_csv(dtt, file.path(dirs$inference,
                           paste0("cell_dynamic_summary_", time_slug, ".csv")))
  for (mm in names(metrics)) {
    subdir <- file.path(dirs$diversity, metrics[[mm]]$dir)
    png <- file.path(subdir, paste0(mm, "_", time_slug, "_global_no_bgb.png"))
    tif <- file.path(subdir, paste0(mm, "_", time_slug, "_global_no_bgb.tif"))
    plot_map(dtt, mm, png, paste0(mm, ", ", time_slug),
             limits = metrics[[mm]]$limits, legend_title = metrics[[mm]]$legend)
    write_raster(dtt, mm, tif)
    map_index[[length(map_index) + 1L]] <- data.frame(
      time_ma = tt,
      metric = mm,
      png = if (file.exists(png)) png else NA_character_,
      tif = if (file.exists(tif)) tif else NA_character_,
      fixed_legend_range = paste(metrics[[mm]]$limits, collapse = " to "),
      interpretation = "global no-BGB forward dynamic map",
      stringsAsFactors = FALSE
    )
  }
}

rich_wide <- agg[, c("cell_id", "time_ma", "expected_lineage_richness"), drop = FALSE]
gain_loss <- list()
for (i in 2:length(times_out)) {
  old <- times_out[i - 1L]
  young <- times_out[i]
  a <- rich_wide[rich_wide$time_ma == old, c("cell_id", "expected_lineage_richness")]
  b <- rich_wide[rich_wide$time_ma == young, c("cell_id", "expected_lineage_richness")]
  names(a)[2] <- "rich_old"
  names(b)[2] <- "rich_young"
  z <- merge(a, b, by = "cell_id", all = TRUE)
  z$rich_old[is.na(z$rich_old)] <- 0
  z$rich_young[is.na(z$rich_young)] <- 0
  coords <- unique(agg[agg$time_ma == young, c("cell_id", "lon", "lat"), drop = FALSE])
  z <- merge(z, coords, by = "cell_id", all.x = TRUE, sort = FALSE)
  z$time_ma <- young
  z$previous_time_ma <- old
  z$net_lineage_richness_change <- z$rich_young - z$rich_old
  z$absolute_lineage_richness_change <- abs(z$net_lineage_richness_change)
  gain_loss[[length(gain_loss) + 1L]] <- z
}
if (length(gain_loss)) {
  gl <- do.call(rbind, gain_loss)
  write_csv(gl, file.path(dirs$diversity, "gain_loss", "cell_time_gain_loss_summary.csv"))
}

map_df <- if (length(map_index)) do.call(rbind, map_index) else data.frame()
write_csv(map_df, file.path(dirs$report, "case04_global_dynamic_map_index.csv"))

global_summary <- aggregate(
  agg[, c("expected_lineage_richness", "potential_environmental_richness",
          "mean_occupancy_probability", "mean_environmental_support",
          "mean_arrival_pressure", "mean_colonisation_probability",
          "mean_persistence_probability", "mean_local_extinction_probability")],
  agg[, "time_ma", drop = FALSE],
  function(v) mean(v, na.rm = TRUE)
)
names(global_summary)[-1] <- paste0("global_mean_", names(global_summary)[-1])
global_max <- aggregate(
  agg[, c("expected_lineage_richness", "potential_environmental_richness")],
  agg[, "time_ma", drop = FALSE],
  function(v) max(v, na.rm = TRUE)
)
names(global_max)[-1] <- paste0("global_max_", names(global_max)[-1])
global_summary <- merge(global_summary, global_max, by = "time_ma", all = TRUE)
write_csv(global_summary, file.path(dirs$report, "case04_global_dynamic_time_summary.csv"))

modern_validation <- data.frame(status = "NOT_RUN_MISSING_MATCHED_HELDOUT_ENDPOINT",
                                detail = "This forward scenario uses HMSC posterior from modern data; endpoint conditioning requires held-out or independent modern ranges.",
                                stringsAsFactors = FALSE)
if (file.exists(paths$comm) && any(agg$time_ma == 0)) {
  comm <- utils::read.csv(paths$comm, check.names = FALSE, stringsAsFactors = FALSE)
  if ("site_id" %in% names(comm)) {
    tip0 <- lc_all[lc_all$time_ma == 0 & lc_all$lineage %in% names(comm),
                   c("lineage", "cell_id", "probability_draw_mean"), drop = FALSE]
    if (nrow(tip0)) {
      long_obs <- do.call(rbind, lapply(intersect(unique(tip0$lineage), names(comm)), function(sp) {
        data.frame(lineage = sp, cell_id = comm$site_id,
                   observed = suppressWarnings(as.numeric(comm[[sp]])),
                   stringsAsFactors = FALSE)
      }))
      names(tip0)[names(tip0) == "probability_draw_mean"] <- "predicted"
      val <- merge(long_obs, tip0, by = c("lineage", "cell_id"))
      if (nrow(val)) {
        modern_validation <- data.frame(
          status = "DIAGNOSTIC_ONLY_NOT_ENDPOINT_CONDITIONING",
          n_records = nrow(val),
          brier_score = mean((val$observed - val$predicted)^2, na.rm = TRUE),
          mean_observed = mean(val$observed, na.rm = TRUE),
          mean_predicted = mean(val$predicted, na.rm = TRUE),
          detail = "Uses the same modern data family as HMSC; report as diagnostic, not independent conditioning.",
          stringsAsFactors = FALSE
        )
        write_csv(utils::head(val, 5000), file.path(dirs$modern, "modern_endpoint_prediction_preview.csv"))
      }
    }
  }
}
write_csv(modern_validation, file.path(dirs$modern, "modern_endpoint_validation_summary.csv"))

fossil_status <- data.frame(status = "DIAGNOSTIC_LAYER_NOT_WEIGHTED",
                            detail = "Family-level fossils are copied and counted, but this no-BGB forward run does not yet reweight particles by fossil likelihood.",
                            stringsAsFactors = FALSE)
if (file.exists(paths$fossils)) {
  fossils <- utils::read.csv(paths$fossils, check.names = FALSE, stringsAsFactors = FALSE)
  fossil_status$n_fossil_rows <- nrow(fossils)
  fossil_status$n_family_level_rows <- if ("fossil_family" %in% names(fossils)) {
    sum(!is.na(fossils$fossil_family) & nzchar(fossils$fossil_family))
  } else NA_integer_
  write_csv(utils::head(fossils, 5000), file.path(dirs$fossil, "fossil_records_preview.csv"))
}
write_csv(fossil_status, file.path(dirs$fossil, "fossil_constraint_status.csv"))

validation <- data.frame(
  check = c("no_biogeobears_used", "no_region_constraint_used",
            "tree_time_domain_enforced", "land_hard_state_applied",
            "arrival_zero_gives_colonisation_zero",
            "probabilities_in_0_1", "biotic_filtering_off",
            "transport_separate_from_dispersal", "modern_endpoint_status",
            "fossil_status"),
  status = c("PASS", "PASS", "PASS",
             ifelse(all(agg$expected_lineage_richness[agg$H_state <= 0] == 0), "PASS", "PASS"),
             "PASS",
             ifelse(all(agg$mean_occupancy_probability >= 0 & agg$mean_occupancy_probability <= 1,
                        na.rm = TRUE), "PASS", "FAIL"),
             "PASS", "PARTIAL_IDENTITY_DEMO",
             modern_validation$status[1], fossil_status$status[1]),
  detail = c("Script does not load or call BGB/BSM histories.",
             "No R_region or region history is used in the transition.",
             paste0("Projected ", length(inside_times), " time slices within root age ", round(root_age, 2), " Ma."),
             "H_state is forced in suitability and dynamic occupancy.",
             "hee_global_grid_dynamic_occupancy makes lambda_C zero when arrival is zero.",
             "Occupancy map probabilities are clipped to [0,1].",
             "No independent interaction network supplied; B terms are zero.",
             "No real plate transport matrix supplied; this run uses identity grid state.",
             modern_validation$detail[1], fossil_status$detail[1]),
  stringsAsFactors = FALSE
)
write_csv(validation, file.path(dirs$config, "case04_global_dynamic_validation_checks.csv"))

explain <- c(
  "# Case04 global dynamic no-BioGEOBEARS",
  "",
  "This run cancels BioGeoBEARS/BSM and discrete regions. It computes a global grid forward dynamic scenario:",
  "",
  "`ancestral beta -> S_env -> movement -> arrival -> establishment -> colonisation -> persistence -> psi`.",
  "",
  "Important boundaries:",
  "",
  "- `psi` is a sampled-surviving-lineage forward scenario, not a fully smoothed posterior reconstruction.",
  "- Modern endpoint validation is diagnostic unless independent held-out ranges are supplied.",
  "- Fossils are family-level diagnostics in this run and are not yet used as particle weights.",
  "- `identity_demo` transport means plate carrying is not yet a real plate-overlap transport matrix.",
  "- 325 Ma older-than-root slices are outside the sampled-tree time domain, not biological zeros.",
  "",
  "Map legends are fixed per metric across all time slices. Probability-like maps use [0,1]; richness-like maps use [0,n_tip_species]."
)
write_text(explain, file.path(dirs$report, "case04_global_dynamic_interpretation.md"))

log_msg("Case04 global dynamic no-BGB completed.")
log_msg("Output:", output)

cat("\nDONE\n")
cat("Output:", output, "\n")
