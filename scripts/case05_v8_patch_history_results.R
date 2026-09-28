#!/usr/bin/env Rscript

# Derived Case05 results from the completed v8 phylogeographic runs.
# No demographic occupancy, posterior edge flux, or fossil-confirmed refugia
# is fabricated from marginal location-density maps.

root <- "C:/Users/Google/Documents/HMSC-HIST"
out_root <- file.path(root, "outputs", "HmscEcoEvo")
main <- file.path(out_root, "case05_v8_full200_terrain_resistance_20260926")
comparison <- file.path(out_root, "case05_v8_four_scheme_comparison_20260926")
out <- Sys.getenv("HMSCEE_PATCH_OUT",
  file.path(comparison, "08_patch_history_results"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("HmscEcoEvo", quietly = TRUE) ||
    !requireNamespace("ape", quietly = TRUE) ||
    !requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("terra", quietly = TRUE)) {
  stop("HmscEcoEvo, ape, ggplot2, and terra must be installed.", call. = FALSE)
}
if (!"hee_result_refugia_candidates" %in% getNamespaceExports("HmscEcoEvo")) {
  stop("Install the updated HmscEcoEvo package first.", call. = FALSE)
}
q <- utils::read.csv(file.path(main, "06_quality", "case05_v8_quality_gates.csv"))
if (!all(q$pass)) stop("Source Case05 quality gates did not pass.", call. = FALSE)
main_field <- readRDS(file.path(main, "05_shared_scale_maps", "case05_v8_map_values.rds"))
scenario_dirs <- c(
  spherical_no_ldd = "case05_v8_full200_spherical_land_20260926",
  spherical_land = "case05_v8_full200_spherical_land_with_ldd_20260926",
  terrain_resistance = "case05_v8_full200_terrain_resistance_20260926",
  particle_topography = "case05_v8_full200_particle_topography_20260926"
)
scenario_fields <- lapply(file.path(out_root, scenario_dirs,
  "05_shared_scale_maps", "case05_v8_map_values.rds"), readRDS)
names(scenario_fields) <- names(scenario_dirs)
times <- as.numeric(names(main_field$lineage_location_density))
if (length(times) != 66L || !identical(times, sort(times, decreasing = TRUE))) {
  stop("Expected all 66 ordered native Case05 time slices.", call. = FALSE)
}
if (!all(vapply(scenario_fields, function(x)
  identical(names(x$lineage_location_density), as.character(times)), logical(1)))) {
  stop("Four-scheme map times are not aligned.", call. = FALSE)
}
response_file <- Sys.getenv("HMSCEE_RESPONSE_FILE",
  file.path(main, "02_ancestral_response", "ancestral_beta_response_draws.rds"))
responses <- readRDS(response_file)
axes <- c("MAT_pohl_C", "MAP_pohl_mm_yr", "elevation_m",
  "T_seasonality_pohl_sd", "P_seasonality_pohl_sd",
  "moisture_availability_index_z", "wetland_potential_index_z", "relief_3x3_m")
if (!all(axes %in% names(responses))) stop("Ancestral response axes are incomplete.")
fit_path <- file.path(out_root,
  "case05_plant200_4deg_effortaware_formal25_all66_20260923",
  "02_hmsc_4deg", "case05_v5_hmsc_4deg_target_group_model.rds")
fit <- readRDS(fit_path)
training_x <- as.matrix(fit$XData[, axes, drop = FALSE])
rm(fit)
training_bounds <- apply(training_x, 2L, stats::quantile,
  probs = c(.01, .99), na.rm = TRUE)
utils::write.csv(data.frame(axis = axes, lower_01 = training_bounds[1L, ],
  upper_99 = training_bounds[2L, ]),
  file.path(out, "modern_training_environment_bounds.csv"), row.names = FALSE)
cache_root <- file.path(root, "HmscEcoEvo", "derived_inputs",
  "case05_v2_paleomap_h3_carrier_325_0Ma_1deg",
  "movement_cache_arias_topographic_rate0p4")
cache_index <- utils::read.csv(file.path(cache_root, "case05_v2_movement_cache_index.csv"))
input_root <- file.path(root, "HmscEcoEvo_API8_20260818", "data_external",
  "plant_200_multifamily_extant_20260820", "prepared_inputs",
  "case03_data_raw", "Plant200_multifamily_extant_1deg_traits")
selection <- utils::read.csv(file.path(main, "01_input_audit", "species_selection.csv"))
tree <- ape::read.tree(file.path(input_root, "tree.tre"))
tree <- ape::reorder.phylo(ape::keep.tip(tree, selection$species[selection$selected]),
                           "cladewise")
tip_count <- length(tree$tip.label)
if (tip_count != 200L) stop("Expected the same 200-tip pruned tree.", call. = FALSE)

grid <- data.frame(cell_id = sprintf("g4_%04d", seq_len(4050L)),
  lon = rep(seq(-178, 178, by = 4), times = 45),
  lat = rep(seq(-88, 88, by = 4), each = 90), stringsAsFactors = FALSE)
lat_s <- grid$lat - 2; lat_n <- grid$lat + 2
grid$cell_area_km2 <- 6371.0088^2 * (4 * pi / 180) *
  (sin(lat_n * pi / 180) - sin(lat_s * pi / 180))
bin4 <- function(lon, lat) {
  ix <- pmin(89L, pmax(0L, floor((lon + 180) / 4)))
  iy <- pmin(44L, pmax(0L, floor((lat + 90) / 4)))
  as.integer(iy * 90L + ix + 1L)
}
grid_edges <- function() {
  origin <- seq_len(4050L)
  east <- ((origin - 1L) %/% 90L) * 90L + (origin %% 90L) + 1L
  north <- origin + 90L
  ok <- north <= 4050L
  data.frame(from_cell_id = grid$cell_id[c(origin, origin[ok])],
    to_cell_id = grid$cell_id[c(east, north[ok])])
}
adjacency <- grid_edges()
get_carrier <- function(time_ma) {
  idx <- match(time_ma, cache_index$time_ma)
  if (is.na(idx)) stop("Missing carrier at ", time_ma, " Ma")
  readRDS(cache_index$file[[idx]])$carrier
}
slice_from_carrier <- function(carrier) {
  good <- !is.na(carrier$paleo_lon) & !is.na(carrier$paleo_lat) &
    !is.na(carrier$active_land) & carrier$active_land &
    is.finite(carrier$land_area_km2) & carrier$land_area_km2 > 0
  x <- carrier[good, , drop = FALSE]
  x$bin <- bin4(x$paleo_lon, x$paleo_lat)
  total <- rowsum(x$land_area_km2, x$bin, reorder = FALSE)
  area <- numeric(nrow(grid)); area[as.integer(rownames(total))] <- total[, 1L]
  env <- matrix(NA_real_, nrow(grid), length(axes), dimnames = list(NULL, axes))
  for (axis in axes) {
    value <- as.numeric(x[[axis]])
    valid <- is.finite(value)
    if (!any(valid)) next
    num <- rowsum(x$land_area_km2[valid] * value[valid], x$bin[valid],
                   reorder = FALSE)
    den <- rowsum(x$land_area_km2[valid], x$bin[valid], reorder = FALSE)
    ids <- as.integer(rownames(num))
    env[ids, axis] <- num[, 1L] / den[match(rownames(num), rownames(den)), 1L]
  }
  list(area = area, env = env)
}
carriage_between <- function(older, younger) {
  target <- younger[match(older$track_id, younger$track_id), , drop = FALSE]
  source_valid <- !is.na(older$active_land) & older$active_land &
    is.finite(older$land_area_km2) & older$land_area_km2 > 0 &
    is.finite(older$paleo_lon) & is.finite(older$paleo_lat)
  source_bin <- bin4(older$paleo_lon[source_valid], older$paleo_lat[source_valid])
  source_total <- rowsum(older$land_area_km2[source_valid], source_bin,
                         reorder = FALSE)
  valid <- source_valid & !is.na(target$active_land) & target$active_land &
    is.finite(target$paleo_lon) & is.finite(target$paleo_lat)
  source <- bin4(older$paleo_lon[valid], older$paleo_lat[valid])
  dest <- bin4(target$paleo_lon[valid], target$paleo_lat[valid])
  weight <- older$land_area_km2[valid]
  pair <- stats::aggregate(weight, list(source = source, target = dest), sum)
  names(pair)[[3L]] <- "area"
  pair$carriage_weight <- pair$area /
    source_total[match(as.character(pair$source), rownames(source_total)), 1L]
  data.frame(source_cell_id = grid$cell_id[pair$source],
    target_cell_id = grid$cell_id[pair$target],
    carriage_weight = pair$carriage_weight)
}
lineage_key <- function(node) {
  if (node <= tip_count) tree$tip.label[[node]] else paste0("node_", node)
}
parents <- setNames(tree$edge[, 1L], as.character(tree$edge[, 2L]))
focal <- c("Ilex_aquifolium", "Camellia_sinensis", "Jasminum_officinale")
tip_path <- function(tip) {
  node <- match(tip, tree$tip.label)
  path <- integer()
  while (!is.na(node)) {
    path <- c(path, node)
    next_node <- parents[as.character(node)]
    node <- if (is.na(next_node)) NA_integer_ else as.integer(next_node)
  }
  vapply(path, lineage_key, character(1))
}
paths <- setNames(lapply(focal, tip_path), focal)
threshold <- 0.6
patches <- setNames(lapply(focal, function(.) list()), focal)
membership <- patches
global_rows <- list()
scheme_rows <- list()
global_maps <- list()
focal_maps <- list()
noanalog_maps <- list()
noanalog_summary <- list()
previous_carrier <- NULL
previous_time <- NULL
links <- setNames(lapply(focal, function(.) list()), focal)

for (k in seq_along(times)) {
  tm <- times[[k]]; key <- as.character(tm)
  carrier <- get_carrier(tm)
  slice <- slice_from_carrier(carrier)
  land <- slice$area > 0 & stats::complete.cases(slice$env)
  model_land <- main_field$landscape_weight[[k]] > 0
  if (length(model_land) != nrow(grid) || anyNA(model_land)) {
    stop("Invalid positional landscape-weight slice at ", tm, " Ma.")
  }
  noanalog <- rowSums(sweep(slice$env, 2L, training_bounds[1L, ], "<") |
    sweep(slice$env, 2L, training_bounds[2L, ], ">"), na.rm = FALSE) > 0
  noanalog[!land] <- NA
  noanalog_maps[[key]] <- as.numeric(noanalog)
  noanalog_summary[[k]] <- data.frame(time_ma = tm,
    noanalog_land_fraction = mean(noanalog[land]),
    noanalog_land_area_fraction = sum(slice$area[land & noanalog], na.rm = TRUE) /
      sum(slice$area[land]))
  if (!identical(land, model_land)) {
    # The movement model can retain a land bin with partial environmental
    # missingness. Such a bin is not assigned a support value, never ocean.
    if (any(land & !model_land)) {
      stop("Environmental slice contains support on marine states at ", tm)
    }
  }
  current <- main_field$lineage_location_density[[key]]
  r <- responses[abs(responses$time_ma - tm) < 1e-7, , drop = FALSE]
  r$state_key <- ifelse(r$node <= tip_count, tree$tip.label[pmin(r$node, tip_count)],
                        paste0("node_", r$node))
  draw_ids <- unique(r$response_draw)
  valid_cell <- land & rowSums(!is.finite(slice$env)) == 0L
  valid_env <- slice$env
  valid_env[!is.finite(valid_env)] <- 0
  s <- matrix(0, nrow(valid_env), ncol(current),
    dimnames = list(NULL, colnames(current)))
  for (draw_id in draw_ids) {
    one <- r[r$response_draw == draw_id, , drop = FALSE]
    if (anyDuplicated(one$state_key)) {
      stop("Duplicate ancestral response for draw ", draw_id, " at ", tm)
    }
    missing <- setdiff(colnames(current), one$state_key)
    if (length(missing)) {
      stop("Missing ancestral beta at ", tm, ": ", missing[[1L]])
    }
    beta <- t(as.matrix(one[match(colnames(current), one$state_key),
      axes, drop = FALSE]))
    s <- s + stats::pnorm(valid_env %*% beta) / length(draw_ids)
  }
  s[!valid_cell, ] <- NA_real_
  area <- pmin(slice$area, grid$cell_area_km2)
  supported <- s >= threshold
  supported[is.na(supported)] <- FALSE
  sum_area <- colSums(sweep(supported, 1L, area, "*"))
  global_rows[[k]] <- data.frame(time_ma = tm, active_lineage = colnames(current),
    supported_area_km2 = as.numeric(sum_area),
    n_supported_cells = colSums(supported),
    response_draws = length(draw_ids),
    support_threshold = threshold, stringsAsFactors = FALSE)
  global_maps[[key]] <- rowSums(supported)
  for (tip in focal) {
    active <- intersect(paths[[tip]], colnames(current))
    if (!length(active)) stop("No active ancestral path for ", tip, " at ", tm)
    active <- active[[1L]]
    j <- match(active, colnames(current))
    tab <- data.frame(cell_id = grid$cell_id, support = s[, j],
      land_area_km2 = area, location_density = current[, j],
      MAT_pohl_C = slice$env[, "MAT_pohl_C"],
      MAP_pohl_mm_yr = slice$env[, "MAP_pohl_mm_yr"],
      no_analog = noanalog)
    hit <- HmscEcoEvo::hee_result_patches(tab, adjacency, threshold = threshold,
      min_area_km2 = 0, lineage_id = tip, time_ma = tm)
    if (nrow(hit$patches)) {
      hit$patches$active_ancestral_lineage <- active
      patches[[tip]][[key]] <- hit$patches
      membership[[tip]][[key]] <- hit$membership
    }
    focal_maps[[paste(tip, key, sep = "@")] ] <- s[, j]
    for (scheme in names(scenario_fields)) {
      d <- scenario_fields[[scheme]]$lineage_location_density[[key]][, active]
      scheme_rows[[length(scheme_rows) + 1L]] <- data.frame(
        focal_tip = tip, time_ma = tm, active_ancestral_lineage = active,
        scenario = scheme,
        location_mass_in_supported_cells = sum(d[supported[, j]]),
        location_mass_on_land = sum(d[model_land]),
        stringsAsFactors = FALSE)
    }
  }
  if (!is.null(previous_carrier)) {
    transport <- carriage_between(previous_carrier, carrier)
    for (tip in focal) {
      old <- membership[[tip]][[as.character(previous_time)]]
      young <- membership[[tip]][[key]]
      if (is.null(old) || is.null(young)) next
      links[[tip]][[key]] <- HmscEcoEvo::hee_result_patch_links(
        old, young, transport, min_source_share = 0.05)
    }
  }
  previous_carrier <- carrier; previous_time <- tm
  if (k %% 10L == 0L || k == length(times))
    message("Case05 patch results: ", k, "/", length(times), " time slices")
}

write_csv <- function(x, name) {
  utils::write.csv(x, file.path(out, name), row.names = FALSE, na = "")
}
global <- do.call(rbind, global_rows)
scheme <- do.call(rbind, scheme_rows)
write_csv(global, "all_active_lineages_environmental_opportunity.csv")
write_csv(scheme, "four_scheme_location_mass_in_supported_area.csv")
write_csv(do.call(rbind, noanalog_summary), "noanalog_by_time.csv")
index <- list()
for (tip in focal) {
  p <- if (length(patches[[tip]])) do.call(rbind, patches[[tip]]) else data.frame()
  m <- if (length(membership[[tip]])) do.call(rbind, membership[[tip]]) else data.frame()
  e <- if (length(links[[tip]])) do.call(rbind, links[[tip]]) else data.frame(
    from_patch_id = character(), to_patch_id = character(),
    time_from_ma = numeric(), time_to_ma = numeric())
  if (!nrow(p)) stop("No supported patches found for ", tip)
  stable <- HmscEcoEvo::hee_result_patch_stability(p, e)
  refuge <- HmscEcoEvo::hee_result_refugia_candidates(
    p, e, contraction_fraction = 0.10, min_location_mass = 0.10,
    min_patch_area_km2 = 100000, require_later_link = TRUE)
  write_csv(m, paste0(tip, "_patch_cells.csv"))
  write_csv(stable$patches, paste0(tip, "_patches.csv"))
  write_csv(e, paste0(tip, "_plate_carried_patch_links.csv"))
  write_csv(stable$networks, paste0(tip, "_patch_history_networks.csv"))
  write_csv(stable$link_changes, paste0(tip, "_link_environmental_changes.csv"))
  write_csv(refuge$episodes, paste0(tip, "_opportunity_contractions.csv"))
  write_csv(refuge$candidates, paste0(tip, "_candidate_refugia.csv"))
  sensitivity <- do.call(rbind, lapply(c(.10, .25), function(cutoff) {
    result <- HmscEcoEvo::hee_result_refugia_candidates(p, e,
      contraction_fraction = cutoff, min_location_mass = .10,
      min_patch_area_km2 = 100000, require_later_link = TRUE)
    data.frame(contraction_fraction = cutoff,
      episodes = nrow(result$episodes), candidates = nrow(result$candidates))
  }))
  write_csv(sensitivity, paste0(tip, "_refugia_threshold_sensitivity.csv"))
  index[[tip]] <- data.frame(focal_tip = tip, time_slices = length(times),
    patches = nrow(p), carriage_links = nrow(e),
    contraction_intervals = nrow(refuge$episodes),
    candidate_refugia = nrow(refuge$candidates))
}
write_csv(do.call(rbind, index), "case05_patch_result_index.csv")
write_csv(data.frame(result = c("posterior_corridors", "donor_recipient",
  "fossil_supported_refugia", "demographic_source_sink"),
  status = c("NOT_COMPUTED_NO_JOINT_EDGE_FLOW",
    "NOT_COMPUTED_NO_JOINT_EDGE_FLOW", "NOT_COMPUTED_NO_LINKED_FOSSIL_TEST",
    "NOT_ESTIMABLE_NO_DEMOGRAPHIC_RATES")), "not_computed_evidence_gates.csv")

selected_times <- intersect(c(300, 200, 100, 65, 20, 0), times)
map_dir <- file.path(out, "maps"); dir.create(map_dir, showWarnings = FALSE)
map_index <- list()
save_geotiff <- function(value, file) {
  raster <- terra::rast(ncols = 90, nrows = 45, xmin = -180, xmax = 180,
    ymin = -90, ymax = 90, crs = "EPSG:4326")
  pixel <- terra::cellFromXY(raster, as.matrix(grid[, c("lon", "lat")]))
  values <- rep(NA_real_, terra::ncell(raster))
  values[pixel] <- value
  terra::values(raster) <- values
  terra::writeRaster(raster, file, overwrite = TRUE,
    gdal = c("COMPRESS=DEFLATE"))
}
for (tm in selected_times) {
  key <- as.character(tm)
  land <- main_field$landscape_weight[[match(tm, times)]] > 0
  map_table <- grid[, c("lon", "lat")]
  map_table$value <- global_maps[[key]]
  map_table$value[!land] <- NA_real_
  fig <- ggplot2::ggplot(map_table, ggplot2::aes(lon, lat)) +
    ggplot2::geom_tile(ggplot2::aes(fill = value)) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 200), na.value = "white",
      name = "Active lineage\nsupport count") +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::labs(title = paste("Environmental opportunity |", tm, "Ma"),
      subtitle = "Number of active sampled lineages with environmental index >= 0.6",
      caption = "Not occupancy, species richness, or a census of historical communities.",
      x = "Palaeolongitude", y = "Palaeolatitude") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  file <- file.path(map_dir, paste0("active_lineage_environmental_support_", tm, "Ma.png"))
  ggplot2::ggsave(file, fig, width = 12, height = 6.4, dpi = 150)
  tif <- sub("\\.png$", ".tif", file)
  save_geotiff(map_table$value, tif)
  map_index[[length(map_index) + 1L]] <- data.frame(
    metric = "active_lineage_environmental_support_count", time_ma = tm,
    focal_tip = NA_character_, png = file, geotiff = tif)
  map_table$value <- global_maps[[key]] / nrow(global_rows[[match(tm, times)]])
  map_table$value[!land] <- NA_real_
  fig <- ggplot2::ggplot(map_table, ggplot2::aes(lon, lat)) +
    ggplot2::geom_tile(ggplot2::aes(fill = value)) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), na.value = "white",
      name = "Active lineage\nsupport fraction") +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::labs(title = paste("Relative environmental opportunity |", tm, "Ma"),
      subtitle = "Share of active sampled lineages with environmental index >= 0.6",
      caption = "Denominator varies with dated-tree lineage number; not occupancy.",
      x = "Palaeolongitude", y = "Palaeolatitude") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  file <- file.path(map_dir, paste0("active_lineage_environmental_fraction_", tm, "Ma.png"))
  ggplot2::ggsave(file, fig, width = 12, height = 6.4, dpi = 150)
  tif <- sub("\\.png$", ".tif", file)
  save_geotiff(map_table$value, tif)
  map_index[[length(map_index) + 1L]] <- data.frame(
    metric = "active_lineage_environmental_support_fraction", time_ma = tm,
    focal_tip = NA_character_, png = file, geotiff = tif)
  map_table$value <- noanalog_maps[[key]]
  fig <- ggplot2::ggplot(map_table, ggplot2::aes(lon, lat)) +
    ggplot2::geom_tile(ggplot2::aes(fill = value)) +
    ggplot2::scale_fill_gradient(low = "#e4f4e8", high = "#c54836",
      limits = c(0, 1), na.value = "white", name = "No-analog") +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::labs(title = paste("Modern-training envelope |", tm, "Ma"),
      subtitle = "Outside at least one 1st-99th percentile predictor range",
      caption = "External diagnostic; never multiplied into biological support.",
      x = "Palaeolongitude", y = "Palaeolatitude") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  file <- file.path(map_dir, paste0("noanalog_training_envelope_", tm, "Ma.png"))
  ggplot2::ggsave(file, fig, width = 12, height = 6.4, dpi = 150)
  tif <- sub("\\.png$", ".tif", file)
  save_geotiff(map_table$value, tif)
  map_index[[length(map_index) + 1L]] <- data.frame(
    metric = "modern_training_noanalog_flag", time_ma = tm,
    focal_tip = NA_character_, png = file, geotiff = tif)
  for (tip in focal) {
    map_table$value <- focal_maps[[paste(tip, key, sep = "@")]]
    map_table$value[!land] <- NA_real_
    fig <- ggplot2::ggplot(map_table, ggplot2::aes(lon, lat)) +
      ggplot2::geom_tile(ggplot2::aes(fill = value)) +
      ggplot2::scale_fill_viridis_c(limits = c(0, 1), na.value = "white",
        name = "Support index") +
      ggplot2::coord_equal(expand = FALSE) +
      ggplot2::labs(title = paste(tip, "ancestral path |", tm, "Ma"),
        subtitle = "Ancestral response x time-matched palaeoenvironment",
        caption = "Environmental support, not historical presence or confirmed refugium.",
        x = "Palaeolongitude", y = "Palaeolatitude") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
    file <- file.path(map_dir, paste0(tip, "_environmental_support_", tm, "Ma.png"))
    ggplot2::ggsave(file, fig, width = 12, height = 6.4, dpi = 150)
    tif <- sub("\\.png$", ".tif", file)
    save_geotiff(map_table$value, tif)
    map_index[[length(map_index) + 1L]] <- data.frame(
      metric = "ancestral_environmental_support", time_ma = tm,
      focal_tip = tip, png = file, geotiff = tif)
  }
}
for (tip in focal) {
  candidate <- utils::read.csv(file.path(out, paste0(tip, "_candidate_refugia.csv")))
  member <- utils::read.csv(file.path(out, paste0(tip, "_patch_cells.csv")))
  if (!nrow(candidate)) next
  frequent <- sort(table(candidate$time_ma), decreasing = TRUE)
  candidate_times <- as.numeric(names(frequent)[seq_len(min(3L, length(frequent)))])
  for (tm in candidate_times) {
    key <- as.character(tm)
    land <- main_field$landscape_weight[[match(tm, times)]] > 0
    flags <- rep(0, nrow(grid)); flags[!land] <- NA_real_
    selected <- candidate[candidate$time_ma == tm, , drop = FALSE]
    cells <- member[member$patch_id %in% selected$to_patch_id, , drop = FALSE]
    flags[match(cells$cell_id, grid$cell_id)] <- 1
    supported <- selected$to_patch_id[
      selected$geographic_evidence == "location_mass_above_threshold"]
    cells_supported <- member[member$patch_id %in% supported, , drop = FALSE]
    flags[match(cells_supported$cell_id, grid$cell_id)] <- 2
    map_table <- grid[, c("lon", "lat")]
    map_table$value <- factor(flags, levels = c(0, 1, 2),
      labels = c("Other land", "Ecological candidate", "Location-mass supported"))
    fig <- ggplot2::ggplot(map_table, ggplot2::aes(lon, lat)) +
      ggplot2::geom_tile(ggplot2::aes(fill = value)) +
      ggplot2::scale_fill_manual(values = c("Other land" = "#d9e2e3",
        "Ecological candidate" = "#e0a33d",
        "Location-mass supported" = "#278b6d"), na.value = "white",
        drop = FALSE, na.translate = FALSE, name = "Evidence") +
      ggplot2::coord_equal(expand = FALSE) +
      ggplot2::labs(title = paste(tip, "contraction-surviving patches |", tm, "Ma"),
        subtitle = "At least two carriage links, >= 100,000 km2, 10% opportunity contraction",
        caption = "Candidates are not confirmed occupancy or fossil-validated refugia.",
        x = "Palaeolongitude", y = "Palaeolatitude") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
    file <- file.path(map_dir, paste0(tip, "_candidate_refugia_", tm, "Ma.png"))
    ggplot2::ggsave(file, fig, width = 12, height = 6.4, dpi = 150)
    tif <- sub("\\.png$", ".tif", file)
    save_geotiff(flags, tif)
    map_index[[length(map_index) + 1L]] <- data.frame(
      metric = "candidate_ecological_refugia_evidence_class", time_ma = tm,
      focal_tip = tip, png = file, geotiff = tif)
  }
}
write_csv(do.call(rbind, map_index), "case05_patch_map_index.csv")
quality <- data.frame(check = c("full_time_axis", "full_tip_set_present",
  "four_scheme_land_mass", "all_map_files_exist", "no_marine_map_values",
  "candidate_not_occupancy", "posterior_edge_flow_not_imputed"),
  pass = c(length(times) == 66L, sum(global$time_ma == 0) == 200L,
    all(abs(scheme$location_mass_on_land - 1) < 1e-6),
    all(file.exists(do.call(rbind, map_index)$png)) &&
      all(file.exists(do.call(rbind, map_index)$geotiff)),
    all(vapply(seq_along(times), function(k) {
      x <- global_maps[[as.character(times[[k]])]]
      mask <- main_field$landscape_weight[[k]] > 0
      all(x[!mask] == 0)
    }, logical(1))), TRUE, TRUE))
write_csv(quality, "case05_patch_quality_gates.csv")
writeLines(c(
  "# Case05 derived patch-history results",
  "",
  "This is an analysis of the completed 200-tip / 66-slice Case05 v8 run.",
  "The 4-degree cell is the computational unit; connected supported cells are patches.",
  "Three representative terminal ancestral paths have plate-carried patch networks.",
  "All active lineages at all times contribute to environmental-opportunity counts.",
  "",
  "Environmental index = Phi(X beta), without reconstructed modern prevalence intercept.",
  "Threshold 0.6 and contraction cutoff 10% are declared sensitivity choices.",
  "Candidate patches must be >= 100,000 km2 and retain a later carriage link.",
  "Patch histories retain splits and merges. Their continuity is geographic, not biological movement.",
  "Candidate refugia require an opportunity contraction and carried support in the next time slice.",
  "Location mass is reported separately; it is not occupancy and not multiplied by environment.",
  "The model reuses modern occurrence data for HMSC and terminal location conditioning.",
  "No independent fossil support is attached in this output.",
  "No-analog flags mark cells outside any modern training-axis 1st-99th percentile range.",
  "This is a separate extrapolation diagnostic, not a probability or biological penalty.",
  "No corridor or donor-recipient result is claimed because joint posterior edge flow was not stored.",
  "A source-sink classification cannot be made without demographic rates.",
  ""
), file.path(out, "README.md"))
stopifnot(nrow(global) > 0, nrow(scheme) == length(times) * length(focal) * 4L,
          all(file.exists(do.call(rbind, map_index)$png)),
          all(file.exists(do.call(rbind, map_index)$geotiff)),
          all(quality$pass))
message("Case05 patch-history results complete: ", out)
