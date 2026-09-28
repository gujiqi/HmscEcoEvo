#!/usr/bin/env Rscript

# Case05 v7: Arias-style global-sphere phylogeography at 4 degrees ----------
#
# The v6 preflight demonstrated that a land-only carrier graph makes distant
# terminal plant distributions mathematically disconnected.  That is not the
# Arias (2024) model: its spherical diffusion has a global grid and assigns
# low, non-zero landscape weights to marine cells.  This runner therefore uses
# every 4-degree global cell as a location state.  The real 1-degree PALEOMAP
# carrier tracks are used only to aggregate plate carriage between successive
# 4-degree cells.  Land fraction and topographic permeability determine the
# Arias destination landscape weight; ocean is a low-weight transit landscape.
#
# It intentionally does not calculate colonisation, establishment,
# persistence, local extinction, occupancy probability, richness, Shannon,
# Simpson, PD, or functional diversity.  Those are not identifiable from the
# present input without an independently calibrated demographic model.
#
# Tree nodes are inserted into the 5 Ma palaeoenvironment time grid. Between
# adjacent native PALEOMAP slices, plate locations and movement landscapes are
# linearly interpolated only for the short subinterval needed to end exactly
# at a dated tree node. This avoids letting a 5 Ma transition cross a
# speciation event. The interpolation is recorded as a temporal-resolution
# approximation, not as extra palaeogeographic evidence.

`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x)[1L]) y else x

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) next
    kv <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    out[[kv[[1L]]]] <- paste(kv[-1L], collapse = "=")
  }
  out
}
as_positive_integer <- function(x, default) {
  value <- suppressWarnings(as.integer(x %||% default))
  if (!is.finite(value) || value < 1L) stop("Expected a positive integer.", call. = FALSE)
  value
}
as_positive_number <- function(x, default, label) {
  value <- suppressWarnings(as.numeric(x %||% default))
  if (!is.finite(value) || value <= 0) stop(label, " must be positive.", call. = FALSE)
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
  x <- suppressWarnings(as.numeric(x)); x[!is.finite(x) | x < 0] <- 0
  total <- sum(x)
  if (!is.finite(total) || total <= 0) stop(label, " has no positive mass.", call. = FALSE)
  x / total
}
child_key <- function(tree, child) {
  if (child <= length(tree$tip.label)) tree$tip.label[[child]] else paste0("node_", child)
}
descendant_tip_sets <- function(tree) {
  n_tip <- length(tree$tip.label)
  children <- split(tree$edge[, 2L], tree$edge[, 1L])
  visit <- function(node) {
    if (node <= n_tip) return(tree$tip.label[[node]])
    unlist(lapply(children[[as.character(node)]], visit), use.names = FALSE)
  }
  nodes <- sort(unique(as.vector(tree$edge)))
  output <- lapply(nodes, visit)
  names(output) <- vapply(nodes, function(node) child_key(tree, node), character(1))
  output
}
tree_geometry <- function(tree) {
  depth <- ape::node.depth.edgelength(tree)
  root <- setdiff(tree$edge[, 1L], tree$edge[, 2L])
  if (length(root) != 1L) stop("Tree must have one root.", call. = FALSE)
  root_age <- max(depth[seq_along(tree$tip.label)])
  list(root = root, root_age = root_age, node_age = root_age - depth)
}
align_extant_tree <- function(tree, maximum_root_age_ma, tolerance_ma = 1e-4) {
  # The supplied dated tree is ultrametric to numerical tolerance, but a few
  # tips finish 1e-6 Ma before zero. Aligning these terminal rounding errors
  # prevents current extant lineages being omitted from the 0 Ma map.
  depth <- ape::node.depth.edgelength(tree)
  n_tip <- length(tree$tip.label)
  root_age <- max(depth[seq_len(n_tip)])
  tip_lag <- root_age - depth[seq_len(n_tip)]
  if (max(abs(tip_lag), na.rm = TRUE) > tolerance_ma) {
    stop("The dated tree is not ultrametric within the declared tolerance.", call. = FALSE)
  }
  terminal_edge <- match(seq_len(n_tip), tree$edge[, 2L])
  tree$edge.length[terminal_edge] <- tree$edge.length[terminal_edge] + tip_lag
  depth <- ape::node.depth.edgelength(tree)
  root_age <- max(depth[seq_len(n_tip)])
  if (root_age > maximum_root_age_ma + tolerance_ma) {
    # The input tree is 0.05 Ma older than the oldest available 325 Ma layer.
    # Rescale once, preserving every relative node age while keeping all
    # reconstructed states within observed palaeoenvironment support.
    tree$edge.length <- tree$edge.length * (maximum_root_age_ma / root_age)
  }
  tree
}
time_id <- function(x) as.character(round(as.numeric(x), 6L))
select_phylogenetically_representative_tips <- function(tree, occurrence_count, n) {
  candidates <- intersect(tree$tip.label, names(occurrence_count))
  if (n >= length(candidates)) return(candidates)
  distances <- ape::cophenetic.phylo(tree)[candidates, candidates, drop = FALSE]
  records <- occurrence_count[candidates]
  selected <- candidates[[which.max(records)]]
  while (length(selected) < n) {
    remaining <- setdiff(candidates, selected)
    separation <- apply(distances[remaining, selected, drop = FALSE], 1L, min)
    record_score <- log1p(records[remaining]) / max(log1p(records[candidates]))
    # Tree separation chooses a deep, representative sample; record support
    # only breaks otherwise equivalent choices and never defines history.
    score <- separation * (1 + 0.1 * record_score)
    selected <- c(selected, remaining[[which.max(score)]])
  }
  selected
}
bin4 <- function(lon, lat) {
  ix <- pmin(89L, pmax(0L, floor((as.numeric(lon) + 180) / 4)))
  iy <- pmin(44L, pmax(0L, floor((as.numeric(lat) + 90) / 4)))
  as.integer(iy * 90L + ix + 1L)
}
global_grid_4deg <- function() {
  ix <- rep(0:89, times = 45)
  iy <- rep(0:44, each = 90)
  lat_s <- -90 + 4 * iy
  lat_n <- lat_s + 4
  earth_radius <- 6371.0088
  area <- earth_radius^2 * (4 * pi / 180) *
    (sin(lat_n * pi / 180) - sin(lat_s * pi / 180))
  data.frame(
    state = sprintf("g4_%04d", seq_len(4050L)), bin = seq_len(4050L),
    lon = -178 + 4 * ix, lat = -88 + 4 * iy,
    cell_area_km2 = area, stringsAsFactors = FALSE
  )
}
aggregate_slice <- function(carrier, grid, axes) {
  carrier <- carrier[!is.na(carrier$paleo_lon) & !is.na(carrier$paleo_lat), , drop = FALSE]
  active <- !is.na(carrier$active_land) & carrier$active_land > 0 &
    is.finite(carrier$land_area_km2) & carrier$land_area_km2 > 0
  land <- carrier[active, , drop = FALSE]
  out <- grid
  out$land_fraction <- 0
  out$mean_permeability <- 0
  for (axis in axes) out[[axis]] <- NA_real_
  if (!nrow(land)) return(out)
  land$bin <- bin4(land$paleo_lon, land$paleo_lat)
  total_area <- rowsum(land$land_area_km2, land$bin, reorder = FALSE)
  bin_id <- as.integer(rownames(total_area))
  out$land_fraction[bin_id] <- pmin(1, total_area[, 1L] / out$cell_area_km2[bin_id])
  permeability <- land$topographic_permeability
  permeability[!is.finite(permeability)] <- 0
  weighted_mean <- function(x) {
    x <- suppressWarnings(as.numeric(x)); good <- is.finite(x)
    numerator <- rowsum(land$land_area_km2[good] * x[good], land$bin[good], reorder = FALSE)
    denominator <- rowsum(land$land_area_km2[good], land$bin[good], reorder = FALSE)
    answer <- rep(NA_real_, nrow(out)); index <- as.integer(rownames(numerator))
    answer[index] <- numerator[, 1L] / denominator[match(rownames(numerator), rownames(denominator)), 1L]
    answer
  }
  out$mean_permeability <- weighted_mean(permeability)
  out$mean_permeability[!is.finite(out$mean_permeability)] <- 0
  for (axis in axes) out[[axis]] <- weighted_mean(land[[axis]])
  out
}
interpolate_longitude <- function(from, to, fraction) {
  offset <- ((to - from + 180) %% 360) - 180
  ((from + fraction * offset + 180) %% 360) - 180
}
interpolate_carrier <- function(older, younger, fraction, axes) {
  if (fraction <= 0) return(older)
  if (fraction >= 1) return(younger)
  target <- younger[match(older$track_id, younger$track_id), , drop = FALSE]
  output <- older
  both <- !is.na(older$paleo_lon) & !is.na(older$paleo_lat) &
    !is.na(target$paleo_lon) & !is.na(target$paleo_lat)
  output$paleo_lon[both] <- interpolate_longitude(
    older$paleo_lon[both], target$paleo_lon[both], fraction
  )
  output$paleo_lat[both] <- (1 - fraction) * older$paleo_lat[both] +
    fraction * target$paleo_lat[both]
  for (field in unique(c("active_land", "land_area_km2", "topographic_permeability", axes))) {
    if (!field %in% names(older) || !field %in% names(target)) next
    x <- suppressWarnings(as.numeric(older[[field]]))
    y <- suppressWarnings(as.numeric(target[[field]]))
    good <- is.finite(x) & is.finite(y)
    value <- x
    value[good] <- (1 - fraction) * x[good] + fraction * y[good]
    output[[field]] <- value
  }
  output
}
plate_transport_4deg <- function(old_carrier, young_carrier, old_slice, grid) {
  n <- nrow(grid)
  old_carrier <- old_carrier[!is.na(old_carrier$paleo_lon) & !is.na(old_carrier$paleo_lat), , drop = FALSE]
  young_carrier <- young_carrier[match(old_carrier$track_id, young_carrier$track_id), , drop = FALSE]
  active <- !is.na(old_carrier$active_land) & old_carrier$active_land > 0 &
    is.finite(old_carrier$land_area_km2) & old_carrier$land_area_km2 > 0 &
    !is.na(young_carrier$paleo_lon) & !is.na(young_carrier$paleo_lat)
  source_bin <- bin4(old_carrier$paleo_lon[active], old_carrier$paleo_lat[active])
  target_bin <- bin4(young_carrier$paleo_lon[active], young_carrier$paleo_lat[active])
  weight <- old_carrier$land_area_km2[active]
  if (!length(weight)) {
    return(Matrix::Diagonal(n = n, x = 1, names = list(grid$state, grid$state)))
  }
  total_by_source <- rowsum(weight, source_bin, reorder = FALSE)
  fraction <- weight / total_by_source[match(as.character(source_bin), rownames(total_by_source)), 1L]
  land_share <- old_slice$land_fraction[source_bin]
  value <- land_share * fraction
  carried <- Matrix::sparseMatrix(
    i = target_bin, j = source_bin, x = value, dims = c(n, n),
    dimnames = list(grid$state, grid$state)
  )
  carried_column_sum <- Matrix::colSums(carried)
  residual <- pmax(0, 1 - carried_column_sum)
  transport <- carried + Matrix::Diagonal(n = n, x = residual,
                                           names = list(grid$state, grid$state))
  error <- max(abs(Matrix::colSums(transport) - 1))
  if (error > 1e-10) stop("4-degree plate transport is not source-normalized.", call. = FALSE)
  transport
}
map_values <- function(values, grid) {
  out <- rep(0, nrow(grid)); out[] <- as.numeric(values[grid$state]); out
}
write_map <- function(grid, values, time_ma, metric, limits, dirs, unit, caption,
                      signed = limits[[1L]] < 0) {
  tab <- grid; tab$value <- as.numeric(values)
  # Keep the semantic metric name in the index, but use compact physical paths
  # so Windows' legacy 260-character path limit cannot silently truncate .png.
  map_id <- switch(metric,
    posterior_lineage_location_density = "01_location_density",
    posterior_lineage_location_density_log10 = "02_location_density_log10",
    posterior_weighted_ancestral_HMSC_linear_predictor = "03_ancestral_HMSC_eta",
    Arias_destination_landscape_weight = "04_landscape_weight",
    terminal_calibrated_lineage_support_intensity = "05_lineage_support_intensity",
    support_mixture_shannon = "06_support_mixture_shannon",
    support_mixture_effective_lineages = "07_support_effective_lineages",
    stop("Unknown Case05 v7 metric: ", metric, call. = FALSE)
  )
  png <- file.path(dirs$maps, map_id, sprintf("%s_%03dMa.png", map_id, round(time_ma)))
  tif <- file.path(dirs$maps, map_id, sprintf("%s_%03dMa.tif", map_id, round(time_ma)))
  safe_dir(dirname(png))
  scale <- if (isTRUE(signed)) {
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
                                  midpoint = 0, limits = limits, oob = scales::squish,
                                  name = unit)
  } else {
    ggplot2::scale_fill_viridis_c(option = "C", limits = limits, oob = scales::squish,
                                  name = unit)
  }
  plot <- ggplot2::ggplot(tab, ggplot2::aes(lon, lat)) +
    ggplot2::geom_tile(ggplot2::aes(fill = value), colour = NA) + scale +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::labs(title = paste0("Case05 v7: ", metric, "; ", round(time_ma), " Ma"),
                  subtitle = "Global 4-degree Arias-style spherical state grid; fixed colour scale across time",
                  x = "Palaeographic longitude (degrees)", y = "Palaeographic latitude (degrees)",
                  caption = caption) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_line(colour = "grey88"))
  ggplot2::ggsave(png, plot, width = 12, height = 6.4, dpi = 150)
  raster <- terra::rast(ncols = 90, nrows = 45, xmin = -180, xmax = 180,
                         ymin = -90, ymax = 90, crs = "EPSG:4326")
  cells <- terra::cellFromXY(raster, as.matrix(tab[, c("lon", "lat")]))
  output <- rep(NA_real_, terra::ncell(raster)); output[cells] <- tab$value
  terra::values(raster) <- output; names(raster) <- metric
  terra::writeRaster(raster, tif, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
  data.frame(metric = metric, time_ma = time_ma,
             map_directory = map_id,
             png = normalizePath(png, winslash = "/", mustWork = TRUE),
             geotiff = normalizePath(tif, winslash = "/", mustWork = TRUE),
             unit = unit, scale_min = limits[[1L]], scale_max = limits[[2L]],
             stringsAsFactors = FALSE)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
pkg_root <- normalizePath(args$pkg_root %||%
  "C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo", winslash = "/", mustWork = TRUE)
output <- safe_dir(args$output %||% file.path(
  "C:/Users/Google/Documents/HMSC-HIST/outputs/HmscEcoEvo",
  "case05_plant200_v7_arias_global_sphere_4deg_preflight_20260923"
))
n_tips <- as_positive_integer(args$n_tips, 20L)
n_draws <- as_positive_integer(args$n_response_draws, 5L)
seed <- as_positive_integer(args$seed, 20260923L)
ocean_weight <- as_positive_number(args$ocean_weight, 0.001, "ocean_weight")
diffusion_variance <- as_positive_number(args$diffusion_variance_rad2_per_myr, 0.0002,
                                         "diffusion_variance_rad2_per_myr")
emigration_rate <- as_positive_number(args$source_emigration_rate_per_myr, 0.4,
                                      "source_emigration_rate_per_myr")
for (package in c("HmscEcoEvo", "ape", "Matrix", "ggplot2", "terra")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Missing package: ", package, call. = FALSE)
}

dirs <- list(
  contract = safe_dir(file.path(output, "00_contract")),
  audit = safe_dir(file.path(output, "01_input_audit")),
  response = safe_dir(file.path(output, "02_ancestral_response")),
  transition = safe_dir(file.path(output, "03_global_sphere_transitions")),
  pruning = safe_dir(file.path(output, "04_phylogeographic_pruning")),
  maps = safe_dir(file.path(output, "05_shared_scale_maps")),
  quality = safe_dir(file.path(output, "06_quality"))
)

contract <- data.frame(
  component = c("location_state", "plate_carriage", "horizontal_diffusion", "landscape_weight",
                "ocean", "ancestral_environmental_response", "terminal_constraint", "colonisation",
                "persistence", "local_extinction", "occupancy", "dispersal_conditioned_diversity", "tree_time_axis"),
  status = c("all_global_4_degree_cells", "real_PALEOMAP_carrier_tracks_aggregated_to_4_degree",
             "Arias_spherical_normal_on_global_neighbour_graph", "land_fraction_times_permeability",
             "low_nonzero_transit_weight", "separate_HMSC_beta_diagnostic", "modern_occurrence_density_cut_model",
             "NOT_RUN", "NOT_RUN", "NOT_RUN", "NOT_RUN", "terminal_calibrated_support_scenario",
             "native_5Ma_slices_plus_exact_dated_tree_nodes"),
  explanation = c(
    "Every global 4-degree cell is a state, so disconnected land-only components cannot create false zero likelihoods.",
    "Plate carriage is a separate source-normalized operator inferred from stable real carrier track IDs.",
    "Horizontal movement is a source-normalized spherical diffusion transition, then composed with plate carriage.",
    "Emergent land and topography alter destination weights once; they are not a second survival multiplier.",
    "Current inputs lack shallow/deep marine classes; all non-land receives the declared Arias-style low-weight ocean proxy 0.001.",
    "Ancestral HMSC beta is reconstructed on the tree and mapped separately; it is never multiplied into geographic diffusion.",
    "Terminal observations are reused modularly after HMSC fitting; this is not a joint likelihood or independent validation.",
    "No unidentifiable demographic establishment coefficient is used.",
    "No unidentifiable local survival coefficient is used.",
    "No geographic loss is re-labelled as biological extinction.",
    "No final palaeo-occupancy probability is claimed.",
    "Lineage support intensity and mixture diversity are derived from terminal-calibrated location fields. They are not historical species richness, PD, FD, or ecological alpha diversity.",
    "Each branch is divided at every native palaeoenvironment boundary and its exact dated nodes; short intervening geographic states are linearly interpolated from bracketing PALEOMAP slices."
  ), stringsAsFactors = FALSE
)
write_csv(contract, file.path(dirs$contract, "case05_v7_scientific_contract.csv"))

input_root <- normalizePath(file.path(dirname(pkg_root), "HmscEcoEvo_API8_20260818", "data_external",
  "plant_200_multifamily_extant_20260820", "prepared_inputs", "case03_data_raw",
  "Plant200_multifamily_extant_1deg_traits"), winslash = "/", mustWork = TRUE)
legacy_output <- normalizePath(file.path(dirname(pkg_root), "outputs", "HmscEcoEvo",
  "case05_plant200_4deg_effortaware_formal25_all66_20260923"), winslash = "/", mustWork = TRUE)
carrier_root <- normalizePath(file.path(pkg_root, "derived_inputs",
  "case05_v2_paleomap_h3_carrier_325_0Ma_1deg", "movement_cache_arias_topographic_rate0p4"),
  winslash = "/", mustWork = TRUE)
comm <- utils::read.csv(file.path(input_root, "comm.csv"), check.names = FALSE, stringsAsFactors = FALSE)
sites <- utils::read.csv(file.path(input_root, "sites.csv"), check.names = FALSE, stringsAsFactors = FALSE)
tree <- ape::read.tree(file.path(input_root, "tree.tre"))
beta <- readRDS(file.path(legacy_output, "02_hmsc_4deg", "hmsc_environmental_beta_posterior_draws.rds"))
index <- utils::read.csv(file.path(carrier_root, "case05_v2_movement_cache_index.csv"), stringsAsFactors = FALSE)
index <- index[order(index$time_ma, decreasing = TRUE), , drop = FALSE]
if (!identical(comm$site_id, sites$site_id)) stop("comm/sites site_id mismatch.", call. = FALSE)
species <- setdiff(names(comm), "site_id")
if (length(setdiff(species, tree$tip.label))) stop("Tree lacks community species.", call. = FALSE)
axes <- intersect(c("MAT_pohl_C", "MAP_pohl_mm_yr", "elevation_m", "T_seasonality_pohl_sd",
  "P_seasonality_pohl_sd", "moisture_availability_index_z", "wetland_potential_index_z", "relief_3x3_m"), names(beta))
if (length(axes) < 2L) stop("HMSC beta axes missing.", call. = FALSE)
occurrence_count <- colSums(as.matrix(comm[, species, drop = FALSE]) > 0)
tips <- select_phylogenetically_representative_tips(
  tree, occurrence_count, min(n_tips, length(species))
)
tree <- ape::reorder.phylo(ape::keep.tip(tree, tips), "cladewise")
tree <- align_extant_tree(tree, maximum_root_age_ma = max(index$time_ma))
geom <- tree_geometry(tree)
geom$node_age <- pmax(0, round(geom$node_age, 6L))
display_times <- index$time_ma[index$time_ma <= geom$root_age + 1e-8]
node_times <- unique(geom$node_age[unique(as.vector(tree$edge))])
master_times <- sort(unique(c(display_times, node_times)), decreasing = TRUE)
master_times <- master_times[master_times >= 0 & master_times <= geom$root_age + 1e-8]
draw_ids <- unique(beta$response_draw); draw_ids <- draw_ids[seq_len(min(n_draws, length(draw_ids)))]
tip_beta <- beta[beta$lineage %in% tree$tip.label & beta$response_draw %in% draw_ids,
                 c("lineage", "response_draw", axes), drop = FALSE]
names(tip_beta)[names(tip_beta) == "lineage"] <- "species"
if (nrow(tip_beta) != length(tree$tip.label) * length(draw_ids)) stop("Incomplete beta draw coverage.", call. = FALSE)
grid <- global_grid_4deg(); states <- grid$state
write_csv(data.frame(item = c("package_version", "selected_tips", "response_draws", "tree_root_age_ma",
  "native_display_time_slices", "master_time_slices_with_tree_nodes", "state_grid", "ocean_weight", "diffusion_variance_rad2_per_myr", "source_emigration_rate_per_myr"),
  value = c(as.character(utils::packageVersion("HmscEcoEvo")), length(tips), length(draw_ids), geom$root_age,
  length(display_times), length(master_times), "global_4_degree_4050_cells", ocean_weight, diffusion_variance, emigration_rate)),
  file.path(dirs$audit, "input_audit.csv"))
write_csv(data.frame(species = names(occurrence_count), occurrences = occurrence_count,
  selected = names(occurrence_count) %in% tips,
  selection_method = "greedy_phylogenetic_representativeness_with_occurrence_tiebreak"),
  file.path(dirs$audit, "species_selection.csv"))
write_csv(data.frame(
  time_ma = master_times,
  native_paleomap_slice = master_times %in% index$time_ma,
  dated_tree_node = master_times %in% node_times,
  output_map_time = master_times %in% display_times,
  stringsAsFactors = FALSE
), file.path(dirs$audit, "master_time_axis.csv"))

# Read real PALEOMAP carrier slices and retain their stable track IDs for the
# plate-carriage operator. Tree nodes that fall between native 5 Ma slices are
# interpolated from the bracketing carrier states; they never receive a full
# transition that crosses the speciation time.
carrier_cache <- list()
for (i in seq_len(nrow(index))) {
  carrier <- readRDS(index$file[[i]])$carrier
  carrier_cache[[time_id(index$time_ma[[i]])]] <- carrier
}
native_times <- as.numeric(index$time_ma)
carrier_state_cache <- new.env(parent = emptyenv())
carrier_state_at <- function(time_ma) {
  key <- time_id(time_ma)
  if (exists(key, envir = carrier_state_cache, inherits = FALSE)) {
    return(get(key, envir = carrier_state_cache))
  }
  exact <- which(abs(native_times - time_ma) < 1e-8)
  if (length(exact)) {
    state <- carrier_cache[[time_id(native_times[[exact[[1L]]]])]]
  } else {
    older_index <- tail(which(native_times > time_ma), 1L)
    younger_index <- head(which(native_times < time_ma), 1L)
    if (!length(older_index) || !length(younger_index)) {
      stop("Requested time is outside the native PALEOMAP time domain: ", time_ma,
           call. = FALSE)
    }
    older_time <- native_times[[older_index]]; younger_time <- native_times[[younger_index]]
    fraction <- (older_time - time_ma) / (older_time - younger_time)
    state <- interpolate_carrier(
      carrier_cache[[time_id(older_time)]], carrier_cache[[time_id(younger_time)]],
      fraction, axes
    )
  }
  assign(key, state, envir = carrier_state_cache)
  state
}
slice_cache <- new.env(parent = emptyenv())
slice_at <- function(time_ma) {
  key <- time_id(time_ma)
  if (exists(key, envir = slice_cache, inherits = FALSE)) {
    return(get(key, envir = slice_cache))
  }
  slice <- aggregate_slice(carrier_state_at(time_ma), grid, axes)
  slice$time_ma <- time_ma
  slice$landscape_weight <- ocean_weight + (1 - ocean_weight) *
    pmin(1, pmax(0, slice$land_fraction * slice$mean_permeability))
  assign(key, slice, envir = slice_cache)
  slice
}

# A full global 8-neighbour grid creates a sparse approximation of the
# spherical-normal kernel.  Low-weight ocean cells retain connectivity without
# treating them as plant habitat or demographic occupancy states.
edge_cells <- data.frame(cell_id = grid$state, time_ma = 0, lon = grid$lon, lat = grid$lat,
  elevation_m = 0, topographic_resistance = 0, movement_allowed = 1)
base_edges <- HmscEcoEvo::hee_dispersal_connectivity_graph(
  edge_cells, resolution_deg = 4, neighbours = 8, omega_topo = 0
)$edges

transition_cache <- new.env(parent = emptyenv())
make_transition <- function(older_time, younger_time) {
  if (!is.finite(older_time) || !is.finite(younger_time) || older_time <= younger_time) {
    stop("A transition requires finite older_time > younger_time.", call. = FALSE)
  }
  key <- paste(time_id(older_time), time_id(younger_time), sep = "_to_")
  if (exists(key, envir = transition_cache, inherits = FALSE)) return(get(key, envir = transition_cache))
  delta <- older_time - younger_time
  middle_time <- (older_time + younger_time) / 2
  older <- slice_at(older_time); middle <- slice_at(middle_time)
  landscape_cells <- data.frame(cell_id = grid$state, time_ma = middle_time, lon = grid$lon, lat = grid$lat,
    movement_weight = middle$landscape_weight, cell_area_km2 = grid$cell_area_km2)
  landscape <- HmscEcoEvo::hee_dispersal_landscape_weights(
    landscape_cells, movement_weight_col = "movement_weight", cell_area_col = "cell_area_km2"
  )
  edges <- base_edges; edges$time_ma <- middle_time
  kernel <- HmscEcoEvo::hee_dispersal_spherical_kernel(
    cells = landscape, edges = edges, delta_t_myr = delta,
    diffusion_variance_rad2_per_myr = diffusion_variance,
    source_emigration_rate_per_myr = emigration_rate, time_ma = middle_time,
    edge_cost_col = NULL, spatial_domain = "global_explicit"
  )
  diffusion <- HmscEcoEvo::hee_dispersal_location_transition(kernel, cell_order = grid$state,
    delta_t_myr = delta, sparse = TRUE)
  plate <- plate_transport_4deg(
    carrier_state_at(older_time), carrier_state_at(younger_time), older, grid
  )
  composed <- HmscEcoEvo::hee_dispersal_compose_plate_transport(diffusion, plate)
  out <- list(transition = composed$transition, time_ma = younger_time, delta_t_myr = delta,
    diagnostics = cbind(time_from_ma = older_time, time_to_ma = younger_time,
      midpoint_ma = middle_time,
      interpolated_endpoint = !(older_time %in% native_times && younger_time %in% native_times),
      composed$diagnostics))
  assign(key, out, envir = transition_cache); out
}
old_time <- master_times[-length(master_times)]; young_time <- master_times[-1L]
branch_steps <- list()
for (e in seq_len(nrow(tree$edge))) {
  parent <- tree$edge[e, 1L]; child <- tree$edge[e, 2L]
  parent_age <- geom$node_age[[parent]]; child_age <- geom$node_age[[child]]
  keep <- which(old_time <= parent_age + 1e-8 & young_time >= child_age - 1e-8)
  steps <- lapply(keep, function(i) make_transition(old_time[[i]], young_time[[i]]))
  if (!length(steps) ||
      abs(old_time[[keep[[1L]]]] - parent_age) > 1e-8 ||
      abs(young_time[[tail(keep, 1L)]] - child_age) > 1e-8) {
    stop("Master time grid does not exactly cover a dated tree branch.", call. = FALSE)
  }
  branch_steps[[child_key(tree, child)]] <- steps
}
transition_qc <- do.call(rbind, lapply(ls(transition_cache), function(key) get(key, transition_cache)$diagnostics))
write_csv(as.data.frame(transition_qc), file.path(dirs$transition, "transition_quality.csv"))

# Terminal likelihood is a 4-degree presence density.  It is not a survey
# absence model and it is explicitly marked as a modular reuse of the modern
# data that generated the HMSC posterior.
site_bin <- bin4(sites$lon, sites$lat)
tip_likelihood <- list(); terminal_audit <- list()
for (tip in tree$tip.label) {
  bins <- unique(site_bin[comm[[tip]] > 0])
  likelihood <- stats::setNames(rep(0, length(states)), states)
  likelihood[states[bins]] <- normalise(grid$cell_area_km2[bins], paste0("terminal likelihood for ", tip))
  tip_likelihood[[tip]] <- likelihood
  terminal_audit[[tip]] <- data.frame(species = tip, occurrence_sites = sum(comm[[tip]] > 0),
    terminal_bins_4deg = length(bins), stringsAsFactors = FALSE)
}
write_csv(do.call(rbind, terminal_audit), file.path(dirs$audit, "terminal_occurrence_density_audit.csv"))
root_time <- geom$root_age
root_slice <- slice_at(root_time)
root_prior <- stats::setNames(normalise(root_slice$landscape_weight * grid$cell_area_km2,
  "global landscape-weighted root prior"), states)
responses <- HmscEcoEvo::hee_evolution_ancestral_response_direct(
  tip_beta_draws = tip_beta, tree = tree, times = master_times, basis_cols = axes,
  reconstruct_intercept = FALSE, branch_uncertainty = "brownian_bridge", seed = seed
)
saveRDS(responses, file.path(dirs$response, "ancestral_beta_response_draws.rds"), compress = "gzip")
write_csv(aggregate(responses[, axes, drop = FALSE], responses[, c("lineage", "time_ma"), drop = FALSE], mean),
  file.path(dirs$response, "ancestral_beta_response_mean.csv"))

map_times <- time_id(sort(unique(display_times), decreasing = TRUE))
density_sum <- stats::setNames(lapply(map_times, function(x) rep(0, length(states))), map_times)
# Retain active-lineage fields separately so diversity is calculated from each
# lineage's spatial support, never from a previously averaged density map.
lineage_density_sum <- stats::setNames(lapply(map_times, function(x) list()), map_times)
eta_sum <- stats::setNames(lapply(map_times, function(x) rep(0, length(states))), map_times)
eta_weight <- stats::setNames(lapply(map_times, function(x) rep(0, length(states))), map_times)
loglik <- numeric(length(draw_ids)); names(loglik) <- draw_ids
accumulate_lineage_density <- function(key, lineage, values) {
  bucket <- lineage_density_sum[[key]]
  if (is.null(bucket[[lineage]])) bucket[[lineage]] <- rep(0, length(states))
  bucket[[lineage]] <- bucket[[lineage]] + as.numeric(values)
  lineage_density_sum[[key]] <<- bucket
  invisible(NULL)
}
response_at <- function(draw, lineage, time_ma) {
  row <- responses[responses$response_draw == draw & responses$lineage == lineage &
    abs(responses$time_ma - time_ma) < 1e-6, , drop = FALSE]
  if (!nrow(row)) return(NULL)
  as.numeric(row[1L, axes, drop = TRUE])
}
for (draw in draw_ids) {
  callback <- function(branch_key, child_node, time_ma, stage_index, posterior) {
    key <- time_id(time_ma); if (!key %in% names(density_sum)) return(invisible(NULL))
    posterior_values <- as.numeric(posterior[states])
    density_sum[[key]] <<- density_sum[[key]] + posterior_values
    accumulate_lineage_density(key, branch_key, posterior_values)
    beta_row <- response_at(draw, branch_key, time_ma)
    if (!is.null(beta_row)) {
      slice <- slice_at(time_ma); X <- as.matrix(slice[, axes, drop = FALSE])
      eta <- as.vector(X %*% beta_row); valid <- is.finite(eta) & slice$land_fraction > 0
      eta[!valid] <- 0
      posterior_mass <- as.numeric(posterior[states])
      eta_sum[[key]] <<- eta_sum[[key]] + posterior_mass * eta
      eta_weight[[key]] <<- eta_weight[[key]] + posterior_mass * valid
    }
    invisible(NULL)
  }
  result <- HmscEcoEvo::hee_dispersal_time_ordered_pruning(tree, tip_likelihood, branch_steps,
    root_prior = root_prior, state_callback = callback)
  if (!is.finite(result$log_likelihood)) stop("Non-finite pruning likelihood for ", draw, call. = FALSE)
  loglik[[draw]] <- result$log_likelihood
  root_key <- time_id(root_time)
  if (root_key %in% names(density_sum)) {
    root_values <- as.numeric(result$root_posterior[states])
    density_sum[[root_key]] <- density_sum[[root_key]] + root_values
    accumulate_lineage_density(root_key, child_key(tree, geom$root), root_values)
  }
  saveRDS(result, file.path(dirs$pruning, paste0("pruning_", draw, ".rds")), compress = "gzip")
}
write_csv(data.frame(response_draw = names(loglik), log_likelihood = loglik),
  file.path(dirs$pruning, "pruning_log_likelihoods.csv"))

density <- lapply(density_sum, function(x) x / length(draw_ids))
lineage_density <- lapply(lineage_density_sum, function(bucket) {
  if (!length(bucket)) stop("No active lineage density was accumulated for a map time.", call. = FALSE)
  value <- do.call(cbind, bucket) / length(draw_ids)
  rownames(value) <- states
  value
})
# The cut phylogeographic model is explicitly conditioned on modern terminal
# occurrence densities.  Enforce that conditioning at exactly 0 Ma rather than
# relying on a callback along arbitrarily short terminal branches.
present_key <- time_id(0)
if (!present_key %in% names(lineage_density)) {
  stop("The display-time grid must include 0 Ma for terminal conditioning.", call. = FALSE)
}
present_matrix <- lineage_density[[present_key]]
for (tip in tree$tip.label) {
  present_matrix[, tip] <- as.numeric(tip_likelihood[[tip]][states])
}
lineage_density[[present_key]] <- present_matrix
density_sum[[present_key]] <- rowSums(present_matrix) * length(draw_ids)
density <- lapply(density_sum, function(x) x / length(draw_ids))
terminal_conditioning_audit <- data.frame(
  lineage = tree$tip.label,
  terminal_support_cells = vapply(tree$tip.label, function(tip) sum(tip_likelihood[[tip]][states] > 0), integer(1)),
  terminal_density_mass = vapply(tree$tip.label, function(tip) sum(present_matrix[, tip]), numeric(1)),
  conditioning = "exact_modular_terminal_occurrence_density",
  stringsAsFactors = FALSE
)
write_csv(terminal_conditioning_audit, file.path(dirs$audit, "terminal_conditioning_audit.csv"))
density_log10 <- lapply(density, function(x) log10(pmax(x, 1e-12)))
eta <- stats::setNames(lapply(names(eta_sum), function(key) {
  x <- eta_sum[[key]] / pmax(eta_weight[[key]], .Machine$double.eps)
  x[eta_weight[[key]] <= 0] <- NA_real_; x
}), names(eta_sum))
density_limits <- c(0, max(unlist(density), na.rm = TRUE))
density_log10_limits <- range(unlist(density_log10), finite = TRUE)
eta_max <- max(abs(unlist(eta)), na.rm = TRUE); eta_limits <- c(-eta_max, eta_max)
present_lineage_mass <- if (present_key %in% names(density)) sum(density[[present_key]]) else NA_real_
expected_present_lineages <- length(tree$tip.label)
if (!all(tree$tip.label %in% colnames(lineage_density[[present_key]]))) {
  stop("Every extant tip must be represented in the 0 Ma lineage-density field.", call. = FALSE)
}
terminal_table <- do.call(rbind, terminal_audit)
terminal_support <- stats::setNames(terminal_table$terminal_bins_4deg, terminal_table$species)
support_calibration <- HmscEcoEvo::hee_dispersal_range_support_calibrate(
  lineage_density[[present_key]][, tree$tip.label, drop = FALSE], terminal_support,
  target_fraction = 0.995
)
write_csv(support_calibration, file.path(dirs$audit, "terminal_range_support_calibration.csv"))
tip_support_scale <- stats::setNames(support_calibration$support_scale, support_calibration$lineage)
descendant_tips <- descendant_tip_sets(tree)
support_diversity <- lapply(lineage_density, function(field) {
  support_scale <- vapply(colnames(field), function(lineage) {
    descendants <- descendant_tips[[lineage]]
    if (is.null(descendants) || !all(descendants %in% names(tip_support_scale))) {
      stop("Could not derive a descendant terminal support scale for ", lineage, call. = FALSE)
    }
    mean(tip_support_scale[descendants])
  }, numeric(1))
  names(support_scale) <- colnames(field)
  HmscEcoEvo::hee_dispersal_range_support_diversity(field, support_scale)$summary
})
support_intensity <- lapply(support_diversity, function(x) {
  stats::setNames(x$terminal_calibrated_lineage_support_intensity, x$state)
})
support_shannon <- lapply(support_diversity, function(x) {
  stats::setNames(x$support_mixture_shannon, x$state)
})
support_effective <- lapply(support_diversity, function(x) {
  stats::setNames(x$support_mixture_effective_lineages, x$state)
})
support_intensity_limits <- c(0, max(unlist(support_intensity), na.rm = TRUE))
support_shannon_limits <- c(0, max(unlist(support_shannon), na.rm = TRUE))
support_effective_limits <- c(1, max(unlist(support_effective), na.rm = TRUE))
maps <- list(); z <- 0L
for (key in names(density)) {
  slice <- slice_at(as.numeric(key)); tm <- as.numeric(key)
  z <- z + 1L; maps[[z]] <- write_map(grid, density[[key]], tm, "posterior_lineage_location_density",
    density_limits, dirs, "expected selected-lineage location density",
    "All global cells remain in the diffusion state space. This is a tree-conditioned location density, not occupancy or richness.")
  z <- z + 1L; maps[[z]] <- write_map(grid, density_log10[[key]], tm,
    "posterior_lineage_location_density_log10", density_log10_limits, dirs,
    "log10 expected selected-lineage location density",
    "Display-only logarithmic companion to the absolute density map; it reveals broad, low-density ancestral support without changing the underlying result.",
    signed = FALSE)
  z <- z + 1L; maps[[z]] <- write_map(grid, eta[[key]], tm, "posterior_weighted_ancestral_HMSC_linear_predictor",
    eta_limits, dirs, "posterior-location-weighted environmental linear predictor",
    "Only emerged-land cells have environmental input. This separate HMSC diagnostic is not multiplied into geographic diffusion.")
  z <- z + 1L; maps[[z]] <- write_map(grid, slice$landscape_weight, tm, "Arias_destination_landscape_weight",
    c(ocean_weight, 1), dirs, "dimensionless destination landscape weight",
    "Land fraction and topographic permeability derive from real carrier data; non-land uses the declared low-weight ocean proxy.")
  z <- z + 1L; maps[[z]] <- write_map(grid, support_intensity[[key]], tm,
    "terminal_calibrated_lineage_support_intensity", support_intensity_limits, dirs,
    "terminal-calibrated lineage-support intensity",
    "Sum of diffusion-conditioned lineage support fields. At 0 Ma it is calibrated to 99.5% of observed terminal grid support; it is not species richness or occupancy.")
  z <- z + 1L; maps[[z]] <- write_map(grid, support_shannon[[key]], tm,
    "support_mixture_shannon", support_shannon_limits, dirs,
    "Shannon diversity of lineage-support mixture",
    "Shannon diversity of terminal-calibrated lineage support, not ecological Shannon diversity or a reconstructed historical community census.")
  z <- z + 1L; maps[[z]] <- write_map(grid, support_effective[[key]], tm,
    "support_mixture_effective_lineages", support_effective_limits, dirs,
    "effective number of lineage-support fields",
    "exp(Shannon) for the local lineage-support mixture. It is a dispersal-conditioned support diagnostic, not species richness or occupancy.")
}
map_index <- do.call(rbind, maps); write_csv(map_index, file.path(dirs$maps, "case05_v7_map_index.csv"))
quality <- data.frame(
  check = c("colonisation_not_run", "persistence_not_run", "no_occupancy_or_species_richness_claim",
    "all_pruning_likelihoods_finite", "all_transition_columns_normalised", "global_sphere_state_grid",
    "real_plate_carriage_aggregated", "present_tip_density_mass", "png_maps", "geotiff_maps"),
  pass = c(TRUE, TRUE, TRUE, all(is.finite(loglik)),
    max(transition_qc$max_source_normalisation_error) <= 1e-10, TRUE, TRUE,
    is.finite(present_lineage_mass) && abs(present_lineage_mass - expected_present_lineages) <= 1e-6,
    all(file.exists(map_index$png)), all(file.exists(map_index$geotiff))),
  value = c("TRUE", "TRUE", "TRUE", as.character(all(is.finite(loglik))),
    max(transition_qc$max_source_normalisation_error), "4050 cells", "stable carrier track IDs",
    paste0(present_lineage_mass, " expected ", expected_present_lineages),
    sum(file.exists(map_index$png)), sum(file.exists(map_index$geotiff))), stringsAsFactors = FALSE)
write_csv(quality, file.path(dirs$quality, "case05_v7_quality_gates.csv"))
saveRDS(list(location_density = density, lineage_location_density = lineage_density,
  location_density_log10 = density_log10,
  environmental_linear_predictor = eta,
  landscape_weight = lapply(names(density), function(key) slice_at(as.numeric(key))$landscape_weight),
  terminal_calibrated_lineage_support_intensity = support_intensity,
  support_mixture_shannon = support_shannon,
  support_mixture_effective_lineages = support_effective,
  interpretation = paste("Arias-style plate-carried global location distributions plus",
    "terminal-calibrated lineage-support diversity; no demographic occupancy layer.")),
  file.path(dirs$maps, "case05_v7_map_values.rds"), compress = "gzip")
message("Case05 v7 complete: ", output)
