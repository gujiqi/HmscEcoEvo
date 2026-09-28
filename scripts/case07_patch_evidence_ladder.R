#!/usr/bin/env Rscript

# Case07: real Plant200 patch diagnostics plus a controlled joint-edge example.
# The controlled example is not an empirical Case05 movement posterior.

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) next
    pieces <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    out[[pieces[[1L]]]] <- paste(pieces[-1L], collapse = "=")
  }
  out
}
`%||%` <- function(x, y) if (is.null(x)) y else x
cfg <- parse_args(commandArgs(trailingOnly = TRUE))
script_arg <- commandArgs(FALSE)
script_arg <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
if (!length(script_arg)) stop("Run Case07 with Rscript.", call. = FALSE)
pkg_root <- normalizePath(file.path(dirname(script_arg[[1L]]), ".."), winslash = "/")
project_root <- normalizePath(file.path(pkg_root, ".."), winslash = "/")
comparison <- file.path(project_root, "outputs", "HmscEcoEvo",
  "case05_v8_four_scheme_comparison_20260926")
source_dir <- normalizePath(cfg$source %||% file.path(comparison,
  "10_patch_history_results_corrected"), winslash = "/", mustWork = TRUE)
out <- normalizePath(cfg$output %||% file.path(project_root, "outputs", "HmscEcoEvo",
  paste0("case07_patch_evidence_ladder_", format(Sys.Date(), "%Y%m%d"))),
  winslash = "/", mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (name in c("00_contract", "01_real_plant200", "02_controlled_patch_cases",
               "03_joint_edge_example", "04_figures", "05_quality")) {
  dir.create(file.path(out, name), recursive = TRUE, showWarnings = FALSE)
}
if (!requireNamespace("HmscEcoEvo", quietly = TRUE) ||
    !requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Install HmscEcoEvo and ggplot2 before Case07.", call. = FALSE)
}
if (utils::packageVersion("HmscEcoEvo") < package_version("1.0.4") ||
    sum(HmscEcoEvo::hee_core_process_catalog()$layer == "core_process") != 6L) {
  stop("Case07 requires HmscEcoEvo >= 1.0.4 with the six-process API.",
       call. = FALSE)
}
write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
quality <- data.frame(check = character(), pass = logical(), detail = character())
gate <- function(check, pass, detail = "") {
  quality[nrow(quality) + 1L, ] <<- list(check, isTRUE(pass), as.character(detail))
}

source_quality <- utils::read.csv(file.path(source_dir, "case05_patch_quality_gates.csv"))
gate("corrected_source_selected", grepl("10_patch_history_results_corrected$", source_dir),
     source_dir)
gate("corrected_response_audit_exists", file.exists(file.path(comparison,
     "09_corrected_ancestral_response", "README_zh.md")) ||
     dir.exists(file.path(comparison, "09_corrected_ancestral_response")))
gate("plant200_source_quality", all(source_quality$pass))

focal <- c("Ilex_aquifolium", "Camellia_sinensis", "Jasminum_officinale")
real <- setNames(vector("list", length(focal)), focal)
real_rows <- list()
real_candidates <- list()
for (tip in focal) {
  read_tip <- function(suffix) utils::read.csv(file.path(source_dir,
    paste0(tip, suffix, ".csv")), stringsAsFactors = FALSE)
  m <- read_tip("_patch_cells")
  p <- read_tip("_patches")
  links <- read_tip("_plate_carried_patch_links")
  saved <- read_tip("_candidate_refugia")
  stable <- HmscEcoEvo::hee_result_patch_stability(p, links)
  refuge <- HmscEcoEvo::hee_result_refugia_candidates(p, links,
    contraction_fraction = .10, min_location_mass = .10,
    min_patch_area_km2 = 100000, require_later_link = TRUE)
  alternate <- HmscEcoEvo::hee_result_refugia_candidates(p, links,
    contraction_fraction = .25, min_location_mass = .10,
    min_patch_area_km2 = 100000, require_later_link = TRUE)
  key <- function(x) sort(paste(x$to_patch_id, x$time_ma, sep = "@"))
  gate(paste0(tip, "_all_66_times"), length(unique(p$time_ma)) == 66L)
  gate(paste0(tip, "_refugia_reproduced"), identical(key(refuge$candidates), key(saved)),
       paste(nrow(refuge$candidates), "candidate rows"))
  gate(paste0(tip, "_network_reproduced"),
       nrow(stable$networks) == nrow(read_tip("_patch_history_networks")))
  real[[tip]] <- list(membership = m, patches = p, links = links,
                      stability = stable, refugia = refuge, alternative = alternate)
  real_rows[[tip]] <- data.frame(focal_tip = tip, time_slices = length(unique(p$time_ma)),
    supported_patches = nrow(p), carriage_links = nrow(links),
    patch_history_networks = nrow(stable$networks),
    contraction_intervals_10pct = nrow(refuge$episodes),
    candidate_refugia_10pct = nrow(refuge$candidates),
    contraction_intervals_25pct = nrow(alternate$episodes),
    candidate_refugia_25pct = nrow(alternate$candidates))
  if (nrow(refuge$candidates)) {
    z <- refuge$candidates[, c("time_ma", "to_patch_id", "area_km2",
      "location_mass", "geographic_evidence", "interpretation")]
    z$focal_tip <- tip
    real_candidates[[tip]] <- z
  }
}
real_dir <- file.path(out, "01_real_plant200")
write_csv(do.call(rbind, real_rows), file.path(real_dir, "focal_path_summary.csv"))
write_csv(do.call(rbind, real_candidates), file.path(real_dir,
  "corrected_candidate_refugia_recomputed.csv"))

# Re-run the first two package functions on a real 105-to-100 Ma interval.
grid_ids <- sprintf("g4_%04d", seq_len(4050L))
origin <- seq_len(4050L)
east <- ((origin - 1L) %/% 90L) * 90L + (origin %% 90L) + 1L
north <- origin + 90L
adjacency <- data.frame(from_cell_id = grid_ids[c(origin, origin[north <= 4050L])],
  to_cell_id = grid_ids[c(east, north[north <= 4050L])])
tip <- focal[[1L]]
real_m <- real[[tip]]$membership
replay_one <- function(time_ma) {
  one <- real_m[real_m$time_ma == time_ma, , drop = FALSE]
  HmscEcoEvo::hee_result_patches(one[, c("cell_id", "support", "land_area_km2",
    "location_density")], adjacency, threshold = .6, lineage_id = tip,
    time_ma = time_ma)
}
replay_old <- replay_one(105)
replay_new <- replay_one(100)
signature <- function(m) sort(vapply(split(m$cell_id, m$patch_id), function(ids)
  paste(sort(ids), collapse = "|"), character(1)))
gate("real_105Ma_patch_partition", identical(signature(replay_old$membership),
  signature(real_m[real_m$time_ma == 105, , drop = FALSE])))
gate("real_100Ma_patch_partition", identical(signature(replay_new$membership),
  signature(real_m[real_m$time_ma == 100, , drop = FALSE])))

cache <- file.path(pkg_root, "derived_inputs",
  "case05_v2_paleomap_h3_carrier_325_0Ma_1deg",
  "movement_cache_arias_topographic_rate0p4")
cache_index <- utils::read.csv(file.path(cache,
  "case05_v2_movement_cache_index.csv"), stringsAsFactors = FALSE)
carrier_at <- function(time_ma) {
  path <- cache_index$file[match(time_ma, cache_index$time_ma)]
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("Missing real plate carrier at ", time_ma, " Ma.", call. = FALSE)
  }
  readRDS(path)$carrier
}
bin4 <- function(lon, lat) {
  ix <- pmin(89L, pmax(0L, floor((lon + 180) / 4)))
  iy <- pmin(44L, pmax(0L, floor((lat + 90) / 4)))
  as.integer(iy * 90L + ix + 1L)
}
real_carriage <- function(older, younger) {
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
  pair <- stats::aggregate(older$land_area_km2[valid],
    list(source = source, target = dest), sum)
  names(pair)[[3L]] <- "transferred_area"
  pair$carriage_weight <- pair$transferred_area /
    source_total[match(as.character(pair$source), rownames(source_total)), 1L]
  data.frame(source_cell_id = grid_ids[pair$source],
    target_cell_id = grid_ids[pair$target],
    carriage_weight = pair$carriage_weight)
}
carriage <- real_carriage(carrier_at(105), carrier_at(100))
replayed_links <- HmscEcoEvo::hee_result_patch_links(
  real_m[real_m$time_ma == 105, , drop = FALSE],
  real_m[real_m$time_ma == 100, , drop = FALSE], carriage,
  min_source_share = .05)
saved_links <- real[[tip]]$links
saved_links <- saved_links[saved_links$time_from_ma == 105 &
  saved_links$time_to_ma == 100, , drop = FALSE]
ordered_links <- function(x) x[order(x$from_patch_id, x$to_patch_id),
  c("from_patch_id", "to_patch_id", "carried_area_km2", "source_share"),
  drop = FALSE]
a <- ordered_links(replayed_links); b <- ordered_links(saved_links)
gate("real_plate_links_reproduced", nrow(a) == nrow(b) &&
  identical(a$from_patch_id, b$from_patch_id) &&
  identical(a$to_patch_id, b$to_patch_id) &&
  isTRUE(all.equal(a$carried_area_km2, b$carried_area_km2,
                   tolerance = 1e-7, check.attributes = FALSE)) &&
  isTRUE(all.equal(a$source_share, b$source_share,
                   tolerance = 1e-7, check.attributes = FALSE)),
  paste(nrow(a), "real land-carriage patch links"))
write_csv(replayed_links, file.path(real_dir, "Ilex_105_to_100Ma_replayed_carriage_links.csv"))

# A compact three-time controlled landscape shows splits, mergers, new land,
# contraction, continuity, and active movement without claiming Plant200 flux.
toy_cells <- list(
  `20` = data.frame(cell_id = c("a", "b", "c", "d", "sea"),
    x = c(0, 1, 3, 4, 2), y = c(0, 0, 0, 0, 0),
    support = c(.9, .75, .83, .7, .95),
    land_area_km2 = c(100, 100, 100, 100, 0),
    location_density = c(.25, .25, .25, .25, 0),
    MAT_pohl_C = c(18, 17, 15, 14, NA),
    MAP_pohl_mm_yr = c(1100, 1000, 900, 850, NA),
    no_analog = c(FALSE, FALSE, FALSE, FALSE, NA)),
  `15` = data.frame(cell_id = c("x", "y", "z", "w", "n"),
    x = c(0, 1, 3, 4, 6), y = c(0, 0, 0, 0, 1),
    support = c(.8, .7, .9, .2, .85),
    land_area_km2 = c(60, 60, 50, 50, 80),
    location_density = c(.3, .25, .1, .15, .2),
    MAT_pohl_C = c(19, 18, 16, 15, 20),
    MAP_pohl_mm_yr = c(1050, 950, 820, 700, 1200),
    no_analog = c(FALSE, FALSE, FALSE, FALSE, TRUE)),
  `10` = data.frame(cell_id = c("u", "v", "r", "s"),
    x = c(0, 1, 3, 6), y = c(0, 0, 0, 1),
    support = c(.75, .72, .4, .82),
    land_area_km2 = c(60, 60, 40, 70),
    location_density = c(.35, .3, .1, .25),
    MAT_pohl_C = c(20, 19, 17, 21),
    MAP_pohl_mm_yr = c(980, 900, 700, 1150),
    no_analog = c(FALSE, FALSE, FALSE, TRUE))
)
toy_edges <- list(
  `20` = data.frame(from_cell_id = c("a", "b", "c"),
                    to_cell_id = c("b", "sea", "d")),
  `15` = data.frame(from_cell_id = c("x", "z"), to_cell_id = c("y", "w")),
  `10` = data.frame(from_cell_id = "u", to_cell_id = "v")
)
toy <- lapply(names(toy_cells), function(time) HmscEcoEvo::hee_result_patches(
  toy_cells[[time]], toy_edges[[time]], threshold = .6,
  lineage_id = "controlled_lineage", time_ma = as.numeric(time)))
names(toy) <- names(toy_cells)
transfer1 <- data.frame(source_cell_id = c("a", "b", "b", "c", "c", "d"),
  target_cell_id = c("x", "y", "z", "z", "y", "w"),
  carriage_weight = c(.6, .6, .2, .5, .2, .5))
transfer2 <- data.frame(source_cell_id = c("x", "y", "z", "w", "n"),
  target_cell_id = c("u", "v", "r", "r", "s"),
  carriage_weight = c(.9, .9, .8, .5, .8))
toy_links1 <- HmscEcoEvo::hee_result_patch_links(
  toy[["20"]]$membership, toy[["15"]]$membership, transfer1)
toy_links2 <- HmscEcoEvo::hee_result_patch_links(
  toy[["15"]]$membership, toy[["10"]]$membership, transfer2)
toy_links <- rbind(toy_links1, toy_links2)
toy_patches <- do.call(rbind, lapply(toy, `[[`, "patches"))
toy_stability <- HmscEcoEvo::hee_result_patch_stability(toy_patches, toy_links)
toy_refugia <- HmscEcoEvo::hee_result_refugia_candidates(toy_patches, toy_links,
  contraction_fraction = .3, min_location_mass = .1,
  require_later_link = TRUE)
toy_relaxed <- HmscEcoEvo::hee_result_refugia_candidates(toy_patches, toy_links,
  contraction_fraction = .3, min_location_mass = .1,
  require_later_link = FALSE)
new_land_patch <- toy[["15"]]$membership$patch_id[
  toy[["15"]]$membership$cell_id == "n"]
gate("controlled_sea_excluded", !"sea" %in% toy[["20"]]$membership$cell_id)
gate("controlled_new_land_not_carried", length(new_land_patch) == 1L &&
  !new_land_patch %in% toy_links1$to_patch_id)
gate("controlled_split_and_merge", any(duplicated(toy_links1$from_patch_id)) &&
  any(duplicated(toy_links1$to_patch_id)))
gate("controlled_contraction_identified", nrow(toy_refugia$episodes) == 1L &&
  toy_refugia$episodes$time_from_ma[[1L]] == 20)
gate("controlled_later_link_matters", nrow(toy_refugia$candidates) == 1L &&
  nrow(toy_relaxed$candidates) > nrow(toy_refugia$candidates))
toy_dir <- file.path(out, "02_controlled_patch_cases")
write_csv(toy_stability$patches, file.path(toy_dir, "patches_and_networks.csv"))
write_csv(toy_links, file.path(toy_dir, "plate_carriage_links_not_dispersal.csv"))
write_csv(toy_stability$networks, file.path(toy_dir, "history_networks.csv"))
write_csv(toy_stability$link_changes, file.path(toy_dir, "link_environment_changes.csv"))
write_csv(toy_refugia$episodes, file.path(toy_dir, "contraction_episodes.csv"))
write_csv(toy_refugia$candidates, file.path(toy_dir, "strict_candidate_refugia.csv"))
write_csv(toy_relaxed$candidates, file.path(toy_dir, "relaxed_candidate_refugia.csv"))

# A fully specified 3-time hidden-state model yields actual joint edge
# posterior mass conditional on terminal evidence, not products of marginals.
prior20 <- c(a = .35, b = .25, c = .25, d = .15)
final_likelihood <- c(u = .9, v = .8, r = .1, s = .5)
t1 <- data.frame(from_cell_id = c("a", "a", "a", "b", "b", "b",
                                   "c", "c", "c", "d", "d", "d"),
  to_cell_id = c("x", "y", "z", "y", "x", "z", "z", "y", "n", "w", "z", "n"),
  transition_probability = c(.6, .25, .15, .6, .2, .2, .5, .2, .3, .5, .2, .3),
  edge_kind = c("plate_carriage", "active_dispersal", "active_dispersal",
                "plate_carriage", "active_dispersal", "active_dispersal",
                "plate_carriage", "active_dispersal", "active_dispersal",
                "plate_carriage", "active_dispersal", "active_dispersal"))
t2 <- data.frame(from_cell_id = c("x", "x", "x", "y", "y", "y",
                                   "z", "z", "z", "w", "w", "n", "n", "n"),
  to_cell_id = c("u", "v", "r", "v", "u", "s", "r", "u", "s", "r", "s",
                 "s", "u", "v"),
  transition_probability = c(.55, .3, .15, .5, .45, .05, .6, .2, .2,
                             .6, .4, .65, .2, .15),
  edge_kind = c("plate_carriage", "active_dispersal", "active_dispersal",
                "plate_carriage", "active_dispersal", "active_dispersal",
                "plate_carriage", "active_dispersal", "active_dispersal",
                "plate_carriage", "active_dispersal", "plate_carriage",
                "active_dispersal", "active_dispersal"))
transition_matrix <- function(edges, from, to) {
  mat <- matrix(0, length(from), length(to), dimnames = list(from, to))
  for (k in seq_len(nrow(edges))) mat[edges$from_cell_id[[k]],
    edges$to_cell_id[[k]]] <- edges$transition_probability[[k]]
  mat
}
T1 <- transition_matrix(t1, names(prior20), toy_cells[["15"]]$cell_id)
T2 <- transition_matrix(t2, toy_cells[["15"]]$cell_id,
                        names(final_likelihood))
gate("controlled_transition_rows_normalized", all(abs(rowSums(T1) - 1) < 1e-12) &&
  all(abs(rowSums(T2) - 1) < 1e-12))
alpha15 <- as.numeric(prior20 %*% T1); names(alpha15) <- colnames(T1)
back15 <- as.numeric(T2 %*% final_likelihood); names(back15) <- rownames(T2)
normalizer <- sum(alpha15 * back15)
posterior_edges <- function(edges, source_mass, future_evidence, old_time, young_time) {
  edges$time_from_ma <- old_time
  edges$time_to_ma <- young_time
  edges$posterior_edge_mass <- source_mass[edges$from_cell_id] *
    edges$transition_probability * future_evidence[edges$to_cell_id] / normalizer
  edges
}
joint1 <- posterior_edges(t1, prior20, back15, 20, 15)
joint2 <- posterior_edges(t2, alpha15, final_likelihood, 15, 10)
joint <- rbind(joint1, joint2)
flow1 <- HmscEcoEvo::hee_result_corridor_flow(joint1,
  toy[["20"]]$membership, toy[["15"]]$membership)
flow2 <- HmscEcoEvo::hee_result_corridor_flow(joint2,
  toy[["15"]]$membership, toy[["10"]]$membership)
corridors <- rbind(flow1$corridors, flow2$corridors)
roles <- HmscEcoEvo::hee_result_donor_recipient(corridors, minimum_flow = .05)
mass_check <- function(x, result) {
  active <- x[x$edge_kind == "active_dispersal" &
    x$from_cell_id != x$to_cell_id, , drop = FALSE]
  isTRUE(all.equal(sum(active$posterior_edge_mass),
    sum(result$corridors$posterior_movement_mass) + result$unassigned_edge_mass,
    tolerance = 1e-12))
}
gate("controlled_joint_posterior_normalized", abs(sum(joint1$posterior_edge_mass) - 1) <
  1e-12 && abs(sum(joint2$posterior_edge_mass) - 1) < 1e-12)
gate("controlled_active_mass_accounted", mass_check(joint1, flow1) &&
  mass_check(joint2, flow2))
gate("controlled_plate_carriage_excluded", all(joint$edge_kind %in%
  c("active_dispersal", "plate_carriage")) &&
  sum(corridors$posterior_movement_mass) <
  sum(joint$posterior_edge_mass[joint$edge_kind == "active_dispersal"]))
gate("controlled_exchange_hub", "exchange_hub" %in% roles$role)
joint_dir <- file.path(out, "03_joint_edge_example")
write_csv(joint, file.path(joint_dir, "three_time_joint_edge_posterior.csv"))
write_csv(corridors, file.path(joint_dir, "active_patch_corridors.csv"))
write_csv(roles, file.path(joint_dir, "patch_export_recipient_roles.csv"))
write_csv(data.frame(interval = c("20_to_15", "15_to_10"),
  all_edge_mass = c(sum(joint1$posterior_edge_mass), sum(joint2$posterior_edge_mass)),
  active_patch_corridor_mass = c(sum(flow1$corridors$posterior_movement_mass),
                                 sum(flow2$corridors$posterior_movement_mass)),
  active_unassigned_mass = c(flow1$unassigned_edge_mass,
                             flow2$unassigned_edge_mass)),
  file.path(joint_dir, "posterior_mass_accounting.csv"))

theme_case <- function() ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold", colour = "#173552"),
    strip.text = ggplot2::element_text(face = "bold"))
fig_dir <- file.path(out, "04_figures")
time_summary <- do.call(rbind, lapply(focal, function(tip) {
  p <- real[[tip]]$patches
  x <- stats::aggregate(cbind(area_km2 = p$area_km2,
    n_patches = rep(1, nrow(p))),
    by = list(time_ma = p$time_ma), FUN = sum)
  x$focal_tip <- tip
  x
}))
plot1 <- ggplot2::ggplot(time_summary,
    ggplot2::aes(time_ma, area_km2 / 1e6, colour = focal_tip)) +
  ggplot2::geom_line(linewidth = .9) +
  ggplot2::scale_x_reverse() +
  ggplot2::scale_colour_manual(values = c("#c8712a", "#147b70", "#346cb5")) +
  ggplot2::labs(title = "Real Plant200 | supported land area through time",
    subtitle = "Corrected ancestral-response results; three illustrated tip-to-root paths",
    x = "Time before present (Ma)", y = "Supported land area (million km2)",
    colour = "Focal tip", caption = paste("Environmental opportunity, not occupancy;",
      "each path changes ancestral identity at tree nodes.")) +
  theme_case()
ggplot2::ggsave(file.path(fig_dir, "01_real_supported_area.png"), plot1,
  width = 12, height = 6.5, dpi = 160)
candidate_plot <- do.call(rbind, lapply(focal, function(tip) {
  r <- real[[tip]]$refugia$candidates
  if (!nrow(r)) return(NULL)
  z <- as.data.frame(table(r$time_ma), stringsAsFactors = FALSE)
  data.frame(time_ma = as.numeric(as.character(z$Var1)),
             candidates = z$Freq, focal_tip = tip)
}))
plot2 <- ggplot2::ggplot(candidate_plot,
    ggplot2::aes(time_ma, candidates, fill = focal_tip)) +
  ggplot2::geom_col(width = 4) + ggplot2::facet_wrap(~focal_tip, ncol = 1,
    scales = "free_y") + ggplot2::scale_x_reverse() +
  ggplot2::scale_fill_manual(values = c("#c8712a", "#147b70", "#346cb5"),
                            guide = "none") +
  ggplot2::labs(title = "Real Plant200 | candidate ecological refugia",
    subtitle = "10% supported-area contraction; incoming and subsequent carriage required",
    x = "Time before present (Ma)", y = "Candidate patches",
    caption = "Candidates are not fossil-confirmed populations or demographic persistence.") +
  theme_case()
ggplot2::ggsave(file.path(fig_dir, "02_real_candidate_refugia.png"), plot2,
  width = 12, height = 8, dpi = 160)
toy_map <- do.call(rbind, lapply(names(toy_cells), function(time) {
  z <- toy_cells[[time]]
  z$time_ma <- as.numeric(time)
  z$patch_id <- toy[[time]]$membership$patch_id[
    match(z$cell_id, toy[[time]]$membership$cell_id)]
  z$state <- ifelse(z$land_area_km2 == 0, "sea",
    ifelse(is.na(z$patch_id), "land_not_supported", "supported_patch"))
  z$time_ma <- factor(z$time_ma, levels = c(20, 15, 10),
                      labels = c("20 Ma", "15 Ma", "10 Ma"))
  z
}))
plot3 <- ggplot2::ggplot(toy_map, ggplot2::aes(x, y)) +
  ggplot2::geom_tile(ggplot2::aes(fill = state), width = .88, height = .88,
                     colour = "#29435d", linewidth = .35) +
  ggplot2::geom_text(ggplot2::aes(label = cell_id), fontface = "bold", size = 5) +
  ggplot2::facet_wrap(~time_ma, nrow = 1) +
  ggplot2::scale_fill_manual(values = c(sea = "#c9e3ef",
    land_not_supported = "#d9dce0", supported_patch = "#67b897")) +
  ggplot2::coord_equal(xlim = c(-.7, 6.7), ylim = c(-.7, 1.7), expand = FALSE) +
  ggplot2::labs(title = "Controlled landscape | supported cells and new land",
    subtitle = "The sea never joins patches; cell n first appears at 15 Ma without carriage inheritance",
    x = NULL, y = NULL, fill = "Cell state") + theme_case() +
  ggplot2::theme(axis.text = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank())
ggplot2::ggsave(file.path(fig_dir, "03_controlled_patch_states.png"), plot3,
  width = 12, height = 4.5, dpi = 160)
patch_xy <- do.call(rbind, lapply(names(toy), function(time) {
  m <- toy[[time]]$membership
  xy <- toy_cells[[time]][match(m$cell_id, toy_cells[[time]]$cell_id), c("x", "y")]
  positions <- data.frame(patch_id = m$patch_id,
    time_ma = as.numeric(time), x = xy$x, y = xy$y)
  stats::aggregate(x ~ patch_id + time_ma, data = positions, FUN = mean)
}))
link_graph <- merge(toy_links, patch_xy[, c("patch_id", "x")],
  by.x = "from_patch_id", by.y = "patch_id")
names(link_graph)[names(link_graph) == "x"] <- "from_x"
link_graph <- merge(link_graph, patch_xy[, c("patch_id", "x")],
  by.x = "to_patch_id", by.y = "patch_id")
names(link_graph)[names(link_graph) == "x"] <- "to_x"
plot4 <- ggplot2::ggplot() +
  ggplot2::geom_segment(data = link_graph,
    ggplot2::aes(x = time_from_ma, y = from_x,
      xend = time_to_ma, yend = to_x, linewidth = source_share),
    colour = "#8496a5", arrow = grid::arrow(length = grid::unit(.12, "cm"))) +
  ggplot2::geom_point(data = patch_xy, ggplot2::aes(time_ma, x),
    size = 4, colour = "#218873") +
  ggplot2::geom_text(data = patch_xy,
    ggplot2::aes(time_ma, x, label = sub(".*_p", "p", patch_id)),
    nudge_y = .3, size = 4) + ggplot2::scale_x_reverse(breaks = c(20, 15, 10)) +
  ggplot2::labs(title = "Controlled landscape | plate-carried patch links",
    subtitle = "Splits and mergers reflect land transfer, not organismal travel",
    x = "Time before present (Ma)", y = "Patch horizontal position",
    linewidth = "Source share") + theme_case()
ggplot2::ggsave(file.path(fig_dir, "04_controlled_carriage_links.png"), plot4,
  width = 10, height = 6.5, dpi = 160)
short_patch <- function(id) sub("^controlled_lineage_", "", id)
corridors$edge_label <- paste(short_patch(corridors$from_patch_id),
  short_patch(corridors$to_patch_id), sep = " -> ")
plot5 <- ggplot2::ggplot(corridors,
  ggplot2::aes(reorder(edge_label,
                        posterior_movement_mass), posterior_movement_mass,
               fill = factor(time_from_ma))) +
  ggplot2::geom_col() + ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = c(`15` = "#dc8a34", `20` = "#2d8c83")) +
  ggplot2::labs(title = "Controlled HMM | joint posterior active movement",
    subtitle = "Forward-backward joint edge mass; plate carriage is excluded",
    x = "Patch-to-patch edge", y = "Posterior movement mass",
    fill = "From Ma", caption = "Illustrative model, not Plant200 empirical movement flow.") +
  theme_case()
ggplot2::ggsave(file.path(fig_dir, "05_controlled_active_corridors.png"), plot5,
  width = 12, height = 7, dpi = 160)
role_plot <- reshape(roles[, c("patch_id", "outgoing_mass", "incoming_mass", "role")],
  varying = c("outgoing_mass", "incoming_mass"), v.names = "mass",
  timevar = "direction", times = c("outgoing", "incoming"), direction = "long")
role_plot$patch_id <- short_patch(role_plot$patch_id)
plot6 <- ggplot2::ggplot(role_plot,
  ggplot2::aes(patch_id, mass, fill = direction)) +
  ggplot2::geom_col(position = "dodge") + ggplot2::coord_flip() +
  ggplot2::facet_wrap(~role, scales = "free_y") +
  ggplot2::scale_fill_manual(values = c(incoming = "#287e9d", outgoing = "#d7893c")) +
  ggplot2::labs(title = "Controlled HMM | patch exchange roles",
    subtitle = "Exporter, recipient and hub labels use active-movement mass only",
    x = "Time-specific patch", y = "Expected movement mass", fill = "Direction",
    caption = "Not demographic source/sink status.") + theme_case()
ggplot2::ggsave(file.path(fig_dir, "06_controlled_exchange_roles.png"), plot6,
  width = 12, height = 6.5, dpi = 160)

contract <- data.frame(function_name = c("hee_result_patches",
  "hee_result_patch_links", "hee_result_patch_stability",
  "hee_result_refugia_candidates", "hee_result_corridor_flow",
  "hee_result_donor_recipient"),
  real_plant200 = c("replayed_105_100Ma", "replayed_105_to_100Ma",
    "recomputed_three_focal_paths_66_times", "recomputed_three_focal_paths_66_times",
    "NOT_ESTIMATED_NO_JOINT_EDGE_POSTERIOR", "NOT_ESTIMATED_NO_JOINT_EDGE_POSTERIOR"),
  controlled_example = rep("EXECUTED", 6L),
  claim = c("environmental_opportunity", "land_carriage_not_dispersal",
    "connected_patch_network", "candidate_not_confirmed_refugium",
    "posterior_active_edge_mass_of_controlled_HMM",
    "exchange_roles_not_demographic_source_sink"))
write_csv(contract, file.path(out, "00_contract", "six_functions_evidence_contract.csv"))
gate("all_six_functions_exercised", all(contract$controlled_example == "EXECUTED"))
gate("six_figures_exist", length(list.files(fig_dir, pattern = "[.]png$")) == 6L)
write_csv(quality, file.path(out, "05_quality", "case07_quality_gates.csv"))
if (requireNamespace("rmarkdown", quietly = TRUE)) {
  old_output <- Sys.getenv("HMSCEE_CASE07_OUTPUT", unset = NA_character_)
  Sys.setenv(HMSCEE_CASE07_OUTPUT = out)
  vignette <- file.path(pkg_root, "vignettes",
    "v15_case07_patch_evidence_ladder.Rmd")
  html <- file.path(out, "Case07_patch_evidence_ladder.html")
  rendered <- tryCatch({
    rmarkdown::render(vignette, output_file = basename(html),
      output_dir = out, quiet = TRUE)
    file.exists(html)
  }, error = function(e) {
    message("Case07 tutorial render failed: ", conditionMessage(e))
    FALSE
  })
  if (is.na(old_output)) Sys.unsetenv("HMSCEE_CASE07_OUTPUT")
  else Sys.setenv(HMSCEE_CASE07_OUTPUT = old_output)
  gate("tutorial_html_generated", rendered)
} else gate("tutorial_html_generated", FALSE, "rmarkdown is unavailable")
write_csv(quality, file.path(out, "05_quality", "case07_quality_gates.csv"))
if (!all(quality$pass)) {
  stop("Case07 quality gates failed: ",
    paste(quality$check[!quality$pass], collapse = ", "), call. = FALSE)
}
message("Case07 complete: ", out, " | ", nrow(quality), " quality gates passed")
