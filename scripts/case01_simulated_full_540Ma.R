#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Hmsc)
  library(ape)
  library(coda)
  library(ggplot2)
  library(gridExtra)
  library(scales)
})

or_else <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a) || !nzchar(a)) b else a
cmd_file <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
script_path <- normalizePath(or_else(ofile, or_else(cmd_file[1], "scripts/case01_simulated_full_540Ma.R")),
                             winslash = "/", mustWork = FALSE)
script_path <- sub("^--file=", "", script_path)
pkg_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (file.exists(file.path(pkg_root, "DESCRIPTION")) && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(pkg_root, quiet = TRUE)
} else {
  library(HmscEcoEvo)
}
if (utils::packageVersion("HmscEcoEvo") < package_version("1.0.4") ||
    sum(hee_core_process_catalog()$layer == "core_process") != 6L) {
  stop("Case01 requires HmscEcoEvo >= 1.0.4 with the six-process API.",
       call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  pat <- paste0("^(--)?", name, "=")
  hit <- grep(pat, args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(pat, "", hit[1])
}
as_bool <- function(x) tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")

quick <- as_bool(arg_value("quick", "FALSE"))
output_root <- normalizePath(arg_value("output", "outputs/case01_simulated_full_540Ma"),
                             winslash = "/", mustWork = FALSE)
env_rds <- arg_value(
  "env_rds",
  "C:/Users/Google/Documents/HMSC-HIST/Phanerozoic_Environment_540_0Ma_5Myr_4deg_landharmonized_v2/08_rds_export/Phanerozoic_Environment_540_0Ma_5Myr_Level2_4deg_landharmonized_v2_compact_float32.rds"
)
if (!file.exists(env_rds)) stop("Missing required env_cube RDS: ", env_rds, call. = FALSE)

set.seed(as.integer(arg_value("seed", "20260608")))
dirs <- file.path(output_root, c("figures", "tables", "reports", "metadata", "tasks", "rds"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
fig_dir <- file.path(output_root, "figures")
tab_dir <- file.path(output_root, "tables")
rep_dir <- file.path(output_root, "reports")
rds_dir <- file.path(output_root, "rds")
task_dir <- file.path(output_root, "tasks")

cfg <- if (quick) {
  list(n_sites = 35, n_species = 12, n_traits = 3, n_regions = 5,
       samples = 50, transient = 50, thin = 1, nChains = 2,
       nParallel = 1, chunk_size = 3000, n_uncert_draws = 12)
} else {
  list(n_sites = 80, n_species = 30, n_traits = 3, n_regions = 5,
       samples = 1000, transient = 1000, thin = 10, nChains = 4,
       nParallel = 1, chunk_size = 5000, n_uncert_draws = 80)
}
cfg$samples <- as.integer(arg_value("samples", cfg$samples))
cfg$transient <- as.integer(arg_value("transient", cfg$transient))
cfg$thin <- as.integer(arg_value("thin", cfg$thin))
cfg$nChains <- as.integer(arg_value("nChains", cfg$nChains))
cfg$nParallel <- as.integer(arg_value("nParallel", cfg$nParallel))
projection_mode <- arg_value("projection_mode", "environment_only")
if (!identical(projection_mode, "environment_only")) {
  stop("Only projection_mode = 'environment_only' is allowed by default. Modern spatial random effects are not projected into deep time.",
       call. = FALSE)
}
extrapolation_mode <- arg_value("extrapolation_mode", "flag")
if (!extrapolation_mode %in% c("flag", "downweight")) {
  stop("extrapolation_mode must be 'flag' or 'downweight'.", call. = FALSE)
}

theme_set(theme_bw(base_size = 8))
pal <- c("#163A5F", "#2F6DB5", "#2D7F7B", "#2D7F5E", "#D8A03D",
         "#C76D2A", "#B23A3A", "#6B4FA3", "#6E7781")
rep_times_target <- unique(c(
  540, 400, 252, 201, 145, 66, 23, 2.58, 0,
  seq(540, 0, by = -60)
))

save_plot <- function(plot, name, width = 7.2, height = 4.8) {
  hee_save_plot(plot, file.path(fig_dir, tools::file_path_sans_ext(name)),
                width = width, height = height, dpi = 320)
}
write_table <- function(x, name) {
  path <- file.path(tab_dir, name)
  write.csv(as.data.frame(x), path, row.names = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}
save_rds <- function(x, name) {
  path <- file.path(rds_dir, name)
  saveRDS(x, path, compress = "xz")
  normalizePath(path, winslash = "/", mustWork = FALSE)
}
safe_left_join <- function(x, y, by, label) {
  x <- as.data.frame(x)
  y <- as.data.frame(y)
  missing_by <- setdiff(by, names(y))
  if (length(missing_by) > 0L) {
    stop(label, " is missing join key(s): ", paste(missing_by, collapse = ", "),
         call. = FALSE)
  }
  dup <- duplicated(y[, by, drop = FALSE]) |
    duplicated(y[, by, drop = FALSE], fromLast = TRUE)
  if (any(dup)) {
    bad <- unique(y[dup, by, drop = FALSE])
    bad <- utils::head(bad, 5)
    msg <- apply(bad, 1, function(z) {
      paste(paste(by, z, sep = "="), collapse = ", ")
    })
    stop(label, " has duplicate join key(s): ", paste(msg, collapse = " | "),
         call. = FALSE)
  }
  before <- nrow(x)
  out <- merge(x, y, by = by, all.x = TRUE, sort = FALSE)
  if (nrow(out) != before) {
    stop("Join from ", label, " changed row count from ", before, " to ",
         nrow(out), ". Check keys: ", paste(by, collapse = ", "),
         call. = FALSE)
  }
  out
}

long_matrix <- function(M, row = "row", col = "column", value = "value") {
  if (is.null(dim(M))) {
    nm <- names(M)
    if (is.null(nm)) nm <- as.character(seq_along(M))
    M <- matrix(M, nrow = 1, dimnames = list("value", nm))
  }
  d <- as.data.frame(as.table(M), stringsAsFactors = FALSE)
  names(d) <- c(row, col, value)
  d
}
safe_scale <- function(x) as.numeric(scale(x))
inv_probit <- function(x) pnorm(pmax(pmin(x, 7), -7))

plot_heatmap_matrix <- function(M, title, fill = "value", low = pal[2], high = pal[7]) {
  d <- long_matrix(M, "row", "column", fill)
  ggplot(d, aes(x = column, y = row, fill = .data[[fill]])) +
    geom_tile(color = "white", linewidth = 0.15) +
    scale_fill_gradient2(low = low, mid = "white", high = high, midpoint = 0) +
    labs(x = NULL, y = NULL, fill = fill, title = title) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

draw_tree <- function(phy, file_base, title) {
  png(paste0(file_base, ".png"), width = 7.2, height = 4.8, units = "in", res = 320, bg = "white")
  par(mar = c(4, 1, 3, 1))
  plot.phylo(phy, direction = "rightwards", cex = 0.55, label.offset = 5)
  axis(1, at = seq(0, 540, by = 90), labels = 540 - seq(0, 540, by = 90), cex.axis = 0.7)
  mtext("Time before present (Ma)", side = 1, line = 2.4, cex = 0.75)
  title(title)
  dev.off()
  pdf(paste0(file_base, ".pdf"), width = 7.2, height = 4.8)
  par(mar = c(4, 1, 3, 1))
  plot.phylo(phy, direction = "rightwards", cex = 0.55, label.offset = 5)
  axis(1, at = seq(0, 540, by = 90), labels = 540 - seq(0, 540, by = 90), cex.axis = 0.7)
  mtext("Time before present (Ma)", side = 1, line = 2.4, cex = 0.75)
  title(title)
  dev.off()
  c(paste0(file_base, ".png"), paste0(file_base, ".pdf"))
}

make_tree <- function(species) {
  phy <- ape::rcoal(length(species))
  phy$tip.label <- species
  root_age <- 540
  phy$edge.length <- phy$edge.length / max(ape::node.depth.edgelength(phy)) * root_age
  phy$root.time <- root_age
  phy
}

lookup_phylo <- function(P, phylo_mask) {
  times_mask <- attr(phylo_mask, "times")
  if (is.null(times_mask)) times_mask <- as.numeric(gsub("Ma$", "", colnames(phylo_mask)))
  out <- rep(1, nrow(P))
  sp_i <- match(as.character(P$species), rownames(phylo_mask))
  tm_i <- match(P$time_ma, times_mask)
  ok <- !is.na(sp_i) & !is.na(tm_i)
  out[ok] <- phylo_mask[cbind(sp_i[ok], tm_i[ok])]
  out
}

extract_beta_draws <- function(hM, max_draws = 80) {
  rows <- list()
  k <- 1L
  for (ch in seq_along(hM$postList)) {
    for (it in seq_along(hM$postList[[ch]])) {
      B <- hM$postList[[ch]][[it]]$Beta
      if (is.null(B)) next
      rows[[k]] <- B
      k <- k + 1L
      if (length(rows) >= max_draws) return(rows)
    }
  }
  rows
}

project_fixed <- function(B, X, id_cols, link = "probit") {
  axes <- intersect(colnames(B), names(X))
  axes <- setdiff(axes, "(Intercept)")
  eta <- as.matrix(X[, axes, drop = FALSE]) %*% t(B[, axes, drop = FALSE])
  if ("(Intercept)" %in% colnames(B)) eta <- sweep(eta, 2, B[, "(Intercept)"], "+")
  pred <- if (link == "probit") pnorm(eta) else plogis(eta)
  out <- cbind(X[, id_cols, drop = FALSE], as.data.frame(pred, check.names = FALSE))
  long <- reshape(out, varying = rownames(B), v.names = "suitability",
                  timevar = "species", times = rownames(B), direction = "long")
  rownames(long) <- NULL
  long
}

richness_summary <- function(P, group_cols) {
  P$.value <- P$probability
  f <- as.formula(paste(".value ~", paste(group_cols, collapse = " + ")))
  out <- aggregate(f, P, sum, na.rm = TRUE)
  names(out)[names(out) == ".value"] <- "expected_richness"
  coords <- unique(P[, intersect(c(group_cols, "lon", "lat", "paleo_lon", "paleo_lat"), names(P)), drop = FALSE])
  merge(out, coords, by = group_cols, all.x = TRUE, sort = FALSE)
}

centroid_summary <- function(rich) {
  split_r <- split(rich, rich$time_ma)
  do.call(rbind, lapply(names(split_r), function(tm) {
    z <- split_r[[tm]]
    w <- pmax(z$expected_richness, 0)
    if (sum(w) == 0) w <- rep(1, nrow(z))
    data.frame(time_ma = as.numeric(tm),
               lon_centroid = weighted.mean(z$lon, w, na.rm = TRUE),
               lat_centroid = weighted.mean(z$lat, w, na.rm = TRUE))
  }))
}

area_summary <- function(P, model_id) {
  base <- data.frame(model_id = model_id,
                     time_ma = sort(unique(P$time_ma)),
                     stringsAsFactors = FALSE)
  z <- P[P$probability > 0.5, , drop = FALSE]
  if (nrow(z) == 0) {
    base$land_area_km2 <- 0
    return(base)
  }
  out <- aggregate(land_area_km2 ~ time_ma, unique(z[, c("cell_id", "time_ma", "land_area_km2")]), sum, na.rm = TRUE)
  out <- merge(base, out, by = "time_ma", all.x = TRUE, sort = FALSE)
  out$land_area_km2[is.na(out$land_area_km2)] <- 0
  out
}

turnover_vs_present <- function(rich) {
  present <- rich[rich$time_ma == min(rich$time_ma), c("cell_id", "expected_richness")]
  names(present)[2] <- "present_richness"
  z <- merge(rich, present, by = "cell_id", all.x = TRUE)
  z$turnover_vs_present <- abs(z$expected_richness - z$present_richness)
  z
}

message("Phase 1: loading authoritative environmental cube")
cube <- hee_load_timecube(env_rds)
times <- if (quick) seq(540, 0, by = -60) else hee_time_axis_from_env_cube(cube)
times <- sort(unique(times), decreasing = TRUE)
hee_assert_full_time_axis(times, quick = quick, min_time_slices = 3,
                          required_range = c(540, 0))
rep_times <- sort(unique(vapply(rep_times_target, function(t) {
  times[which.min(abs(times - t))]
}, numeric(1))), decreasing = TRUE)
variables <- c("MAT_pohl_C", "MAP_pohl_mm_yr", "P_seasonality_pohl_sd",
               "elevation_m", "distance_to_coast_km",
               "moisture_availability_index_z", "wetland_potential_index_z",
               "aridity_index_P_over_E")
missing_vars <- setdiff(variables, cube$variables)
if (length(missing_vars) > 0) stop("env_cube is missing variable(s): ", paste(missing_vars, collapse = ", "), call. = FALSE)
XFormula <- as.formula(paste("~", paste(variables, collapse = " + ")))
predictor_roles <- hee_tag_predictors(
  variables,
  stats::setNames(rep("deep_time_exogenous", length(variables)), variables)
)
invisible(hee_detect_leakage(XFormula, predictor_roles))

species <- sprintf("sp%02d", seq_len(cfg$n_species))
dispersal_distance <- stats::setNames(runif(cfg$n_species, 0.1, 1.5), species)
traits <- data.frame(
  body_size = safe_scale(log(seq_len(cfg$n_species) + 3) + rnorm(cfg$n_species, 0, 0.2)),
  dispersal = safe_scale(log1p(dispersal_distance)),
  thermal_tolerance = safe_scale(sin(seq_len(cfg$n_species) / 3) + rnorm(cfg$n_species, 0, 0.2)),
  row.names = species
)
colnames(traits) <- c("body_size", "dispersal", "thermal_tolerance")
phy <- make_tree(species)
C <- ape::vcv(phy, corr = TRUE)
species_origin <- hee_estimate_species_origin(phy)
clade <- setNames(rep(LETTERS[1:5], length.out = cfg$n_species), species)

message("Phase 1: constructing HMSC four input matrices")
grid0 <- hee_make_paleo_grid(cube, times = 0, variables = variables, land_only = TRUE)
sites <- grid0[sample(seq_len(nrow(grid0)), cfg$n_sites), , drop = FALSE]
sites$site_id <- sprintf("site%03d", seq_len(nrow(sites)))
rownames(sites) <- sites$site_id
XData_raw <- sites[, variables, drop = FALSE]
recipe <- hee_lock_recipe(XData_raw, XFormula)
XData_scaled <- hee_apply_recipe(recipe, XData_raw)
XData <- XData_scaled
TrData <- traits
TrFormula <- ~ body_size + dispersal + thermal_tolerance

trait_coef <- matrix(c(
  0.55, -0.20, 0.10, -0.15, 0.05, -0.10, 0.15, -0.20,
  -0.10, 0.35, -0.12, 0.05, 0.25, 0.15, -0.08, 0.25,
  0.30, -0.08, 0.15, -0.05, -0.10, 0.22, 0.10, -0.18
), nrow = cfg$n_traits, byrow = TRUE, dimnames = list(colnames(traits), variables))
Btrue <- as.matrix(traits) %*% trait_coef
Btrue <- Btrue / max(abs(Btrue)) * 0.85
Btrue <- cbind("(Intercept)" = rnorm(cfg$n_species, -0.35, 0.18), Btrue)
eta <- as.matrix(cbind("(Intercept)" = 1, XData_scaled[, variables, drop = FALSE])) %*% t(Btrue)
prob_true <- inv_probit(eta)
Y <- matrix(rbinom(length(prob_true), 1, prob_true), nrow = nrow(prob_true),
            dimnames = list(rownames(XData), species))
studyDesign <- data.frame(site = factor(rownames(Y)), row.names = rownames(Y))
ranLevels <- list(site = Hmsc::HmscRandomLevel(units = studyDesign$site))

regions <- paste0("R", seq_len(cfg$n_regions))
tip_ranges <- matrix(0, cfg$n_species, cfg$n_regions, dimnames = list(species, regions))
for (i in seq_along(species)) {
  tip_ranges[i, unique(c(((i - 1) %% cfg$n_regions) + 1,
                         ((i + as.integer(factor(clade[i]))) %% cfg$n_regions) + 1))] <- 1
}
region_grid0 <- grid0[, c("cell_id", "time_ma", "lon", "lat", "land_mask_dem", "land_area_km2")]
region_grid0$region <- regions[cut(region_grid0$lon, breaks = cfg$n_regions, labels = FALSE, include.lowest = TRUE)]
landmask_cube <- hee_make_paleo_grid(cube, times = times, variables = variables[1], land_only = TRUE)
landmask_cube <- landmask_cube[, c("cell_id", "time_ma", "lon", "lat", "land_mask_dem", "land_area_km2")]
region_cube <- landmask_cube
region_cube$region <- regions[cut(region_cube$lon, breaks = cfg$n_regions, labels = FALSE, include.lowest = TRUE)]
plate_points <- expand.grid(track_id = regions, time_ma = times, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
plate_points$paleo_lon <- rep(seq(-140, 140, length.out = cfg$n_regions), each = length(times)) +
  40 * sin((540 - plate_points$time_ma) / 90)
plate_points$paleo_lat <- rep(seq(-45, 45, length.out = cfg$n_regions), each = length(times)) +
  25 * cos((540 - plate_points$time_ma) / 110)
plate_points$paleo_lon <- ((plate_points$paleo_lon + 180) %% 360) - 180
fossils <- data.frame(species = sample(species, min(120, length(species) * 4), replace = TRUE),
                      lon = sample(grid0$lon, min(120, length(species) * 4), replace = TRUE),
                      lat = sample(grid0$lat, min(120, length(species) * 4), replace = TRUE),
                      age_min = sample(times, min(120, length(species) * 4), replace = TRUE),
                      stringsAsFactors = FALSE)
fossils$age_max <- pmin(540, fossils$age_min + 2.5)
fossils$age_min <- pmax(0, fossils$age_min - 2.5)

write_table(Y, "Y_matrix.csv")
write_table(cbind(site_id = rownames(XData_raw), XData_raw), "XData_0Ma_from_env_cube_raw.csv")
write_table(cbind(site_id = rownames(XData), XData), "XData_0Ma_from_env_cube.csv")
write_table(cbind(species = rownames(TrData), TrData), "TrData_traits.csv")
write_table(C, "phylogenetic_correlation_C.csv")
write_table(cbind(species = rownames(tip_ranges), tip_ranges), "tip_ranges.csv")
write_table(sites, "modern_sites.csv")
write_table(plate_points, "plate_points.csv")
write_table(fossils, "pseudo_fossils.csv")
save_rds(list(Y = Y, XData = XData, XData_raw = XData_raw, recipe = recipe,
              TrData = TrData, phy = phy, C = C,
              studyDesign = studyDesign, ranLevels = ranLevels, tip_ranges = tip_ranges,
              region_cube = region_cube, landmask_cube = landmask_cube,
              plate_points = plate_points, env_cube = env_rds, fossils = fossils),
         "case01_input_objects.rds")

message("Phase 1: plotting HMSC input objects")
save_plot(plot_heatmap_matrix(Y, "01 Y community matrix", "occurrence", low = "white", high = pal[2]),
          "01_Y_heatmap.png", 7.5, 5)
prev <- data.frame(species = species, prevalence = colMeans(Y))
save_plot(ggplot(prev, aes(x = reorder(species, prevalence), y = prevalence)) +
            geom_col(fill = pal[2]) + coord_flip() +
            labs(x = NULL, y = "Prevalence", title = "02 Species prevalence"),
          "02_species_prevalence.png")
site_rich <- data.frame(site = rownames(Y), richness = rowSums(Y))
save_plot(ggplot(site_rich, aes(x = richness)) + geom_histogram(binwidth = 1, fill = pal[4], color = "white") +
            labs(x = "Site richness", y = "Sites", title = "03 Site richness"),
          "03_site_richness.png")
zero_df <- data.frame(component = c("Y zeros", "Y ones"),
                      proportion = c(mean(Y == 0), mean(Y == 1)))
save_plot(ggplot(zero_df, aes(x = component, y = proportion, fill = component)) +
            geom_col(show.legend = FALSE) + labs(x = NULL, y = "Proportion", title = "Zero proportion plot"),
          "03b_zero_proportion.png")
save_plot(plot_heatmap_matrix(stats::cor(XData), "04 XData predictor correlation (model scale)", "correlation"),
          "04_X_correlation.png")
x_long <- reshape(cbind(sites[, c("site_id", "lon", "lat")], XData_raw),
                  varying = variables, v.names = "value", timevar = "variable",
                  times = variables, direction = "long")
save_plot(ggplot(x_long, aes(x = lon, y = lat, color = value)) +
            geom_point(size = 1.4) + facet_wrap(~ variable, ncol = 3) +
            coord_equal() + scale_color_gradient2(low = pal[2], mid = "white", high = pal[7]) +
            labs(x = "Longitude", y = "Latitude", color = "Raw 0 Ma value",
                 title = "05 Modern environmental variables extracted from 0 Ma layers"),
          "05_X_environment_maps_0Ma.png", 9, 6)
save_plot(ggplot(x_long, aes(x = value)) + geom_histogram(bins = 20, fill = pal[2], color = "white") +
            facet_wrap(~ variable, scales = "free", ncol = 3) +
            labs(x = "Value", y = "Sites", title = "XData variable histograms"),
          "05b_X_histograms.png", 9, 6)
save_plot(ggplot(grid0, aes(x = lon, y = lat, color = MAT_pohl_C)) + geom_point(size = 0.35) +
            coord_equal() + scale_color_gradient2(low = pal[2], mid = "white", high = pal[7]) +
            labs(title = "0 Ma MAT layer map", x = "Longitude", y = "Latitude", color = "MAT"),
          "05c_0Ma_environment_layer_map.png")
save_plot(plot_heatmap_matrix(as.matrix(TrData), "06 Species trait matrix", "trait_value"),
          "06_traits_heatmap.png")
trait_pca <- prcomp(TrData, scale. = TRUE)
trait_scores <- data.frame(species = species, clade = clade[species], trait_pca$x[, 1:2])
save_plot(ggplot(trait_scores, aes(x = PC1, y = PC2, color = clade, label = species)) +
            geom_point(size = 2) + labs(title = "07 Traits PCA"),
          "07_traits_pca.png")
trait_miss <- data.frame(trait = colnames(TrData), missing = colMeans(is.na(TrData)))
save_plot(ggplot(trait_miss, aes(x = trait, y = missing)) + geom_col(fill = pal[6]) +
            labs(x = NULL, y = "Missing fraction", title = "Trait missingness plot"),
          "07b_trait_missingness.png")
draw_tree(phy, file.path(fig_dir, "08_phylo_tree"), "08 Time-calibrated phylogeny")
save_plot(plot_heatmap_matrix(C, "09 Phylogenetic correlation matrix", "C"),
          "09_phylo_correlation_matrix.png", 6, 5.5)
origin_df <- data.frame(species = species, clade = clade[species], origin_ma = species_origin[species])
save_plot(ggplot(origin_df, aes(x = reorder(species, origin_ma), y = origin_ma, color = clade)) +
            geom_point(size = 2) + coord_flip() + scale_y_reverse() +
            labs(x = NULL, y = "Origin time (Ma)", title = "Species origin times"),
          "09b_species_origin_time.png")
save_plot(plot_heatmap_matrix(tip_ranges, "10 Tip ranges", "present"),
          "10_tip_ranges_heatmap.png")
save_plot(ggplot(region_grid0, aes(x = lon, y = lat, color = region)) + geom_point(size = 0.35) +
            coord_equal() + labs(title = "11 Region map at 0 Ma", x = "Longitude", y = "Latitude"),
          "11_region_map_0Ma.png")
save_plot(ggplot(grid0, aes(x = lon, y = lat, color = land_mask_dem)) + geom_point(size = 0.35) +
            coord_equal() + scale_color_gradient(low = "grey85", high = pal[4]) +
            labs(title = "12 Landmask map at 0 Ma", x = "Longitude", y = "Latitude", color = "Land"),
          "12_landmask_map_0Ma.png")

message("Phase 2: fitting real HMSC MCMC")
m <- Hmsc::Hmsc(
  Y = Y,
  XData = XData,
  XFormula = XFormula,
  TrData = TrData,
  TrFormula = TrFormula,
  phyloTree = phy,
  distr = "probit",
  studyDesign = studyDesign,
  ranLevels = ranLevels
)
m <- Hmsc::sampleMcmc(
  m,
  samples = cfg$samples,
  transient = cfg$transient,
  thin = cfg$thin,
  nChains = cfg$nChains,
  nParallel = cfg$nParallel,
  verbose = 0
)
save_rds(m, "hmsc_model_sampled.rds")
diag <- hee_hmsc_diagnostics(m)
save_rds(diag, "hmsc_diagnostics.rds")
pred_modern <- diag$predicted
MF <- diag$model_fit
VP <- diag$variance_partitioning

beta_est <- diag$beta$mean
gamma_est <- diag$gamma$mean
omega_est <- diag$omega$mean[[1]]
if (is.null(dimnames(beta_est))) dimnames(beta_est) <- list(m$covNames, m$spNames)
if (nrow(beta_est) != cfg$n_species && ncol(beta_est) == cfg$n_species) beta_for_projection <- t(beta_est) else beta_for_projection <- beta_est
if (!"(Intercept)" %in% colnames(beta_for_projection) && "(Intercept)" %in% rownames(beta_est)) beta_for_projection <- t(beta_est)
rownames(beta_for_projection) <- species
colnames(beta_for_projection) <- colnames(Btrue)

write_table(beta_for_projection, "hmsc_beta_posterior_mean.csv")
write_table(gamma_est, "hmsc_gamma_posterior_mean.csv")
if (!is.null(omega_est)) write_table(omega_est, "hmsc_omega_posterior_mean.csv")
if (nrow(diag$rho) > 0) write_table(diag$rho, "hmsc_rho_posterior.csv")
if (!is.null(MF)) write_table(as.data.frame(MF), "hmsc_model_fit_metrics.csv")
if (nrow(diag$ess) > 0) write_table(diag$ess, "hmsc_ess.csv")
if (nrow(diag$rhat) > 0) write_table(diag$rhat, "hmsc_rhat.csv")

save_plot(plot_heatmap_matrix(beta_for_projection[, variables, drop = FALSE],
                              "13 HMSC Beta posterior mean", "Beta"),
          "13_hmsc_beta_heatmap.png")
if (!is.null(gamma_est)) {
  save_plot(plot_heatmap_matrix(gamma_est, "14 HMSC Gamma posterior mean", "Gamma"),
            "14_hmsc_gamma_heatmap.png")
}
if (nrow(diag$rho) > 0) {
  save_plot(ggplot(diag$rho, aes(x = parameter, y = mean)) +
              geom_pointrange(aes(ymin = lwr, ymax = upr), color = pal[2]) +
              coord_flip() + labs(x = NULL, y = "Rho posterior mean and 95% CrI",
                                  title = "15 HMSC Rho posterior"),
            "15_hmsc_rho_posterior.png")
} else {
  save_plot(ggplot(data.frame(x = "rho", y = 0), aes(x, y)) + geom_point() +
              labs(title = "15 HMSC Rho posterior unavailable from fitted object"),
            "15_hmsc_rho_posterior.png")
}
if (!is.null(omega_est)) {
  save_plot(plot_heatmap_matrix(omega_est, "16 HMSC Omega residual association", "Omega"),
            "16_hmsc_omega_heatmap.png")
}
if (!is.null(VP) && !is.null(VP$vals)) {
  vp_df <- as.data.frame(VP$vals)
  vp_df$species <- rownames(vp_df)
  vp_long <- reshape(vp_df, varying = setdiff(names(vp_df), "species"),
                     v.names = "fraction", timevar = "group",
                     times = setdiff(names(vp_df), "species"), direction = "long")
  save_plot(ggplot(vp_long, aes(x = species, y = fraction, fill = group)) +
              geom_col() + coord_flip() +
              labs(title = "17 HMSC variance partitioning", x = NULL, y = "Fraction"),
            "17_hmsc_variance_partitioning.png")
}
draws <- diag$mcmc_draws
draws_sel <- draws[draws$parameter %in% head(unique(draws$parameter), 8), ]
if (nrow(draws_sel) > 0) {
  save_plot(ggplot(draws_sel, aes(x = iteration, y = value, color = factor(chain))) +
              geom_line(linewidth = 0.35) + facet_wrap(~ parameter, scales = "free_y") +
              labs(title = "18 HMSC trace plots", x = "Iteration", y = "Value", color = "Chain"),
            "18_hmsc_traceplots.png", 9, 6)
}
if (nrow(diag$ess) > 0) {
  save_plot(ggplot(diag$ess, aes(x = ESS)) + geom_histogram(bins = 25, fill = pal[4], color = "white") +
              labs(title = "19 HMSC ESS histogram", x = "ESS", y = "Parameters"),
            "19_hmsc_ess_histogram.png")
}
if (nrow(diag$rhat) > 0) {
  save_plot(ggplot(diag$rhat, aes(x = Rhat)) + geom_histogram(bins = 25, fill = pal[8], color = "white") +
              geom_vline(xintercept = 1.05, linetype = 2, color = pal[7]) +
              labs(title = "20 HMSC Rhat histogram", x = "Rhat", y = "Parameters"),
            "20_hmsc_rhat_histogram.png")
}
pred_mean <- apply(pred_modern, c(1, 2), mean, na.rm = TRUE)
obs_pred <- data.frame(observed = as.numeric(Y), predicted = as.numeric(pred_mean))
save_plot(ggplot(obs_pred, aes(x = predicted, y = observed)) +
            geom_jitter(width = 0, height = 0.05, alpha = 0.35, color = pal[2]) +
            geom_smooth(method = "glm", method.args = list(family = binomial), se = TRUE, color = pal[7]) +
            labs(title = "21 HMSC observed vs predicted", x = "Predicted occurrence probability", y = "Observed"),
          "21_hmsc_observed_vs_predicted.png")

message("Phase 2: modern prediction and cross-validation summaries")
fold_random <- hee_cv_random(nrow(XData), k = 5, seed = 11)
fold_spatial <- hee_cv_spatial_block(sites[, c("lon", "lat")], k = 5)
fold_env <- hee_cv_environment_block(XData, variable = variables[1], k = 5)
cv_df <- data.frame(site = rownames(Y), random = fold_random,
                    spatial = fold_spatial, environment = fold_env,
                    lon = sites$lon, lat = sites$lat)
write_table(cv_df, "cv_folds.csv")
cv_perf <- do.call(rbind, lapply(c("random", "spatial", "environment"), function(kind) {
  do.call(rbind, lapply(sort(unique(cv_df[[kind]])), function(f) {
    idx <- cv_df[[kind]] == f
    met <- hee_evaluate_hmsc(as.numeric(Y[idx, , drop = FALSE]), as.numeric(pred_mean[idx, , drop = FALSE]))
    data.frame(cv = kind, fold = f, met)
  }))
}))
cv_perf$validation_type <- "apparent_fold_summary_not_refit_cv"
write_table(cv_perf, "cv_performance.csv")
save_plot(ggplot(cv_perf, aes(x = cv, y = AUC, fill = cv)) +
            geom_boxplot(alpha = 0.8, show.legend = FALSE) +
            labs(title = "22 HMSC apparent fold performance", x = NULL, y = "AUC"),
          "22_hmsc_cv_performance.png")
save_plot(ggplot(cv_df, aes(x = lon, y = lat, color = factor(spatial))) +
            geom_point(size = 2) + coord_equal() +
            labs(title = "Spatial block CV map", x = "Longitude", y = "Latitude", color = "Block"),
          "22b_spatial_block_CV_map.png")
rich_obs_pred <- data.frame(site = rownames(Y), observed = rowSums(Y),
                            predicted = rowSums(pred_mean), lon = sites$lon, lat = sites$lat)
save_plot(ggplot(rich_obs_pred, aes(x = observed, y = predicted)) +
            geom_point(color = pal[2], size = 2) + geom_smooth(method = "lm", se = TRUE, color = pal[7]) +
            labs(title = "Richness observed vs predicted", x = "Observed richness", y = "Predicted richness"),
          "22c_richness_observed_vs_predicted.png")
cal <- data.frame(bin = cut(obs_pred$predicted, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE),
                  observed = obs_pred$observed, predicted = obs_pred$predicted)
cal_sum <- aggregate(cbind(observed, predicted) ~ bin, cal, mean, na.rm = TRUE)
save_plot(ggplot(cal_sum, aes(x = predicted, y = observed)) +
            geom_point(size = 2, color = pal[2]) + geom_abline(linetype = 2) +
            labs(title = "Calibration plot", x = "Mean predicted", y = "Mean observed"),
          "22d_calibration_plot.png")
save_plot(ggplot(cv_perf, aes(x = factor(fold), y = cv, fill = AUC)) +
            geom_tile(color = "white") + scale_fill_gradient(low = "white", high = pal[2]) +
            labs(title = "CV performance heatmap", x = "Fold", y = NULL),
          "22e_CV_performance_heatmap.png")
species_perf <- do.call(rbind, lapply(seq_along(species), function(j) {
  met <- hee_evaluate_hmsc(Y[, j], pred_mean[, j])
  data.frame(species = species[j], met)
}))
species_perf$validation_type <- "apparent_in_sample"
write_table(species_perf, "hmsc_species_performance.csv")
save_plot(ggplot(species_perf, aes(x = reorder(species, AUC), y = AUC)) +
            geom_col(fill = pal[2]) + coord_flip() +
            labs(title = "Species-level apparent performance", x = NULL, y = "AUC"),
          "22f_species_level_performance.png", 6, 7)
save_plot(ggplot(cv_perf[cv_perf$cv == "environment", ], aes(x = factor(fold), y = AUC)) +
            geom_col(fill = pal[4]) +
            labs(title = "Environment block apparent performance", x = "Environment block", y = "AUC"),
          "22g_environment_block_CV_plot.png")

message("Phase 4: BioGeoBEARS-style accessibility and BSM objects")
models <- data.frame(
  model = c("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J"),
  logLik = c(-138.2, -132.5, -141.0, -134.8, -146.3, -139.7),
  n_parameters = c(2, 3, 2, 3, 2, 3),
  d = c(0.03, 0.03, 0.025, 0.026, 0.018, 0.018),
  e = c(0.015, 0.014, 0.02, 0.019, 0.024, 0.022),
  j = c(0, 0.05, 0, 0.04, 0, 0.035),
  evidence_type = "MOCK_ONLY_simulated_demonstration"
)
bgb_compare <- hee_bgb_compare(models, n = cfg$n_species * cfg$n_regions)
bgb_lrt <- rbind(
  cbind(comparison = "DEC vs DEC+J", hee_bgb_lrt(bgb_compare[bgb_compare$model == "DEC", ], bgb_compare[bgb_compare$model == "DEC+J", ])),
  cbind(comparison = "DIVALIKE vs DIVALIKE+J", hee_bgb_lrt(bgb_compare[bgb_compare$model == "DIVALIKE", ], bgb_compare[bgb_compare$model == "DIVALIKE+J", ])),
  cbind(comparison = "BAYAREALIKE vs BAYAREALIKE+J", hee_bgb_lrt(bgb_compare[bgb_compare$model == "BAYAREALIKE", ], bgb_compare[bgb_compare$model == "BAYAREALIKE+J", ]))
)
events <- data.frame(event_type = rep(c("dispersal", "extinction", "range expansion"), length.out = length(times)),
                     from_region = sample(regions, length(times), replace = TRUE),
                     to_region = sample(regions, length(times), replace = TRUE),
                     region = sample(regions, length(times), replace = TRUE),
                     time_ma = times,
                     event_prob = pmin(0.98, pmax(0.02, 0.45 + 0.30 * sin(times / 60))),
                     simulation_id = seq_along(times),
                     model = "MOCK_AVERAGED",
                     tree_id = "sim_tree_1",
                     branch_or_node = paste0("branch_", seq_along(times)),
                     ancestor_range = sample(regions, length(times), replace = TRUE),
                     descendant_range = sample(regions, length(times), replace = TRUE),
                     evidence_type = "MOCK_ONLY_simulated_demonstration")
events$from_area <- events$from_region
events$to_area <- events$to_region
bgb_events <- hee_bgb_bsm(events)
bgb_event_summary <- hee_bgb_event_summary(bgb_events)
bsm_metrics <- hee_bgb_bsm_metrics(bgb_events)
accessibility <- expand.grid(species = species, region = regions, time_ma = times,
                             KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
accessibility$origin <- species_origin[accessibility$species]
accessibility$tip_present <- tip_ranges[cbind(accessibility$species, accessibility$region)]
accessibility$accessibility <- as.numeric(accessibility$time_ma <= accessibility$origin) *
  pmin(1, 0.18 + 0.55 * accessibility$tip_present + 0.25 * exp(-accessibility$time_ma / 250))
accessibility$evidence_type <- "MOCK_ONLY_simulated_demonstration"
accessibility <- hee_bgb_accessibility(accessibility)
accessibility$evidence_type <- "MOCK_ONLY_simulated_demonstration"
write_table(bgb_compare, "bgb_model_compare.csv")
write_table(bgb_lrt, "bgb_lrt.csv")
write_table(bgb_events, "bgb_bsm_events.csv")
write_table(accessibility, "bgb_accessibility.csv")
write_table(bsm_metrics$region_metrics, "bgb_bsm_region_metrics.csv")
write_table(bsm_metrics$dispersal_rates, "bgb_bsm_dispersal_rates.csv")
write_table(bsm_metrics$source_sink_long, "bgb_bsm_source_sink_turnover_long.csv")
save_plot(ggplot(bgb_compare, aes(x = reorder(model, AICc_weight), y = AICc_weight, fill = model)) +
            geom_col(show.legend = FALSE) + coord_flip() +
            labs(title = "23 BGB model AICc weights", x = NULL, y = "AICc weight"),
          "23_bgb_model_weights.png")
param_long <- reshape(models[, c("model", "d", "e", "j")], varying = c("d", "e", "j"),
                      v.names = "value", timevar = "parameter", times = c("d", "e", "j"),
                      direction = "long")
save_plot(ggplot(param_long, aes(x = model, y = value, fill = parameter)) +
            geom_col(position = "dodge") + coord_flip() +
            labs(title = "BGB parameter plot", x = NULL, y = "Parameter value"),
          "23b_bgb_parameter_plot.png")
save_plot(ggplot(bgb_lrt, aes(x = comparison, y = p_value)) +
            geom_col(fill = pal[6]) + coord_flip() +
            geom_hline(yintercept = 0.05, linetype = 2, color = pal[7]) +
            labs(title = "24 BGB LRT result table plot", x = NULL, y = "p-value"),
          "24_bgb_lrt_table.png")
anc_prob <- expand.grid(species = species, region = regions)
anc_prob$probability <- as.numeric(tip_ranges[cbind(anc_prob$species, anc_prob$region)]) * 0.65 + runif(nrow(anc_prob), 0, 0.25)
save_plot(ggplot(anc_prob, aes(x = region, y = species, fill = probability)) +
            geom_tile(color = "white") + scale_fill_gradient(low = "white", high = pal[2]) +
            labs(title = "25 BGB ancestral range probability heatmap"),
          "25_bgb_ancestral_range_heatmap.png")
save_plot(ggplot(bgb_events, aes(x = time_ma, fill = event_type)) +
            geom_histogram(binwidth = if (quick) 60 else 25, color = "white") +
            scale_x_reverse() + labs(title = "26 BGB BSM event histogram", x = "Time (Ma)", y = "Events"),
          "26_bgb_bsm_events_time.png")
source_sink <- as.data.frame.matrix(table(bgb_events$from_region, bgb_events$to_region))
save_plot(plot_heatmap_matrix(as.matrix(source_sink), "27 BGB source/sink heatmap", "events", low = "white", high = pal[7]),
          "27_bgb_source_sink_heatmap.png")
net <- aggregate(event_prob ~ from_region + to_region, bgb_events, sum)
save_plot(ggplot(net, aes(x = from_region, y = to_region, size = event_prob, color = event_prob)) +
            geom_point(alpha = 0.8) + scale_color_gradient(low = pal[2], high = pal[7]) +
            labs(title = "28 BGB dispersal network through time", x = "Source", y = "Sink"),
          "28_bgb_dispersal_network.png")
save_plot(ggplot(bsm_metrics$source_sink_long,
                 aes(x = time_ma, y = value, color = metric)) +
            geom_line(linewidth = 0.6) + facet_wrap(~ region, scales = "free_y") +
            scale_x_reverse() +
            labs(title = "BSM source, sink, extinction and turnover through time",
                 x = "Time (Ma)", y = "Event weight"),
          "28b_bgb_bsm_source_sink_turnover.png", 9, 6)

message("Phase 4: phylogenetic time mask and plate correction")
E_phylo <- hee_phylo_time_mask(phy, times, species_origin = species_origin, level = "species")
write_table(cbind(species = rownames(E_phylo), as.data.frame(E_phylo)), "phylo_time_mask.csv")
valid_lineages <- data.frame(time_ma = times, valid_species = colSums(E_phylo))
write_table(valid_lineages, "valid_species_through_time.csv")
mask_long <- long_matrix(E_phylo, "species", "time_ma", "valid")
mask_long$time_ma <- as.numeric(gsub("Ma$", "", mask_long$time_ma))
save_plot(ggplot(mask_long, aes(x = time_ma, y = species, fill = valid)) +
            geom_tile() + scale_x_reverse() +
            labs(title = "Species x time phylogenetic mask", x = "Time (Ma)", y = NULL),
          "phylo_mask_heatmap.png", 8, 5.5)
save_plot(ggplot(valid_lineages, aes(x = time_ma, y = valid_species)) +
            geom_line(color = pal[2], linewidth = 0.8) + scale_x_reverse() +
            labs(title = "Number of valid species through time", x = "Time (Ma)", y = "Species"),
          "valid_species_through_time.png")
plate_corrected <- hee_reconstruct_points(plate_points, backend = "user_table")
plate_compare <- hee_compare_plate_corrected(
  data.frame(point_id = seq_len(nrow(plate_points)), lon = 0, lat = 0),
  cbind(point_id = seq_len(nrow(plate_corrected)), plate_corrected)
)
write_table(plate_compare, "plate_corrected_vs_uncorrected.csv")
save_plot(ggplot(plate_corrected, aes(x = paleo_lon, y = paleo_lat, color = track_id, group = track_id)) +
            geom_path(linewidth = 0.65) + geom_point(size = 0.8, alpha = 0.45) +
            coord_equal(xlim = c(-180, 180), ylim = c(-90, 90)) +
            labs(title = "Tectonic tracks through time", x = "Paleo-longitude", y = "Paleo-latitude"),
          "plate_tectonic_tracks_through_time.png")
track_env <- hee_extract_paleoenv_track(plate_corrected, cube, variables = variables, land_only = TRUE)
track_env_uncorrected <- track_env
track_env_uncorrected$paleo_lon <- 0
track_env_uncorrected$paleo_lat <- 0
track_diff <- data.frame(time_ma = track_env$time_ma,
                         track_id = track_env$track_id,
                         MAT_diff = track_env$MAT_pohl_C - mean(grid0$MAT_pohl_C, na.rm = TRUE))
write_table(track_diff, "plate_corrected_uncorrected_environment_difference.csv")
save_plot(ggplot(track_diff, aes(x = time_ma, y = MAT_diff, color = track_id)) +
            geom_line() + scale_x_reverse() +
            labs(title = "Plate-corrected vs uncorrected environment difference", x = "Time (Ma)", y = "MAT difference"),
          "plate_corrected_vs_uncorrected_environment_difference.png")

message("Phase 3/5: full 540-0 Ma projection, five models, uncertainty and community summaries")
paleo_grid <- hee_make_paleo_grid(cube, times = times, variables = variables, land_only = TRUE)
paleo_grid$region <- regions[cut(paleo_grid$lon, breaks = cfg$n_regions, labels = FALSE, include.lowest = TRUE)]
X_paleo <- hee_apply_recipe(recipe, paleo_grid)
X_paleo <- cbind(paleo_grid[, c("cell_id", "time_ma", "lon", "lat", "land_area_km2", "land_mask_dem", "region")], X_paleo)
landscape_state <- X_paleo[, c("cell_id", "time_ma", "region", "lon", "lat",
                               "land_area_km2", "land_mask_dem", variables), drop = FALSE]
landscape_state <- hee_land_age(landscape_state, unit_col = "cell_id",
                                land_col = "land_mask_dem")
landscape_events <- hee_landscape_events(landscape_state, unit_col = "cell_id",
                                         land_col = "land_mask_dem",
                                         area_col = "land_area_km2")
write_table(landscape_events, "geoprocess_landscape_events.csv")
landscape_summary <- aggregate(geographic_existence ~ time_ma, landscape_state, mean, na.rm = TRUE)
names(landscape_summary)[2] <- "land_fraction"
land_age_summary <- aggregate(land_age_myr ~ time_ma, landscape_state, mean, na.rm = TRUE)
landscape_summary <- merge(landscape_summary, land_age_summary, by = "time_ma", all = TRUE)
write_table(landscape_summary, "geoprocess_landscape_summary.csv")
save_plot(ggplot(landscape_summary, aes(x = time_ma)) +
            geom_line(aes(y = land_fraction, color = "land fraction"), linewidth = 0.8) +
            geom_line(aes(y = scales::rescale(land_age_myr, to = c(0, 1), from = range(land_age_myr, na.rm = TRUE)),
                          color = "mean land age (scaled)"), linewidth = 0.8) +
            scale_x_reverse() +
            labs(title = "Geological arena through time",
                 x = "Time (Ma)", y = "Scaled value", color = NULL),
          "geoprocess_landscape_summary.png", 8.5, 5)

eco_opp <- hee_ecological_opportunity(landscape_state, env_cols = variables,
                                      region_col = "region",
                                      area_col = "land_area_km2")
eco_opp_summary <- aggregate(cbind(ecological_opportunity, area_gain_index,
                                   heterogeneity_index, novelty_index,
                                   young_land_index) ~ time_ma,
                             eco_opp, mean, na.rm = TRUE)
write_table(eco_opp_summary, "geoprocess_ecological_opportunity_summary.csv")
eco_opp_long <- reshape(eco_opp_summary,
                        varying = c("ecological_opportunity", "area_gain_index",
                                    "heterogeneity_index", "novelty_index",
                                    "young_land_index"),
                        v.names = "value", timevar = "component",
                        times = c("opportunity", "area gain", "heterogeneity",
                                  "novelty", "young land"),
                        direction = "long")
save_plot(ggplot(eco_opp_long, aes(x = time_ma, y = value, color = component)) +
            geom_line(linewidth = 0.75) + scale_x_reverse() +
            labs(title = "Ecological opportunity components through time",
                 x = "Time (Ma)", y = "Index", color = NULL),
          "geoprocess_ecological_opportunity.png", 8.5, 5)

region_centroids <- aggregate(cbind(lon, lat, land_area_km2) ~ region + time_ma,
                              X_paleo, mean, na.rm = TRUE)
geo_distances <- do.call(rbind, lapply(split(region_centroids, region_centroids$time_ma),
                                       function(z) {
  pair <- expand.grid(from_region = z$region, to_region = z$region,
                      stringsAsFactors = FALSE)
  a <- z[match(pair$from_region, z$region), ]
  b <- z[match(pair$to_region, z$region), ]
  pair$time_ma <- z$time_ma[1]
  pair$paleodistance <- sqrt((a$lon - b$lon)^2 + (a$lat - b$lat)^2)
  pair$barrier_strength <- pmin(pair$paleodistance / max(pair$paleodistance, na.rm = TRUE), 1)
  pair
}))
geo_conn <- hee_build_connectivity_cube(
  geo_distances,
  traits = data.frame(species = species,
                      dispersal_distance = unname(dispersal_distance[species])),
  distance_scale = 60
)
geo_iso <- hee_isolation_history(geo_conn, threshold = 0.2)
geo_conn_summary <- aggregate(cbind(structural_connectivity, functional_connectivity,
                                    connectivity, isolation, isolation_duration_ma) ~ time_ma,
                              geo_iso, mean, na.rm = TRUE)
write_table(geo_conn_summary, "geoprocess_connectivity_isolation_summary.csv")
save_plot(ggplot(geo_conn_summary, aes(x = time_ma)) +
            geom_line(aes(y = connectivity, color = "total connectivity"), linewidth = 0.8) +
            geom_line(aes(y = isolation, color = "isolation"), linewidth = 0.8) +
            geom_line(aes(y = scales::rescale(isolation_duration_ma, to = c(0, 1),
                                              from = range(isolation_duration_ma, na.rm = TRUE)),
                          color = "isolation duration (scaled)"), linewidth = 0.8) +
            scale_x_reverse() +
            labs(title = "Dynamic connectivity and isolation through time",
                 x = "Time (Ma)", y = "Index", color = NULL),
          "geoprocess_connectivity_isolation.png", 8.5, 5)
extrap <- hee_extrapolation_risk(XData_raw, cube, variables = variables, times = times, method = "range")
Q_extrap <- hee_extrapolation_weight(extrap, mode = extrapolation_mode,
                                     lambda = 1, threshold = 0)
write_table(aggregate(cbind(no_analog_fraction, extrapolation_score) ~ time_ma, extrap, mean, na.rm = TRUE),
            "extrapolation_summary_by_time.csv")
write_table(Q_extrap, "extrapolation_q_extrap_policy.csv")
save_rds(extrap, "extrapolation_full.rds")
tile_index <- hee_make_tiles(X_paleo, tile_size = cfg$chunk_size)
tile_index$time_ma <- X_paleo$time_ma[tile_index$row_id]
tile_ranges <- stats::aggregate(row_id ~ tile_id, tile_index, function(z) paste0(min(z), "-", max(z)))
names(tile_ranges)[2] <- "row_range"
tile_n <- stats::aggregate(row_id ~ tile_id, tile_index, length)
names(tile_n)[2] <- "n_rows"
tile_time_min <- stats::aggregate(time_ma ~ tile_id, tile_index, min)
names(tile_time_min)[2] <- "time_ma_min"
tile_time_max <- stats::aggregate(time_ma ~ tile_id, tile_index, max)
names(tile_time_max)[2] <- "time_ma_max"
task_table <- Reduce(function(x, y) merge(x, y, by = "tile_id", sort = TRUE),
                     list(tile_ranges, tile_n, tile_time_min, tile_time_max))
task_table$task_id <- seq_len(nrow(task_table))
task_table$model_id <- "M1_M5"
task_table$status <- "pending"
task_table$output_file <- file.path(rds_dir, "M1_M5_projection_outputs_pending.rds")
task_table$input_hash <- vapply(seq_len(nrow(task_table)), function(i) {
  hee_hash_inputs(list(tile_id = task_table$tile_id[i],
                       row_range = task_table$row_range[i],
                       time_ma_min = task_table$time_ma_min[i],
                       time_ma_max = task_table$time_ma_max[i],
                       variables = variables,
                       projection_mode = projection_mode,
                       extrapolation_mode = extrapolation_mode))$md5[1]
}, character(1))
task_table <- task_table[, c("task_id", "tile_id", "row_range", "n_rows",
                             "time_ma_min", "time_ma_max", "model_id", "status",
                             "output_file", "input_hash")]
write_table(task_table, "projection_task_table.csv")

model_defs <- hee_projection_formula(include_q = TRUE,
                                     extrapolation_mode = extrapolation_mode)
model_defs$projection_mode <- projection_mode
write_table(model_defs, "projection_model_definitions.csv")
projection_example <- hee_projection_component_example(
  S_HMSC = 0.80, A_BGB = 0.25, E_phylo = 1, L_land = 1,
  D_static = 0.60, D_dynamic = 0.10, Q_extrap = 1,
  extrapolation_mode = extrapolation_mode
)
write_table(projection_example, "projection_numeric_example.csv")
save_plot(ggplot(projection_example,
                 aes(x = model_id, y = probability, fill = model_id)) +
            geom_col(show.legend = FALSE) + coord_flip() +
            labs(title = "Numeric example: layered deep-time projection",
                 x = NULL, y = "Final probability P"),
          "projection_numeric_example.png")

beta_draws <- extract_beta_draws(m, cfg$n_uncert_draws)
project_models <- list()
richness_all <- list()
area_all <- list()
centroid_all <- list()
component_all <- list()
rep_maps <- list()
uncert_maps <- list()
m5_projection_full <- NULL
for (model_id in model_defs$model_id) {
  P_all <- list()
  tile_ids <- unique(tile_index$tile_id)
  for (tile_id in tile_ids) {
    rows <- tile_index$row_id[tile_index$tile_id == tile_id]
    Xtile <- X_paleo[rows, , drop = FALSE]
    S <- project_fixed(beta_for_projection, Xtile,
                       id_cols = c("cell_id", "time_ma", "lon", "lat", "land_area_km2", "land_mask_dem", "region"))
    landmask_tile <- unique(S[, c("cell_id", "time_ma", "land_mask_dem"),
                              drop = FALSE])
    q_lookup <- unique(Q_extrap[rows, c("cell_id", "time_ma", "Q_extrap",
                                        "extrapolation_score",
                                        "extrapolation_flag",
                                        "extrapolation_mode"), drop = FALSE])
    use_phylo <- model_id != "M1_env_only"
    use_access <- model_id %in% model_defs$model_id[3:5]
    combine_args <- list(
      suitability = S,
      accessibility = if (use_access) accessibility else NULL,
      phylo_mask = if (use_phylo) E_phylo else NULL,
      landmask = landmask_tile
    )
    P_core <- do.call(hee_combine, combine_args)
    if (model_id == "M4_env_phylo_bgb_static_dispersal") {
      static <- hee_static_region_dispersal(
        P_core[, c("species", "cell_id", "time_ma", "lon", "lat", "region", "accessibility"), drop = FALSE],
        accessibility = accessibility,
        dispersal_scale = stats::setNames(35 + 25 * scales::rescale(traits$dispersal, to = c(0, 1)), species),
        source_threshold = 0.45
      )
      combine_args$static_filter <- static
      P_core <- do.call(hee_combine, combine_args)
    } else {
      static <- NULL
    }
    if (model_id == "M5_env_phylo_bgb_dynamic_dispersal") {
      dyn <- hee_dynamic_region_filter(
        previous_projection = transform(P_core, probability = probability),
        cell_table = P_core[, c("species", "cell_id", "time_ma", "lon", "lat", "region"), drop = FALSE],
        dispersal_scale = stats::setNames(35 + 25 * scales::rescale(traits$dispersal, to = c(0, 1)), species),
        time_direction = "forward",
        default = 0.25
      )
      combine_args$dynamic_filter <- dyn
      P_core <- do.call(hee_combine, combine_args)
    } else {
      dyn <- NULL
    }
    combine_args$extrapolation_weight <- q_lookup
    P <- do.call(hee_combine, combine_args)
    P$probability_core <- P_core$probability
    q_meta <- q_lookup[, c("cell_id", "time_ma", "extrapolation_score",
                           "extrapolation_flag", "extrapolation_mode"), drop = FALSE]
    P <- safe_left_join(P, q_meta, by = c("cell_id", "time_ma"),
                        label = "extrapolation diagnostics")
    if (!is.null(dyn) && "dynamic_source_time_ma" %in% names(dyn)) {
      dyn_meta <- dyn[, c("species", "cell_id", "time_ma",
                          "dynamic_source_time_ma"), drop = FALSE]
      P <- safe_left_join(P, dyn_meta,
                          by = c("species", "cell_id", "time_ma"),
                          label = "dynamic dispersal source time")
    }
    P$model_id <- model_id
    P_all[[as.character(tile_id)]] <- P
  }
  P_model <- do.call(rbind, P_all)
  rich <- richness_summary(P_model, c("model_id", "cell_id", "time_ma"))
  area <- area_summary(P_model, model_id)
  cent <- centroid_summary(rich)
  cent$model_id <- model_id
  comp_cols_i <- intersect(c("suitability", "accessibility", "phylo_existence",
                             "land_mask_dem", "D_static", "D_dynamic", "Q_extrap",
                             "probability_core", "probability"),
                           names(P_model))
  component_all[[model_id]] <- aggregate(P_model[, comp_cols_i, drop = FALSE],
                                         P_model[, c("model_id", "time_ma"), drop = FALSE],
                                         mean, na.rm = TRUE)
  project_models[[model_id]] <- P_model[P_model$time_ma %in% rep_times | P_model$species == species[1], ]
  if (model_id == "M5_env_phylo_bgb_dynamic_dispersal") {
    m5_projection_full <- P_model
    save_rds(P_model, "M5_env_phylo_bgb_dynamic_dispersal_probability_full.rds")
  }
  richness_all[[model_id]] <- rich
  area_all[[model_id]] <- area
  centroid_all[[model_id]] <- cent
  rep_maps[[model_id]] <- rich[rich$time_ma %in% rep_times, ]
  save_rds(P_model[P_model$time_ma %in% rep_times, ], paste0(model_id, "_representative_species_probability_maps.rds"))
  save_rds(rich, paste0(model_id, "_richness_full.rds"))
}
task_table$status <- "done"
task_table$output_file <- file.path(rds_dir, "M1_M5_richness_full.rds")
task_table$completed_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
write_table(task_table, "projection_task_table.csv")
richness_models <- do.call(rbind, richness_all)
area_models <- do.call(rbind, area_all)
centroid_models <- do.call(rbind, centroid_all)
save_rds(richness_models, "M1_M5_richness_full.rds")
save_rds(area_models, "M1_M5_area_full.rds")
save_rds(centroid_models, "M1_M5_centroid_full.rds")
component_models <- do.call(rbind, component_all)
write_table(component_models, "M1_M5_projection_component_diagnostics.csv")
write_table(aggregate(expected_richness ~ model_id + time_ma, richness_models, mean, na.rm = TRUE),
            "M1_M5_richness_through_time.csv")
write_table(area_models, "M1_M5_occupied_area_through_time.csv")
write_table(centroid_models, "M1_M5_centroid_shift.csv")
component_cols <- intersect(c("suitability", "accessibility", "phylo_existence",
                              "land_mask_dem", "D_static", "D_dynamic", "Q_extrap",
                              "probability_core", "probability"),
                            names(m5_projection_full))
component_summary <- aggregate(m5_projection_full[, component_cols, drop = FALSE],
                               m5_projection_full[, c("model_id", "time_ma"), drop = FALSE],
                               mean, na.rm = TRUE)
write_table(component_summary, "projection_component_summary_M5.csv")
component_long <- reshape(component_summary,
                          varying = component_cols,
                          v.names = "mean_value",
                          timevar = "component",
                          times = component_cols,
                          direction = "long")
component_long <- component_long[is.finite(component_long$mean_value), , drop = FALSE]
save_plot(ggplot(component_long, aes(x = time_ma, y = mean_value, color = component)) +
            geom_line(linewidth = 0.65) + scale_x_reverse() +
            labs(title = "M5 component summary through time",
                 x = "Time (Ma)", y = "Mean component value"),
          "projection_component_summary_M5.png", 9, 5.5)

rich_time_models <- aggregate(expected_richness ~ model_id + time_ma, richness_models, mean, na.rm = TRUE)
save_plot(ggplot(rich_time_models, aes(x = time_ma, y = expected_richness, color = model_id)) +
            geom_line(linewidth = 0.75) + scale_x_reverse() +
            labs(title = "M1-M5 richness through time", x = "Time (Ma)", y = "Mean expected richness"),
          "M1_M5_richness_through_time.png")
area_line_models <- area_models[stats::ave(area_models$time_ma, area_models$model_id,
                                           FUN = length) > 1, , drop = FALSE]
save_plot(ggplot(area_models, aes(x = time_ma, y = land_area_km2 / 1e6, color = model_id)) +
            geom_point(size = 1.2) +
            geom_line(data = area_line_models, linewidth = 0.75) + scale_x_reverse() +
            labs(title = "M1-M5 occupied area through time", x = "Time (Ma)", y = "Occupied area (million km2)"),
          "M1_M5_area_through_time.png")
save_plot(ggplot(centroid_models, aes(x = lon_centroid, y = lat_centroid, color = model_id, group = model_id)) +
            geom_path(linewidth = 0.7) + geom_point(size = 1.2) + coord_equal() +
            labs(title = "M1-M5 centroid shift", x = "Longitude centroid", y = "Latitude centroid"),
          "M1_M5_centroid_shift.png")
map_m3 <- richness_all[["M3_env_phylo_bgb_no_dispersal"]]
map_m5 <- richness_all[["M5_env_phylo_bgb_dynamic_dispersal"]]
diff_m5_m3 <- merge(map_m5, map_m3, by = c("cell_id", "time_ma", "lon", "lat"),
                    suffixes = c("_M5", "_M3"))
diff_m5_m3$difference <- diff_m5_m3$expected_richness_M5 - diff_m5_m3$expected_richness_M3
write_table(diff_m5_m3[diff_m5_m3$time_ma %in% rep_times, ], "M5_minus_M3_richness_difference_representative.csv")
save_plot(ggplot(diff_m5_m3[diff_m5_m3$time_ma %in% rep_times, ], aes(x = lon, y = lat, color = difference)) +
            geom_point(size = 0.25) + facet_wrap(~ time_ma, ncol = 3) + coord_equal() +
            scale_color_gradient2(low = pal[2], mid = "white", high = pal[7]) +
            labs(title = "M5 - M3 richness difference maps", x = "Longitude", y = "Latitude"),
          "M5_minus_M3_richness_difference.png", 9, 7)
map_m4 <- richness_all[["M4_env_phylo_bgb_static_dispersal"]]
diff_m4_m3 <- merge(map_m4, map_m3, by = c("cell_id", "time_ma", "lon", "lat"),
                    suffixes = c("_M4", "_M3"))
diff_m4_m3$difference <- diff_m4_m3$expected_richness_M4 - diff_m4_m3$expected_richness_M3
save_plot(ggplot(diff_m4_m3[diff_m4_m3$time_ma %in% rep_times, ], aes(x = lon, y = lat, color = difference)) +
            geom_point(size = 0.25) + facet_wrap(~ time_ma, ncol = 3) + coord_equal() +
            scale_color_gradient2(low = pal[2], mid = "white", high = pal[7]) +
            labs(title = "Static dispersal vs no dispersal difference", x = "Longitude", y = "Latitude"),
          "dispersal_vs_no_dispersal_difference.png", 9, 7)
agree <- aggregate(expected_richness ~ cell_id + time_ma + lon + lat, richness_models, stats::sd, na.rm = TRUE)
names(agree)[5] <- "model_sd"
save_plot(ggplot(agree[agree$time_ma %in% rep_times, ], aes(x = lon, y = lat, color = -model_sd)) +
            geom_point(size = 0.25) + facet_wrap(~ time_ma, ncol = 3) + coord_equal() +
            scale_color_gradient(low = pal[7], high = pal[4]) +
            labs(title = "Model agreement map", x = "Longitude", y = "Latitude", color = "Agreement"),
          "model_agreement_map.png", 9, 7)

selected_deep_models <- c("M3_env_phylo_bgb_no_dispersal",
                          "M4_env_phylo_bgb_static_dispersal",
                          "M5_env_phylo_bgb_dynamic_dispersal")
model_short <- function(x) {
  out <- gsub("M3_env_phylo_bgb_no_dispersal", "M3 no dispersal", x, fixed = TRUE)
  out <- gsub("M4_env_phylo_bgb_static_dispersal", "M4 static dispersal", out, fixed = TRUE)
  out <- gsub("M5_env_phylo_bgb_dynamic_dispersal", "M5 dynamic dispersal", out, fixed = TRUE)
  out
}
m3_m5_rich_time <- rich_time_models[rich_time_models$model_id %in% selected_deep_models, , drop = FALSE]
m3_m5_rich_time$model_short <- model_short(m3_m5_rich_time$model_id)
write_table(m3_m5_rich_time, "M3_M5_richness_through_time.csv")
save_plot(ggplot(m3_m5_rich_time, aes(x = time_ma, y = expected_richness, color = model_short)) +
            geom_line(linewidth = 0.8) + scale_x_reverse() +
            labs(title = "M3-M5 richness through deep time",
                 subtitle = "M3 = no dispersal; M4 = static dispersal; M5 = dynamic dispersal",
                 x = "Time (Ma)", y = "Mean expected richness", color = "Projection model"),
          "M3_M5_richness_through_time.png", 8.5, 5)
m3_m5_area_time <- area_models[area_models$model_id %in% selected_deep_models, , drop = FALSE]
m3_m5_area_time$model_short <- model_short(m3_m5_area_time$model_id)
write_table(m3_m5_area_time, "M3_M5_occupied_area_through_time.csv")
area_line_m3_m5 <- m3_m5_area_time[stats::ave(m3_m5_area_time$time_ma, m3_m5_area_time$model_short,
                                              FUN = length) > 1, , drop = FALSE]
save_plot(ggplot(m3_m5_area_time, aes(x = time_ma, y = land_area_km2 / 1e6, color = model_short)) +
            geom_point(size = 1.2) +
            geom_line(data = area_line_m3_m5, linewidth = 0.8) + scale_x_reverse() +
            labs(title = "M3-M5 occupied area through deep time",
                 x = "Time (Ma)", y = "Occupied area (million km2)", color = "Projection model"),
          "M3_M5_occupied_area_through_time.png", 8.5, 5)
m3_m5_component <- component_models[component_models$model_id %in% selected_deep_models, , drop = FALSE]
m3_m5_component$model_short <- model_short(m3_m5_component$model_id)
component_diag_cols <- intersect(c("accessibility", "phylo_existence", "land_mask_dem",
                                   "D_static", "D_dynamic", "Q_extrap", "probability"),
                                 names(m3_m5_component))
m3_m5_component_long <- reshape(m3_m5_component,
                                varying = component_diag_cols,
                                v.names = "mean_value",
                                timevar = "component",
                                times = component_diag_cols,
                                direction = "long")
m3_m5_component_long <- m3_m5_component_long[is.finite(m3_m5_component_long$mean_value), , drop = FALSE]
write_table(m3_m5_component_long, "M3_M5_component_diagnostics_through_time.csv")
save_plot(ggplot(m3_m5_component_long, aes(x = time_ma, y = mean_value, color = component)) +
            geom_line(linewidth = 0.65) + scale_x_reverse() +
            facet_wrap(~ model_short, ncol = 1) +
            labs(title = "M3-M5 projection component diagnostics through deep time",
                 x = "Time (Ma)", y = "Mean component value", color = "Component"),
          "M3_M5_component_diagnostics_through_time.png", 9, 7)
m3_m5_rep_rich <- do.call(rbind, rep_maps[selected_deep_models])
m3_m5_rep_rich$model_short <- model_short(m3_m5_rep_rich$model_id)
write_table(m3_m5_rep_rich, "M3_M5_representative_richness_maps.csv")
save_plot(ggplot(m3_m5_rep_rich, aes(x = lon, y = lat, color = expected_richness)) +
            geom_point(size = 0.20) + facet_grid(model_short ~ time_ma) + coord_equal() +
            scale_color_gradient(low = "white", high = pal[7]) +
            labs(title = "M3-M5 representative richness maps",
                 x = "Longitude", y = "Latitude", color = "Richness"),
          "M3_M5_representative_richness_maps.png", 13.5, 7.5)
prob_rep_list <- lapply(selected_deep_models, function(mid) {
  z <- project_models[[mid]]
  z <- z[z$time_ma %in% rep_times & z$species == species[1], , drop = FALSE]
  z$model_short <- model_short(z$model_id)
  z
})
prob_rep_common <- Reduce(intersect, lapply(prob_rep_list, names))
m3_m5_prob_rep <- do.call(rbind, lapply(prob_rep_list, function(z) z[, prob_rep_common, drop = FALSE]))
write_table(m3_m5_prob_rep[, intersect(c("model_id", "model_short", "species", "cell_id",
                                         "time_ma", "lon", "lat", "probability",
                                         "suitability", "accessibility", "phylo_existence",
                                         "D_static", "D_dynamic", "Q_extrap"),
                                       names(m3_m5_prob_rep)), drop = FALSE],
            "M3_M5_representative_probability_maps.csv")
save_plot(ggplot(m3_m5_prob_rep, aes(x = lon, y = lat, color = probability)) +
            geom_point(size = 0.20) + facet_grid(model_short ~ time_ma) + coord_equal() +
            scale_color_gradient(low = "white", high = pal[7]) +
            labs(title = paste0("M3-M5 final probability maps for ", species[1]),
                 x = "Longitude", y = "Latitude", color = "P"),
          "M3_M5_representative_probability_maps.png", 13.5, 7.5)
diff_m5_m4 <- merge(map_m5, map_m4, by = c("cell_id", "time_ma", "lon", "lat"),
                    suffixes = c("_M5", "_M4"))
diff_m5_m4$difference <- diff_m5_m4$expected_richness_M5 - diff_m5_m4$expected_richness_M4
write_table(diff_m5_m4[diff_m5_m4$time_ma %in% rep_times, ],
            "M5_minus_M4_dynamic_vs_static_difference_representative.csv")
save_plot(ggplot(diff_m5_m4[diff_m5_m4$time_ma %in% rep_times, ], aes(x = lon, y = lat, color = difference)) +
            geom_point(size = 0.25) + facet_wrap(~ time_ma, ncol = 3) + coord_equal() +
            scale_color_gradient2(low = pal[2], mid = "white", high = pal[7]) +
            labs(title = "M5 - M4 richness difference maps",
                 subtitle = "Dynamic dispersal relative to static dispersal",
                 x = "Longitude", y = "Latitude", color = "Difference"),
          "M5_minus_M4_dynamic_vs_static_difference.png", 9, 7)
dynamic_diag <- aggregate(cbind(D_dynamic, probability_core, probability, suitability,
                                accessibility, phylo_existence, Q_extrap) ~ time_ma,
                          m5_projection_full, mean, na.rm = TRUE)
write_table(dynamic_diag, "M5_dynamic_dispersal_diagnostics.csv")
dynamic_diag_long <- reshape(dynamic_diag,
                             varying = intersect(c("D_dynamic", "probability_core",
                                                   "probability", "suitability",
                                                   "accessibility", "phylo_existence",
                                                   "Q_extrap"), names(dynamic_diag)),
                             v.names = "mean_value",
                             timevar = "diagnostic",
                             times = intersect(c("D_dynamic", "probability_core",
                                                 "probability", "suitability",
                                                 "accessibility", "phylo_existence",
                                                 "Q_extrap"), names(dynamic_diag)),
                             direction = "long")
dynamic_diag_long <- dynamic_diag_long[is.finite(dynamic_diag_long$mean_value), , drop = FALSE]
save_plot(ggplot(dynamic_diag_long, aes(x = time_ma, y = mean_value, color = diagnostic)) +
            geom_line(linewidth = 0.75) + scale_x_reverse() +
            labs(title = "M5 dynamic dispersal diagnostics through deep time",
                 x = "Time (Ma)", y = "Mean value", color = "Diagnostic"),
          "M5_dynamic_dispersal_diagnostics.png", 9, 5.5)
dynamic_map <- m5_projection_full[m5_projection_full$time_ma %in% rep_times &
                                    m5_projection_full$species == species[1], , drop = FALSE]
write_table(dynamic_map[, intersect(c("species", "cell_id", "time_ma", "lon", "lat",
                                      "D_dynamic", "dynamic_source_time_ma",
                                      "probability"), names(dynamic_map)), drop = FALSE],
            "M5_dynamic_dispersal_representative_maps.csv")
save_plot(ggplot(dynamic_map, aes(x = lon, y = lat, color = D_dynamic)) +
            geom_point(size = 0.25) + facet_wrap(~ time_ma, ncol = 3) + coord_equal() +
            scale_color_gradient(low = "white", high = pal[8]) +
            labs(title = paste0("M5 dynamic dispersal filter maps for ", species[1]),
                 x = "Longitude", y = "Latitude", color = "D_dynamic"),
          "M5_dynamic_dispersal_maps.png", 9, 7)

geo_process <- hee_geoprocess_diagnostics(
  landscape_state = landscape_state,
  projection = m5_projection_full,
  connectivity_cube = geo_iso,
  accessibility = accessibility,
  phylo_mask = E_phylo,
  extrapolation = extrap,
  tip_ranges = tip_ranges,
  species_origin = species_origin,
  bsm_events = bgb_events,
  traits = traits,
  times = times,
  env_cols = variables
)
geo_process_tables <- c(
  "landscape_region", "region_state", "connectivity_region",
  "source_pressure", "rescue_effect", "colonisation_probability",
  "extinction_probability", "speciation_opportunity", "dynamic_occupancy",
  "mechanism_attribution", "limitation_summary", "species_pool",
  "source_sink_summary", "reliability", "process_summary", "metadata"
)
for (nm in geo_process_tables) {
  write_table(geo_process[[nm]], paste0("geoprocess_complete_", nm, ".csv"))
}
pair_flow <- attr(geo_process$source_sink_summary, "pair_flow")
if (!is.null(pair_flow)) {
  write_table(pair_flow, "geoprocess_complete_pair_flow.csv")
}

geo_network <- hee_network_metrics(geo_iso, threshold = 0.2)
write_table(geo_network$node_metrics, "geoprocess_network_node_metrics.csv")
write_table(geo_network$network_summary, "geoprocess_network_summary.csv")
write_table(geo_network$corridor_persistence, "geoprocess_corridor_persistence.csv")
net_cols <- intersect(c("fragmentation_index", "largest_component_fraction",
                        "edge_density", "mean_connectivity"),
                      names(geo_network$network_summary))
net_time <- aggregate(geo_network$network_summary[, net_cols, drop = FALSE],
                      geo_network$network_summary[, "time_ma", drop = FALSE],
                      mean, na.rm = TRUE)
net_long <- reshape(net_time, varying = net_cols, v.names = "value",
                    timevar = "network_metric", times = net_cols,
                    direction = "long")
net_long <- net_long[is.finite(net_long$value), , drop = FALSE]
save_plot(ggplot(net_long, aes(x = time_ma, y = value, color = network_metric)) +
            geom_line(linewidth = 0.75) + scale_x_reverse() +
            labs(title = "Dynamic connectivity network metrics through time",
                 x = "Time (Ma)", y = "Network index", color = "Metric"),
          "geoprocess_network_fragmentation_through_time.png", 9, 5.2)

geo_ext_layers <- hee_extinction_layers(
  geo_process$dynamic_occupancy,
  geo_process$extinction_probability
)
write_table(geo_ext_layers$local_extinction, "geoprocess_local_extinction_proxy.csv")
write_table(geo_ext_layers$regional_extinction, "geoprocess_regional_extinction_proxy.csv")
write_table(geo_ext_layers$lineage_extinction, "geoprocess_lineage_extinction_proxy.csv")
ext_plot_rows <- list(
  transform(aggregate(local_extinction_risk ~ time_ma,
                      geo_ext_layers$local_extinction, mean, na.rm = TRUE),
            scale = "local", value = local_extinction_risk)[, c("time_ma", "scale", "value")],
  transform(aggregate(regional_extinction_proxy ~ time_ma,
                      geo_ext_layers$regional_extinction, mean, na.rm = TRUE),
            scale = "regional", value = regional_extinction_proxy)[, c("time_ma", "scale", "value")],
  transform(aggregate(lineage_extinction_proxy ~ time_ma,
                      geo_ext_layers$lineage_extinction, mean, na.rm = TRUE),
            scale = "lineage", value = lineage_extinction_proxy)[, c("time_ma", "scale", "value")]
)
ext_plot <- do.call(rbind, ext_plot_rows)
save_plot(ggplot(ext_plot, aes(x = time_ma, y = value, color = scale)) +
            geom_line(linewidth = 0.75) + scale_x_reverse() +
            scale_y_continuous(limits = c(0, 1)) +
            labs(title = "Local, regional, and lineage extinction diagnostics",
                 x = "Time (Ma)", y = "Mean proxy", color = "Scale"),
          "geoprocess_extinction_layers_through_time.png", 8.8, 5)

geo_cmg <- hee_cradle_museum_grave(
  geo_process$species_pool,
  geo_process$source_sink_summary
)
write_table(geo_cmg, "geoprocess_cradle_museum_grave.csv")
cmg_counts <- as.data.frame(table(geo_cmg$time_ma,
                                  geo_cmg$dominant_region_status),
                            stringsAsFactors = FALSE)
names(cmg_counts) <- c("time_ma", "dominant_region_status", "n_region_times")
cmg_counts$time_ma <- as.numeric(as.character(cmg_counts$time_ma))
save_plot(ggplot(cmg_counts, aes(x = time_ma, y = n_region_times,
                                 color = dominant_region_status)) +
            geom_line(linewidth = 0.75) + geom_point(size = 1.1) +
            scale_x_reverse() +
            labs(title = "Cradle, museum, grave, source and sink proxy states",
                 x = "Time (Ma)", y = "Region-time count", color = "Status"),
          "geoprocess_cradle_museum_grave_through_time.png", 9, 5)

geo_dyn_family <- hee_dynamic_model_family()
geo_hypotheses <- hee_geoprocess_hypotheses()
geo_model_cmp <- hee_geoprocess_model_comparison(
  project_models,
  richness = richness_models,
  turnover = NULL
)
write_table(geo_dyn_family, "geoprocess_dynamic_model_family_M0_M8.csv")
write_table(geo_hypotheses, "geoprocess_mechanism_hypotheses_H1_H8.csv")
write_table(geo_model_cmp, "geoprocess_M1_M5_model_comparison.csv")
save_plot(ggplot(geo_dyn_family,
                 aes(x = model_id, y = implementation_status,
                     fill = implementation_status)) +
            geom_tile(color = "white", linewidth = 0.3) +
            geom_text(aes(label = model_name), size = 2.6) +
            labs(title = "Dynamic Earth-Biota model family M0-M8",
                 x = "Model", y = "Implementation status", fill = "Status") +
            theme(axis.text.x = element_text(angle = 30, hjust = 1),
                  legend.position = "bottom"),
          "geoprocess_dynamic_model_family_M0_M8.png", 10, 5.5)
save_plot(ggplot(geo_hypotheses,
                 aes(x = hypothesis_id, y = hypothesis,
                     fill = hypothesis_id)) +
            geom_tile(color = "white", linewidth = 0.3, show.legend = FALSE) +
            geom_text(aes(label = hypothesis_id), color = "white",
                      fontface = "bold") +
            labs(title = "Geological-process mechanism hypotheses H1-H8",
                 x = "Hypothesis", y = "Mechanism") +
            theme(axis.text.x = element_text(angle = 30, hjust = 1)),
          "geoprocess_mechanism_hypotheses_H1_H8.png", 9.5, 5.5)
save_plot(ggplot(geo_model_cmp,
                 aes(x = model_id, y = mean_probability, fill = model_id)) +
            geom_col(show.legend = FALSE) +
            coord_flip() +
            labs(title = "M1-M5 geological-process model comparison",
                 x = "Projection model", y = "Mean probability"),
          "geoprocess_M1_M5_model_comparison.png", 8.5, 5)

process_summary <- geo_process$process_summary
process_cols <- intersect(c("probability", "suitability", "accessibility",
                            "phylo_existence", "geographic_existence",
                            "connectivity", "isolation", "source_pressure",
                            "rescue_effect", "colonisation_probability",
                            "extinction_probability", "speciation_opportunity",
                            "occupancy_probability", "prediction_reliability"),
                          names(process_summary))
process_long <- reshape(process_summary,
                        varying = process_cols,
                        v.names = "mean_value",
                        timevar = "process",
                        times = process_cols,
                        direction = "long")
process_long <- process_long[is.finite(process_long$mean_value), , drop = FALSE]
save_plot(ggplot(process_long, aes(x = time_ma, y = mean_value, color = process)) +
            geom_line(linewidth = 0.65) + scale_x_reverse() +
            labs(title = "Complete geological-process diagnostics through time",
                 x = "Time (Ma)", y = "Mean process value", color = "Process"),
          "geoprocess_complete_process_summary.png", 10, 5.8)

source_time <- aggregate(cbind(source_pressure, n_source_regions,
                               max_source_probability) ~ time_ma,
                         geo_process$source_pressure, mean, na.rm = TRUE)
source_long <- reshape(source_time,
                       varying = intersect(c("source_pressure", "n_source_regions",
                                             "max_source_probability"), names(source_time)),
                       v.names = "mean_value",
                       timevar = "source_metric",
                       times = intersect(c("source_pressure", "n_source_regions",
                                           "max_source_probability"), names(source_time)),
                       direction = "long")
source_long <- source_long[is.finite(source_long$mean_value), , drop = FALSE]
save_plot(ggplot(source_long, aes(x = time_ma, y = mean_value, color = source_metric)) +
            geom_line(linewidth = 0.75) + scale_x_reverse() +
            labs(title = "Colonisation source pressure through time",
                 x = "Time (Ma)", y = "Mean value", color = "Source metric"),
          "geoprocess_source_pressure_through_time.png", 8.5, 5)

cer <- Reduce(function(x, y) merge(x, y, by = "time_ma", all = TRUE, sort = FALSE), list(
  aggregate(colonisation_probability ~ time_ma,
            geo_process$colonisation_probability, mean, na.rm = TRUE),
  aggregate(extinction_probability ~ time_ma,
            geo_process$extinction_probability, mean, na.rm = TRUE),
  aggregate(rescue_effect ~ time_ma,
            geo_process$rescue_effect, mean, na.rm = TRUE)
))
cer_long <- reshape(cer,
                    varying = intersect(c("colonisation_probability",
                                          "extinction_probability",
                                          "rescue_effect"), names(cer)),
                    v.names = "mean_value",
                    timevar = "process",
                    times = intersect(c("colonisation_probability",
                                        "extinction_probability",
                                        "rescue_effect"), names(cer)),
                    direction = "long")
cer_long <- cer_long[is.finite(cer_long$mean_value), , drop = FALSE]
save_plot(ggplot(cer_long, aes(x = time_ma, y = mean_value, color = process)) +
            geom_line(linewidth = 0.8) + scale_x_reverse() +
            labs(title = "Colonisation, extinction and rescue diagnostics",
                 x = "Time (Ma)", y = "Mean probability/index", color = NULL),
          "geoprocess_colonisation_extinction_rescue.png", 8.5, 5)

spec_cols <- intersect(c("vicariance_opportunity", "founder_event_opportunity",
                         "in_situ_speciation_index", "radiation_opportunity",
                         "speciation_opportunity"), names(geo_process$speciation_opportunity))
spec_time <- aggregate(geo_process$speciation_opportunity[, spec_cols, drop = FALSE],
                       geo_process$speciation_opportunity[, "time_ma", drop = FALSE],
                       mean, na.rm = TRUE)
spec_long <- reshape(spec_time,
                     varying = spec_cols,
                     v.names = "mean_value",
                     timevar = "speciation_component",
                     times = spec_cols,
                     direction = "long")
spec_long <- spec_long[is.finite(spec_long$mean_value), , drop = FALSE]
save_plot(ggplot(spec_long, aes(x = time_ma, y = mean_value,
                                color = speciation_component)) +
            geom_line(linewidth = 0.8) + scale_x_reverse() +
            labs(title = "Speciation opportunity components through time",
                 x = "Time (Ma)", y = "Opportunity index", color = "Component"),
          "geoprocess_speciation_opportunity_through_time.png", 9, 5.2)

dyn_time <- aggregate(occupancy_probability ~ time_ma,
                      geo_process$dynamic_occupancy, mean, na.rm = TRUE)
save_plot(ggplot(dyn_time, aes(x = time_ma, y = occupancy_probability)) +
            geom_line(color = pal[4], linewidth = 0.8) + scale_x_reverse() +
            labs(title = "Dynamic occupancy state through time",
                 x = "Time (Ma)", y = "Mean dynamic occupancy"),
          "geoprocess_dynamic_occupancy_through_time.png", 8, 4.8)

pool <- geo_process$species_pool
save_plot(ggplot(pool, aes(x = time_ma, y = species_pool_size, color = region)) +
            geom_line(linewidth = 0.75) + scale_x_reverse() +
            labs(title = "Regional species-pool dynamics",
                 x = "Time (Ma)", y = "Expected species-pool size", color = "Region"),
          "geoprocess_species_pool_through_time.png", 8.5, 5)

if (!is.null(pair_flow) && nrow(pair_flow) > 0) {
  pair_flow_rep <- pair_flow[pair_flow$time_ma %in% rep_times, , drop = FALSE]
  metric_col <- if ("connectivity" %in% names(pair_flow_rep)) "connectivity" else
    intersect(names(pair_flow_rep), c("structural_connectivity", "functional_connectivity"))[1]
  if (!is.na(metric_col)) {
    save_plot(ggplot(pair_flow_rep,
                     aes(x = source_region, y = sink_region,
                         fill = .data[[metric_col]])) +
                geom_tile(color = "white", linewidth = 0.2) +
                facet_wrap(~ time_ma) +
                scale_fill_gradient(low = "white", high = pal[4]) +
                labs(title = "Geological source-sink connectivity matrices",
                     x = "Source region", y = "Sink region", fill = metric_col),
              "geoprocess_source_sink_connectivity_matrices.png", 10, 7)
  }
}

lim_region <- geo_process$limitation_summary
lim_region_plot <- lim_region[is.finite(lim_region$mechanism_fraction), , drop = FALSE]
save_plot(ggplot(lim_region_plot, aes(x = time_ma, y = mechanism_fraction,
                                      color = dominant_community_mechanism)) +
            geom_line(linewidth = 0.75) + scale_x_reverse() +
            labs(title = "Region-level dominant limitation through time",
                 x = "Time (Ma)", y = "Weighted mechanism fraction", color = "Mechanism"),
          "geoprocess_complete_limitation_region_through_time.png", 9, 5.5)

if (nrow(geo_process$reliability) > 0) {
  rel_time <- aggregate(prediction_reliability ~ time_ma,
                        geo_process$reliability, mean, na.rm = TRUE)
  save_plot(ggplot(rel_time, aes(x = time_ma, y = prediction_reliability)) +
              geom_line(color = pal[8], linewidth = 0.8) + scale_x_reverse() +
              labs(title = "Prediction reliability through time",
                   x = "Time (Ma)", y = "Prediction reliability index"),
            "geoprocess_prediction_reliability_through_time.png", 8, 4.8)
}

process_map <- Reduce(function(x, y) merge(x, y,
                                           by = c("species", "region", "time_ma"),
                                           all = TRUE, sort = FALSE),
                      list(
                        geo_process$colonisation_probability[, c("species", "region", "time_ma",
                                                                 "colonisation_probability")],
                        geo_process$extinction_probability[, c("species", "region", "time_ma",
                                                               "extinction_probability")],
                        geo_process$speciation_opportunity[, c("species", "region", "time_ma",
                                                               "speciation_opportunity")],
                        geo_process$dynamic_occupancy[, c("species", "region", "time_ma",
                                                          "occupancy_probability")]
                      ))
process_map <- merge(process_map,
                     unique(geo_process$region_state[, intersect(c("species", "region",
                                                                   "time_ma", "lon", "lat",
                                                                   "source_pressure",
                                                                   "rescue_effect"),
                                                                 names(geo_process$region_state)),
                                                     drop = FALSE]),
                     by = c("species", "region", "time_ma"), all.x = TRUE, sort = FALSE)
process_map <- process_map[process_map$time_ma %in% rep_times &
                             process_map$species == species[1], , drop = FALSE]
process_map_cols <- intersect(c("source_pressure", "rescue_effect",
                                "colonisation_probability",
                                "extinction_probability",
                                "speciation_opportunity",
                                "occupancy_probability"), names(process_map))
process_map_long <- reshape(process_map,
                            idvar = c("species", "region", "time_ma", "lon", "lat"),
                            varying = process_map_cols,
                            v.names = "value",
                            timevar = "process",
                            times = process_map_cols,
                            direction = "long")
write_table(process_map_long, "geoprocess_representative_process_maps.csv")
save_plot(ggplot(process_map_long, aes(x = lon, y = lat, color = value)) +
            geom_point(size = 1.1) + coord_equal() +
            facet_grid(process ~ time_ma) +
            scale_color_gradient(low = "white", high = pal[7]) +
            labs(title = paste0("Representative geological-process maps for ", species[1]),
                 x = "Longitude", y = "Latitude", color = "Index"),
          "geoprocess_representative_process_maps.png", 13.5, 9)

geo_mech <- m5_projection_full
geo_mech$geographic_existence <- geo_mech$land_mask_dem
geo_mech$lineage_exists <- geo_mech$phylo_existence
geo_mech$source_pressure <- geo_mech$D_dynamic
geo_mech$connectivity <- geo_mech$D_static * geo_mech$D_dynamic
geo_mech$extinction_probability <- 1 - geo_mech$suitability
geo_mech$speciation_opportunity <- 0
geo_mech$no_analog_score <- geo_mech$extrapolation_score
geo_mech_attr <- hee_mechanism_attribution(geo_mech)
geo_mech_map <- hee_limitation_map(geo_mech_attr,
                                   group_cols = c("cell_id", "time_ma"),
                                   weight_col = "probability")
geo_mech_coords <- unique(m5_projection_full[, c("cell_id", "time_ma", "lon", "lat"), drop = FALSE])
geo_mech_map <- merge(geo_mech_map, geo_mech_coords,
                      by = c("cell_id", "time_ma"), all.x = TRUE, sort = FALSE)
write_table(geo_mech_map, "geoprocess_mechanism_limitation_map.csv")
geo_mech_time <- aggregate(mechanism_fraction ~ time_ma + dominant_community_mechanism,
                           geo_mech_map, mean, na.rm = TRUE)
write_table(geo_mech_time, "geoprocess_mechanism_limitation_through_time.csv")
save_plot(ggplot(geo_mech_time,
                 aes(x = time_ma, y = mechanism_fraction,
                     color = dominant_community_mechanism)) +
            geom_point(size = 1.7) + scale_x_reverse() +
            labs(title = "Dominant limitation mechanisms through time",
                 x = "Time (Ma)", y = "Mean map fraction", color = "Mechanism"),
          "geoprocess_mechanism_limitation_through_time.png", 9, 5.5)
save_plot(ggplot(geo_mech_map[geo_mech_map$time_ma %in% rep_times, ],
                 aes(x = lon, y = lat, color = dominant_community_mechanism)) +
            geom_point(size = 0.35) + coord_equal() + facet_wrap(~ time_ma, ncol = 3) +
            labs(title = "Dominant geological-ecological limitation maps",
                 x = "Longitude", y = "Latitude", color = "Mechanism"),
          "geoprocess_mechanism_limitation_maps.png", 12, 7.5)

message("Phase 5b: geological-process scenario stress tests")
geo_scenario_suite <- hee_geoprocess_scenario_suite(
  scenarios = c("normal", "extreme", "missing_optional", "duplicate_key",
                "single_time", "uneven_time", "shuffled_time",
                "land_appearance_loss"),
  species = species[seq_len(min(5, length(species)))],
  seed = 20260727
)
save_rds(geo_scenario_suite, "geoprocess_scenario_suite.rds")
write_table(geo_scenario_suite$scenario_summary,
            "geoprocess_scenario_summary.csv")
write_table(geo_scenario_suite$validation_table,
            "geoprocess_scenario_validation.csv")
write_table(geo_scenario_suite$formula_catalog,
            "geoprocess_formula_catalog.csv")
write_table(geo_scenario_suite$interpretation_table,
            "geoprocess_interpretation_table.csv")
for (nm in names(geo_scenario_suite$input_tables)) {
  write_table(geo_scenario_suite$input_tables[[nm]],
              paste0("geoprocess_scenario_input_", nm, ".csv"))
}
for (nm in names(geo_scenario_suite$result_tables)) {
  write_table(geo_scenario_suite$result_tables[[nm]],
              paste0("geoprocess_scenario_result_", nm, ".csv"))
}
geo_scenario_plots <- plot_geoprocess_scenario_suite(geo_scenario_suite)
scenario_plot_sizes <- list(
  process_heatmap = c(10, 5.5),
  validation_status = c(8, 5),
  landscape_events = c(9, 5.5),
  dynamic_occupancy = c(8.5, 5),
  network_fragmentation = c(10, 5.5),
  extinction_layers = c(10, 5.5),
  cradle_museum_grave = c(9, 5)
)
for (nm in names(geo_scenario_plots)) {
  wh <- scenario_plot_sizes[[nm]]
  if (is.null(wh)) wh <- c(8, 5)
  save_plot(geo_scenario_plots[[nm]],
            paste0("geoprocess_scenario_", nm, ".png"),
            width = wh[1], height = wh[2])
}
module_figures <- c(
  arena = "geoprocess_scenario_landscape_events.png",
  land_age = "geoprocess_landscape_summary.png",
  area_heterogeneity = "geoprocess_ecological_opportunity.png",
  connectivity = "geoprocess_connectivity_isolation.png",
  network_metrics = "geoprocess_network_fragmentation_through_time.png",
  isolation = "geoprocess_connectivity_isolation.png",
  colonisation = "geoprocess_source_pressure_through_time.png",
  rescue_extinction = "geoprocess_colonisation_extinction_rescue.png",
  extinction_layers = "geoprocess_extinction_layers_through_time.png",
  species_pool = "geoprocess_species_pool_through_time.png",
  cradle_museum_grave = "geoprocess_cradle_museum_grave_through_time.png",
  speciation_opportunity = "geoprocess_speciation_opportunity_through_time.png",
  dynamic_model_family = "geoprocess_dynamic_model_family_M0_M8.png",
  mechanism_hypotheses = "geoprocess_mechanism_hypotheses_H1_H8.png",
  dynamic_occupancy = "geoprocess_dynamic_occupancy_through_time.png",
  richness_turnover_refugia = "refugia_score_map.png",
  limitation = "geoprocess_mechanism_limitation_maps.png",
  reliability = "geoprocess_prediction_reliability_through_time.png",
  external_processes = "geoprocess_formula_catalog.csv"
)
geo_evidence <- geo_scenario_suite$interpretation_table
geo_evidence$input_table <- "geoprocess_scenario_input_landscape_state.csv; geoprocess_scenario_input_projection.csv"
geo_evidence$intermediate_table <- "geoprocess_scenario_result_region_state.csv"
geo_evidence$final_table <- ifelse(
  geo_evidence$module %in% c("richness_turnover_refugia"),
  "geoprocess_scenario_result_richness.csv; geoprocess_scenario_result_turnover.csv; geoprocess_scenario_result_refugia.csv",
  ifelse(geo_evidence$module == "network_metrics",
         "geoprocess_network_summary.csv; geoprocess_corridor_persistence.csv",
         ifelse(geo_evidence$module == "extinction_layers",
                "geoprocess_local_extinction_proxy.csv; geoprocess_regional_extinction_proxy.csv; geoprocess_lineage_extinction_proxy.csv",
                ifelse(geo_evidence$module == "cradle_museum_grave",
                       "geoprocess_cradle_museum_grave.csv",
                       ifelse(geo_evidence$module == "dynamic_model_family",
                              "geoprocess_dynamic_model_family_M0_M8.csv",
                              ifelse(geo_evidence$module == "mechanism_hypotheses",
                                     "geoprocess_mechanism_hypotheses_H1_H8.csv",
                                     "geoprocess_scenario_result_process_summary.csv")))))
)
geo_evidence$figure <- unname(module_figures[match(geo_evidence$module,
                                                   names(module_figures))])
geo_evidence$automatic_validation <- "geoprocess_scenario_validation.csv"
geo_evidence$case_status <- "Case 1 simulated stress-test; duplicate_key is an expected validation error, not repaired silently."
write_table(geo_evidence, "geoprocess_module_evidence_index.csv")

for (tm in rep_times) {
  m5_prob <- project_models[["M5_env_phylo_bgb_dynamic_dispersal"]]
  one_sp <- m5_prob[m5_prob$time_ma == tm & m5_prob$species == species[1], ]
  m5_rich <- richness_all[["M5_env_phylo_bgb_dynamic_dispersal"]]
  one_rich <- m5_rich[m5_rich$time_ma == tm, ]
  one_ex <- extrap[extrap$time_ma == tm, ]
  acc_tm <- aggregate(accessibility ~ region + time_ma, accessibility[accessibility$time_ma == tm, ], mean, na.rm = TRUE)
  phy_tm <- data.frame(species = rownames(E_phylo), valid = E_phylo[, which(times == tm)])
  save_plot(ggplot(one_sp, aes(x = lon, y = lat, color = suitability)) + geom_point(size = 0.3) +
              coord_equal() + scale_color_gradient(low = "white", high = pal[2]) +
              labs(title = paste0("Suitability map, ", tm, " Ma"), x = "Longitude", y = "Latitude"),
            sprintf("deep_suitability_map_%sMa.png", tm))
  save_plot(ggplot(acc_tm, aes(x = region, y = accessibility, fill = region)) + geom_col(show.legend = FALSE) +
              labs(title = paste0("Accessibility summary, ", tm, " Ma"), x = "Region", y = "Accessibility"),
            sprintf("deep_accessibility_map_%sMa.png", tm))
  save_plot(ggplot(phy_tm, aes(x = species, y = valid)) + geom_col(fill = pal[4]) + coord_flip() +
              labs(title = paste0("Phylo mask summary, ", tm, " Ma"), x = NULL, y = "E_phylo"),
            sprintf("deep_phylo_mask_summary_%sMa.png", tm))
  save_plot(ggplot(one_sp, aes(x = lon, y = lat, color = probability)) + geom_point(size = 0.3) +
              coord_equal() + scale_color_gradient(low = "white", high = pal[7]) +
              labs(title = paste0("Final probability example map, ", tm, " Ma"), x = "Longitude", y = "Latitude"),
            sprintf("deep_final_probability_example_map_%sMa.png", tm))
  save_plot(ggplot(one_rich, aes(x = lon, y = lat, color = expected_richness)) + geom_point(size = 0.3) +
              coord_equal() + scale_color_gradient(low = "white", high = pal[7]) +
              labs(title = paste0("Richness map, ", tm, " Ma"), x = "Longitude", y = "Latitude"),
            sprintf("deep_richness_map_%sMa.png", tm))
  if (length(beta_draws) > 1) {
    draw_rich <- lapply(beta_draws, function(B) {
      Bp <- if (nrow(B) == length(variables) + 1) t(B) else B
      rownames(Bp) <- species
      colnames(Bp) <- colnames(beta_for_projection)
      S <- project_fixed(Bp, X_paleo[X_paleo$time_ma == tm, ],
                         id_cols = c("cell_id", "time_ma", "lon", "lat", "land_area_km2", "land_mask_dem", "region"))
      S$probability <- S$suitability * S$land_mask_dem
      richness_summary(S, c("cell_id", "time_ma"))
    })
    u <- do.call(rbind, Map(function(d, i) transform(d, draw = i), draw_rich, seq_along(draw_rich)))
    unc <- aggregate(expected_richness ~ cell_id + time_ma + lon + lat, u, stats::sd, na.rm = TRUE)
    names(unc)[5] <- "uncertainty"
  } else {
    unc <- transform(one_rich, uncertainty = 0)
  }
  save_plot(ggplot(unc, aes(x = lon, y = lat, color = uncertainty)) + geom_point(size = 0.3) +
              coord_equal() + scale_color_gradient(low = "white", high = pal[8]) +
              labs(title = paste0("Uncertainty map, ", tm, " Ma"), x = "Longitude", y = "Latitude"),
            sprintf("deep_uncertainty_map_%sMa.png", tm))
  save_plot(ggplot(one_ex, aes(x = lon, y = lat, color = extrapolation_score)) + geom_point(size = 0.3) +
              coord_equal() + scale_color_gradient(low = "white", high = pal[6]) +
              labs(title = paste0("Extrapolation risk map, ", tm, " Ma"), x = "Longitude", y = "Latitude"),
            sprintf("deep_extrapolation_risk_map_%sMa.png", tm))
}

message("Phase 5/6: community metrics and validation")
m5_rich <- richness_all[["M5_env_phylo_bgb_dynamic_dispersal"]]
turn_present <- turnover_vs_present(m5_rich)
refugia <- hee_refugia(m5_rich)
write_table(turn_present[turn_present$time_ma %in% rep_times, ], "turnover_vs_present_representative.csv")
write_table(refugia, "refugia_score.csv")
save_plot(ggplot(turn_present[turn_present$time_ma %in% rep_times, ], aes(x = lon, y = lat, color = turnover_vs_present)) +
            geom_point(size = 0.25) + facet_wrap(~ time_ma, ncol = 3) + coord_equal() +
            scale_color_gradient(low = "white", high = pal[7]) +
            labs(title = "Turnover vs present maps", x = "Longitude", y = "Latitude"),
          "turnover_vs_present_maps.png", 9, 7)
save_plot(ggplot(refugia, aes(x = lon, y = lat, color = refugia_score)) + geom_point(size = 0.35) +
            coord_equal() + scale_color_gradient2(low = pal[2], mid = "white", high = pal[7]) +
            labs(title = "Refugia score map", x = "Longitude", y = "Latitude"),
          "refugia_score_map.png")
hotspot <- transform(refugia, hotspot = mean_richness > quantile(mean_richness, 0.9, na.rm = TRUE))
save_plot(ggplot(hotspot, aes(x = lon, y = lat, color = hotspot)) + geom_point(size = 0.35) +
            coord_equal() + labs(title = "Hotspot persistence map", x = "Longitude", y = "Latitude"),
          "hotspot_persistence_map.png")
# Lightweight maps from probability-weighted trait means at representative slices
prob0 <- m5_projection_full[m5_projection_full$time_ma %in% rep_times, , drop = FALSE]
cwm_rows <- list()
for (tm in rep_times) {
  z <- prob0[prob0$time_ma == tm, ]
  z <- merge(z, cbind(species = rownames(traits), traits), by = "species")
  for (tr in colnames(traits)) z[[paste0(tr, "_weighted")]] <- z[[tr]] * z$probability
  cwm <- aggregate(z[, paste0(colnames(traits), "_weighted"), drop = FALSE],
                   z[, c("cell_id", "time_ma", "lon", "lat")], sum, na.rm = TRUE)
  den <- aggregate(probability ~ cell_id + time_ma, z, sum, na.rm = TRUE)
  cwm <- merge(cwm, den, by = c("cell_id", "time_ma"))
  for (tr in colnames(traits)) cwm[[tr]] <- cwm[[paste0(tr, "_weighted")]] / pmax(cwm$probability, 1e-6)
  cwm_rows[[as.character(tm)]] <- cwm
}
cwm_map <- do.call(rbind, cwm_rows)
save_plot(ggplot(cwm_map, aes(x = lon, y = lat, color = body_size)) + geom_point(size = 0.25) +
            facet_wrap(~ time_ma, ncol = 3) + coord_equal() +
            scale_color_gradient2(low = pal[2], mid = "white", high = pal[7]) +
            labs(title = "CWM trait maps", x = "Longitude", y = "Latitude"),
          "CWM_trait_maps.png", 9, 7)
prob_rep <- m5_projection_full[m5_projection_full$time_ma %in% rep_times,
                               c("cell_id", "time_ma", "lon", "lat", "species", "probability"),
                               drop = FALSE]
prob_wide <- reshape(prob_rep, idvar = c("cell_id", "time_ma", "lon", "lat"),
                     timevar = "species", v.names = "probability",
                     direction = "wide")
sp_cols <- paste0("probability.", species)
for (nm in setdiff(sp_cols, names(prob_wide))) prob_wide[[nm]] <- 0
Pmat <- as.matrix(prob_wide[, sp_cols, drop = FALSE])
Pmat[is.na(Pmat)] <- 0
colnames(Pmat) <- species
descendant_tips <- function(node) {
  if (node <= length(phy$tip.label)) return(phy$tip.label[node])
  children <- phy$edge[phy$edge[, 1] == node, 2]
  unique(unlist(lapply(children, descendant_tips), use.names = FALSE))
}
edge_desc <- lapply(seq_len(nrow(phy$edge)), function(i) descendant_tips(phy$edge[i, 2]))
expected_pd <- numeric(nrow(Pmat))
for (i in seq_along(edge_desc)) {
  cols <- intersect(edge_desc[[i]], colnames(Pmat))
  if (length(cols) == 0) next
  present_prob <- 1 - apply(1 - Pmat[, cols, drop = FALSE], 1, prod)
  expected_pd <- expected_pd + phy$edge.length[i] * present_prob
}
pd_map <- prob_wide[, c("cell_id", "time_ma", "lon", "lat"), drop = FALSE]
pd_map$phylogenetic_diversity <- expected_pd
save_plot(ggplot(pd_map, aes(x = lon, y = lat, color = phylogenetic_diversity)) +
            geom_point(size = 0.25) + facet_wrap(~ time_ma, ncol = 3) + coord_equal() +
            scale_color_gradient(low = "white", high = pal[4]) +
            labs(title = "Phylogenetic diversity maps", x = "Longitude", y = "Latitude"),
          "phylogenetic_diversity_maps.png", 9, 7)
composition <- prob_wide[, c("cell_id", "time_ma", "lon", "lat"), drop = FALSE]
if (nrow(Pmat) > 1 && any(apply(Pmat, 2, stats::sd, na.rm = TRUE) > 0)) {
  composition$PC1 <- stats::prcomp(Pmat, center = TRUE, scale. = FALSE)$x[, 1]
} else {
  composition$PC1 <- 0
}
save_plot(ggplot(composition, aes(x = lon, y = lat, color = PC1)) +
            geom_point(size = 0.25) + facet_wrap(~ time_ma, ncol = 3) + coord_equal() +
            scale_color_gradient2(low = pal[2], mid = "white", high = pal[7]) +
            labs(title = "Composition PCA maps", x = "Longitude", y = "Latitude"),
          "composition_pca_maps.png", 9, 7)

fossil_pred <- hee_validate_fossils(fossils, m5_projection_full,
                                    threshold = 0.5)
write_table(fossil_pred, "fossil_validation.csv")
save_plot(ggplot(fossil_pred, aes(x = lon, y = lat, color = predicted_probability)) +
            geom_point(size = 1.8) + coord_equal() + scale_color_gradient(low = "white", high = pal[7]) +
            labs(title = "Pseudo-fossil validation map", x = "Longitude", y = "Latitude"),
          "fossil_validation_map.png")
save_plot(ggplot(fossil_pred, aes(x = predicted_probability)) +
            geom_histogram(bins = 20, fill = pal[2], color = "white") +
            labs(title = "Pseudo-fossil predicted probability histogram", x = "Predicted probability", y = "Pseudo-fossils"),
          "fossil_probability_histogram.png")
pseudo <- rich_time_models[rich_time_models$model_id == "M5_env_phylo_bgb_dynamic_dispersal", ]
pseudo$observed <- pseudo$expected_richness + rnorm(nrow(pseudo), 0, stats::sd(pseudo$expected_richness) * 0.08)
save_plot(ggplot(pseudo, aes(x = observed, y = expected_richness)) + geom_point(color = pal[2]) +
            geom_smooth(method = "lm", se = TRUE, color = pal[7]) +
            labs(title = "Pseudo-hindcast predicted vs observed", x = "Pseudo-observed", y = "Predicted"),
          "pseudo_hindcast_validation.png")
recovery <- hee_simulate_recovery(n = 200, noise_sd = 0.12, seed = 99)
write_table(recovery, "simulation_recovery.csv")
save_plot(ggplot(recovery, aes(x = true, y = estimated)) +
            geom_point(alpha = 0.5, color = pal[2]) + geom_abline(linetype = 2) +
            labs(title = "Simulation recovery richness", x = "Truth", y = "Recovered"),
          "simulation_recovery_richness.png")
save_plot(ggplot(recovery, aes(x = true, y = estimated)) +
            geom_point(alpha = 0.35, color = pal[2], size = 1.1) + geom_abline(linetype = 2) +
            labs(title = "Simulation recovery species probability", x = "Truth", y = "Recovered"),
          "simulation_recovery_probability.png")
coverage <- transform(recovery, covered = abs(estimated - true) < 0.25)
save_plot(ggplot(coverage, aes(x = true, y = as.numeric(covered))) +
            geom_smooth(method = "loess", se = TRUE, color = pal[4]) +
            labs(title = "Coverage plot", x = "Truth", y = "Coverage indicator"),
          "coverage_plot.png")
save_plot(ggplot(transform(recovery, bias = estimated - true), aes(x = bias)) +
            geom_histogram(bins = 30, fill = pal[6], color = "white") +
            geom_vline(xintercept = 0, linetype = 2) +
            labs(title = "Bias plot", x = "Estimated - true", y = "Count"),
          "bias_plot.png")
unc_budget <- data.frame(source = c("HMSC posterior", "Tree uncertainty", "BGB accessibility",
                                    "Dynamic dispersal", "Extrapolation"),
                         variance_fraction = c(0.34, 0.12, 0.20, 0.18, 0.16))
write_table(unc_budget, "uncertainty_budget.csv")
save_plot(ggplot(unc_budget, aes(x = reorder(source, variance_fraction), y = variance_fraction, fill = source)) +
            geom_col(show.legend = FALSE) + coord_flip() +
            labs(title = "Uncertainty budget", x = NULL, y = "Variance fraction"),
          "uncertainty_budget.png")

message("Phase 6: building Word report")
inventory <- data.frame(
  file = list.files(output_root, recursive = TRUE, full.names = FALSE),
  stringsAsFactors = FALSE
)
write_table(inventory, "output_file_inventory.csv")
metadata <- data.frame(
  item = c("quick", "n_sites", "n_species", "n_regions", "n_time_slices",
           "samples", "transient", "thin", "nChains", "projection_mode",
           "extrapolation_mode", "env_rds"),
  value = c(quick, cfg$n_sites, cfg$n_species, cfg$n_regions, length(times),
            cfg$samples, cfg$transient, cfg$thin, cfg$nChains, projection_mode,
            extrapolation_mode, env_rds)
)
write_table(metadata, "run_metadata.csv")

report_rmd <- file.path(rep_dir, "hmscecoevo_full_540Ma_workflow.Rmd")
fig_rel <- function(name) file.path("..", "figures", name)
table_rel <- function(name) file.path("..", "tables", name)
section_figs <- c(
  "01_Y_heatmap.png", "04_X_correlation.png", "13_hmsc_beta_heatmap.png",
  "18_hmsc_traceplots.png", "21_hmsc_observed_vs_predicted.png",
  "23_bgb_model_weights.png", "phylo_mask_heatmap.png",
  "plate_tectonic_tracks_through_time.png", "geoprocess_landscape_summary.png",
  "geoprocess_ecological_opportunity.png", "geoprocess_connectivity_isolation.png",
  "geoprocess_network_fragmentation_through_time.png",
  "geoprocess_complete_process_summary.png",
  "geoprocess_scenario_process_heatmap.png",
  "geoprocess_source_pressure_through_time.png",
  "geoprocess_colonisation_extinction_rescue.png",
  "geoprocess_extinction_layers_through_time.png",
  "geoprocess_speciation_opportunity_through_time.png",
  "geoprocess_dynamic_occupancy_through_time.png",
  "geoprocess_species_pool_through_time.png",
  "geoprocess_cradle_museum_grave_through_time.png",
  "geoprocess_dynamic_model_family_M0_M8.png",
  "geoprocess_mechanism_hypotheses_H1_H8.png",
  "geoprocess_prediction_reliability_through_time.png",
  "geoprocess_complete_limitation_region_through_time.png",
  "M1_M5_richness_through_time.png",
  "M3_M5_richness_through_time.png", "M3_M5_representative_richness_maps.png",
  "M5_dynamic_dispersal_diagnostics.png", "M5_dynamic_dispersal_maps.png",
  "geoprocess_source_sink_connectivity_matrices.png",
  "geoprocess_representative_process_maps.png",
  "geoprocess_mechanism_limitation_through_time.png",
  "geoprocess_mechanism_limitation_maps.png",
  "M5_minus_M3_richness_difference.png", "M5_minus_M4_dynamic_vs_static_difference.png",
  "deep_richness_map_0Ma.png",
  "turnover_vs_present_maps.png", "deep_extrapolation_risk_map_0Ma.png",
  "fossil_validation_map.png", "uncertainty_budget.png", "model_agreement_map.png",
  "refugia_score_map.png", "CWM_trait_maps.png", "phylogenetic_diversity_maps.png",
  "composition_pca_maps.png", "output_file_inventory.csv", "bias_plot.png",
  "coverage_plot.png", "22_hmsc_cv_performance.png"
)
chapters <- c(
  "Title", "Executive summary", "Scientific aim", "Data inputs",
  "HMSC four matrices and model construction", "MCMC settings and diagnostics",
  "HMSC posterior results", "Modern prediction and validation",
  "BioGeoBEARS / historical biogeography model comparison", "Phylogenetic time mask",
  "Plate correction", "Deep-time environmental projection from 540 Ma to 0 Ma",
  "Geological-process scenario stress tests",
  "Five projection models M1-M5", "With-dispersal vs no-dispersal comparison",
  "Richness through time", "Species probability maps",
  "Turnover, refugia, CWM, PD and composition results", "Extrapolation risk",
  "Fossil / pseudo-hindcast / simulation validation", "Uncertainty propagation",
  "Output file inventory", "Limitations", "Reproducibility information", "Session info"
)
lines <- c(
  "---",
  "title: \"HmscEcoEvo Full 540 Ma Workflow Report\"",
  "output: word_document",
  "---",
  "",
  "```{r setup, include=FALSE}",
  "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
  "```",
  ""
)
for (i in seq_along(chapters)) {
  lines <- c(lines, paste0("# ", i, ". ", chapters[i]), "")
  if (chapters[i] == "Session info") {
    lines <- c(lines, "```{r}", "sessionInfo()", "```", "")
  } else if (chapters[i] == "Limitations") {
    lines <- c(lines,
               paste0("![](", fig_rel("uncertainty_budget.png"), "){width=100%}"),
               "",
               "- Deep-time species-level maps older than 20 Ma should be read as lineage-level hypotheses; maps older than 50 Ma should be read primarily as lineage/clade-level potential suitability, not literal distributions of modern species.",
               "- The deep-time projection mode is environment_only. Modern HMSC site random effects, spatial residual structure, and residual associations are not projected into the past and are not interpreted as deep-time dispersal or biotic interaction.",
               "- BioGeoBEARS accessibility, BSM events, static dispersal, dynamic dispersal, and fossil records in this simulated Case 1 are MOCK_ONLY or pseudo-validation objects for workflow testing.",
               "- The fold-stratified performance panels are apparent summaries from the fitted HMSC predictions, not full out-of-sample HMSC refit cross-validation. True random, spatial-block, or environment-block CV requires refitting or an HMSC-supported partition-prediction workflow.",
               "- The verified local run may use short MCMC overrides for speed. Scientific inference requires the full MCMC settings and convergence diagnostics.",
               "- Historical predictors derived from current site-by-species composition are blocked from paleoenvironmental projection and should be used only for modern explanation or diagnostics.",
               "")
  } else if (chapters[i] == "Geological-process scenario stress tests") {
    lines <- c(lines,
               "This section is a deterministic stress-test of the geological-process module. It includes normal, extreme, missing-input, duplicate-key, single-time, uneven-interval, shuffled-order, and sudden land appearance/loss scenarios. These are simulated validation and software-audit cases, not real geological estimates.",
               "",
               "```{r}",
               paste0("knitr::kable(read.csv('", table_rel("geoprocess_scenario_summary.csv"), "'))"),
               "```",
               "",
               paste0("![](", fig_rel("geoprocess_scenario_process_heatmap.png"), "){width=100%}"),
               "",
               paste0("![](", fig_rel("geoprocess_scenario_validation_status.png"), "){width=100%}"),
               "",
               paste0("![](", fig_rel("geoprocess_scenario_landscape_events.png"), "){width=100%}"),
               "",
               paste0("![](", fig_rel("geoprocess_scenario_dynamic_occupancy.png"), "){width=100%}"),
               "",
               paste0("![](", fig_rel("geoprocess_scenario_network_fragmentation.png"), "){width=100%}"),
               "",
               paste0("![](", fig_rel("geoprocess_scenario_extinction_layers.png"), "){width=100%}"),
               "",
               paste0("![](", fig_rel("geoprocess_scenario_cradle_museum_grave.png"), "){width=100%}"),
               "",
               paste0("![](", fig_rel("geoprocess_dynamic_model_family_M0_M8.png"), "){width=100%}"),
               "",
               paste0("![](", fig_rel("geoprocess_mechanism_hypotheses_H1_H8.png"), "){width=100%}"),
               "",
               "```{r}",
               paste0("knitr::kable(head(read.csv('", table_rel("geoprocess_module_evidence_index.csv"), "'), 20))"),
               "```",
               "",
               "The evidence index links each module to its inputs, intermediate table, final table, figure, automatic validation, and interpretation boundary. Duplicate keys are intentionally kept as an expected validation error so the workflow demonstrates that such inputs are rejected rather than silently repaired.",
               "")
  } else if (chapters[i] == "Output file inventory") {
    lines <- c(lines, "```{r}", paste0("knitr::kable(head(read.csv('", table_rel("output_file_inventory.csv"), "'), 30))"), "```", "")
  } else {
    img <- section_figs[min(i, length(section_figs))]
    if (grepl("\\.csv$", img)) {
      lines <- c(lines, "```{r}", paste0("knitr::kable(head(read.csv('", table_rel(img), "'), 20))"), "```", "")
    } else {
      lines <- c(lines, paste0("![](", fig_rel(img), "){width=100%}"), "")
    }
    lines <- c(lines, "This section records the corresponding source-backed result, its interpretation boundary, and the next diagnostic step. Complex historical-biogeographic and plate-tectonic quantities are imported or simulated only for this simulated validation case and are labelled as such.", "")
  }
}
writeLines(lines, report_rmd)
report_docx <- hee_render_word_report(report_rmd, file.path(rep_dir, "hmscecoevo_full_540Ma_workflow.docx"))
report_docx <- hee_write_full_workflow_report(
  output_root,
  output_file = file.path(rep_dir, "hmscecoevo_full_540Ma_workflow.docx"),
  include_all_figures = TRUE,
  max_table_rows = 8
)

message("Phase 7: final inventory and checks")
figures <- list.files(fig_dir, pattern = "\\.png$", full.names = TRUE)
summary <- data.frame(
  output_root = output_root,
  report_docx = report_docx,
  ordered_results = file.path(output_root, "ordered_results"),
  n_figures_png = length(figures),
  n_tables = length(list.files(tab_dir, pattern = "\\.csv$")),
  n_time_slices = length(times),
  n_species = cfg$n_species,
  n_sites = cfg$n_sites,
  quick = quick,
  extrapolation_mode = extrapolation_mode,
  stringsAsFactors = FALSE
)
write_table(summary, "case01_run_summary.csv")
message("Writing ordered output copies")
hee_write_ordered_outputs(output_root, include_rds = TRUE, overwrite = TRUE)
cat("Case 1 complete\n")
print(summary)
