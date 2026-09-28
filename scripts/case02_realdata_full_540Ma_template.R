#!/usr/bin/env Rscript

# Case 2: real-data full 540-0 Ma workflow template.
# This script never fabricates real-data inputs. Missing required files stop
# immediately with a clear message.

arg_value <- function(key, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  pat <- paste0("^(--)?", key, "=")
  hit <- grep(pat, args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(pat, "", hit[1])
}

arg_flag <- function(key, default = FALSE) {
  args <- commandArgs(trailingOnly = TRUE)
  if (any(args %in% c(key, paste0("--", key)))) return(TRUE)
  value <- arg_value(key, NA_character_)
  if (is.na(value)) return(default)
  tolower(value) %in% c("true", "t", "1", "yes", "y")
}

cmd_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
script_path <- normalizePath(if (!is.na(cmd_file)) cmd_file else "scripts/case02_realdata_full_540Ma_template.R",
                             winslash = "/", mustWork = FALSE)
pkg_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
setwd(pkg_root)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(pkg_root, quiet = TRUE)
} else {
  library(HmscEcoEvo)
}
if (utils::packageVersion("HmscEcoEvo") < package_version("1.0.4") ||
    sum(hee_core_process_catalog()$layer == "core_process") != 6L) {
  stop("Case02 requires HmscEcoEvo >= 1.0.4 with the six-process API.",
       call. = FALSE)
}

need <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required for the real-data workflow.", call. = FALSE)
  }
}
need("ggplot2")
need("scales")
need("ape")
need("Hmsc")

`%||%` <- function(x, y) if (is.null(x)) y else x

data_root <- normalizePath(arg_value("data_root", "data_raw"), winslash = "/", mustWork = FALSE)
output_root <- normalizePath(arg_value("output", "outputs/case02_realdata_full_540Ma_template"),
                             winslash = "/", mustWork = FALSE)
quick <- arg_flag("quick", FALSE)
verbose <- arg_flag("verbose", TRUE)
chunk_size <- as.integer(arg_value("chunk_size", "5000"))
projection_mode <- arg_value("projection_mode", "environment_only")
if (!identical(projection_mode, "environment_only")) {
  stop("Only projection_mode = 'environment_only' is allowed by default. Modern spatial random effects are not projected into deep time.",
       call. = FALSE)
}

required_files <- c(
  "comm.csv",
  "sites.csv",
  "traits.csv",
  "tree.tre",
  "tip_ranges.csv",
  "env_cube_index.csv",
  "region_cube_index.csv",
  "landmask_cube_index.csv",
  "plate_points.csv",
  "fossils.csv"
)

for (nm in required_files) {
  path <- file.path(data_root, nm)
  if (!file.exists(path)) {
    stop("Missing required file: ", file.path("data_raw", nm), call. = FALSE)
  }
}

fig_dir <- file.path(output_root, "figures")
tab_dir <- file.path(output_root, "tables")
rds_dir <- file.path(output_root, "rds")
rep_dir <- file.path(output_root, "reports")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)

theme_hee <- function() {
  ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"),
                   strip.background = ggplot2::element_rect(fill = "#eef3fb"))
}

save_plot <- function(plot, name, width = 7, height = 5, dpi = 300) {
  base <- tools::file_path_sans_ext(name)
  path <- file.path(fig_dir, paste0(base, ".png"))
  ggplot2::ggsave(path, plot = plot + theme_hee(), width = width, height = height,
                  dpi = dpi, bg = "white")
  invisible(path)
}

write_table <- function(x, name) {
  path <- file.path(tab_dir, name)
  utils::write.csv(x, path, row.names = FALSE)
  invisible(path)
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

read_keyed_csv <- function(path, key = NULL) {
  x <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(key) && key %in% names(x)) {
    rownames(x) <- make.unique(as.character(x[[key]]))
    x[[key]] <- NULL
  } else if (ncol(x) > 1 && !is.numeric(x[[1]]) && !anyDuplicated(x[[1]])) {
    rownames(x) <- make.unique(as.character(x[[1]]))
    x[[1]] <- NULL
  }
  x
}

resolve_index_path <- function(index_file, label) {
  idx <- utils::read.csv(index_file, check.names = FALSE, stringsAsFactors = FALSE)
  candidate_cols <- intersect(c("cube_path", "rds_path", "path", "file"), names(idx))
  if (length(candidate_cols) == 0) {
    stop(label, " must contain one of these columns: cube_path, rds_path, path, file.",
         call. = FALSE)
  }
  p <- idx[[candidate_cols[1]]][1]
  if (!nzchar(p)) stop(label, " contains an empty cube path.", call. = FALSE)
  if (!grepl("^[A-Za-z]:|^/|^\\\\\\\\", p)) p <- file.path(dirname(index_file), p)
  normalizePath(p, winslash = "/", mustWork = TRUE)
}

extract_env_at_sites <- function(env_cube, sites, variables) {
  if (!all(c("lon", "lat") %in% names(sites))) {
    stop("sites.csv must contain lon and lat columns.", call. = FALSE)
  }
  grid0 <- hee_get_time_slice(hee_load_timecube(env_cube), 0, variables)
  grid0$cell_id <- seq_len(nrow(grid0))
  nearest <- vapply(seq_len(nrow(sites)), function(i) {
    which.min((grid0$lon - sites$lon[i])^2 + (grid0$lat - sites$lat[i])^2)
  }, integer(1))
  out <- grid0[nearest, variables, drop = FALSE]
  rownames(out) <- rownames(sites)
  out
}

message("Phase 1: reading real input files")
comm <- read_keyed_csv(file.path(data_root, "comm.csv"), key = "site_id")
Y <- as.matrix(comm)
storage.mode(Y) <- "numeric"
sites <- read_keyed_csv(file.path(data_root, "sites.csv"), key = "site_id")
traits <- read_keyed_csv(file.path(data_root, "traits.csv"), key = "species")
phy <- ape::read.tree(file.path(data_root, "tree.tre"))
tip_ranges <- read_keyed_csv(file.path(data_root, "tip_ranges.csv"), key = "species")
plate_points <- utils::read.csv(file.path(data_root, "plate_points.csv"), stringsAsFactors = FALSE)
fossils <- utils::read.csv(file.path(data_root, "fossils.csv"), stringsAsFactors = FALSE)

if (is.null(rownames(Y)) || any(!rownames(Y) %in% rownames(sites))) {
  stop("comm.csv site rows must match sites.csv site_id values.", call. = FALSE)
}
sites <- sites[rownames(Y), , drop = FALSE]

species <- colnames(Y)
missing_traits <- setdiff(species, rownames(traits))
if (length(missing_traits) > 0) {
  stop("traits.csv is missing species: ", paste(missing_traits, collapse = ", "), call. = FALSE)
}
traits <- traits[species, , drop = FALSE]

missing_tips <- setdiff(species, phy$tip.label)
if (length(missing_tips) > 0) {
  stop("tree.tre is missing species: ", paste(missing_tips, collapse = ", "), call. = FALSE)
}
phy <- ape::keep.tip(phy, species)
tip_ranges <- tip_ranges[species, , drop = FALSE]

env_path <- resolve_index_path(file.path(data_root, "env_cube_index.csv"), "env_cube_index.csv")
region_path <- resolve_index_path(file.path(data_root, "region_cube_index.csv"), "region_cube_index.csv")
landmask_path <- resolve_index_path(file.path(data_root, "landmask_cube_index.csv"), "landmask_cube_index.csv")
env_cube <- hee_load_timecube(env_path)
region_cube <- hee_load_timecube(region_path)
landmask_cube <- hee_load_timecube(landmask_path)
env_index <- utils::read.csv(file.path(data_root, "env_cube_index.csv"), stringsAsFactors = FALSE)
env_variables <- if ("variable" %in% names(env_index)) unique(env_index$variable) else env_cube$variables
env_variables <- setdiff(env_variables, c("land_mask_dem", "region", "region_id"))
if (length(env_variables) < 2) {
  stop("env_cube_index.csv must identify at least two environmental variables.", call. = FALSE)
}

times <- hee_time_axis_from_env_cube(env_cube)
hee_assert_full_time_axis(times, quick = quick, min_time_slices = 3,
                          required_range = c(540, 0))
rep_times_target <- unique(c(
  540, 400, 252, 201, 145, 66, 23, 2.58, 0,
  seq(540, 0, by = -60)
))
rep_times <- sort(unique(vapply(rep_times_target, function(t) {
  times[which.min(abs(times - t))]
}, numeric(1))), decreasing = TRUE)
if (!all(c(540, 0) %in% times)) {
  warning("The supplied env cube does not include both 540 Ma and 0 Ma exactly; nearest available ages will be used.")
}

XData_raw <- extract_env_at_sites(env_path, sites, env_variables)
XFormula <- stats::as.formula(paste("~", paste(env_variables, collapse = " + ")))
predictor_roles <- hee_tag_predictors(
  env_variables,
  stats::setNames(rep("deep_time_exogenous", length(env_variables)), env_variables)
)
invisible(hee_detect_leakage(XFormula, predictor_roles))
TrFormula <- if (ncol(traits) > 0) stats::as.formula(paste("~", paste(colnames(traits), collapse = " + "))) else ~ 1
recipe <- hee_lock_recipe(XData_raw, XFormula)
XData <- hee_apply_recipe(recipe, XData_raw)
C <- ape::vcv(phy, corr = TRUE)
studyDesign <- data.frame(site = rownames(Y), row.names = rownames(Y))
ranLevels <- list(site = Hmsc::HmscRandomLevel(units = studyDesign$site))

write_table(cbind(site_id = rownames(Y), as.data.frame(Y)), "Y_matrix.csv")
write_table(cbind(site_id = rownames(XData_raw), XData_raw), "XData_0Ma_from_env_cube_raw.csv")
write_table(cbind(site_id = rownames(XData), XData), "XData_0Ma_from_env_cube.csv")
write_table(cbind(species = rownames(traits), traits), "TrData_traits.csv")
write_table(cbind(species = rownames(C), as.data.frame(C, check.names = FALSE)), "phylogenetic_correlation_C.csv")
write_table(cbind(species = rownames(tip_ranges), tip_ranges), "tip_ranges.csv")
write_table(data.frame(time_ma = times), "full_time_axis_540_0Ma.csv")
saveRDS(list(Y = Y, XData = XData, XData_raw = XData_raw, recipe = recipe,
             TrData = traits, phy = phy, C = C, sites = sites,
             tip_ranges = tip_ranges, env_cube = env_path, region_cube = region_path,
             landmask_cube = landmask_path),
        file.path(rds_dir, "case02_input_objects.rds"))

message("Phase 1: plotting real-data input matrices")
Y_long <- as.data.frame(as.table(Y))
names(Y_long) <- c("site", "species", "presence")
save_plot(ggplot2::ggplot(Y_long, ggplot2::aes(species, site, fill = presence)) +
            ggplot2::geom_tile() + ggplot2::scale_fill_gradient(low = "white", high = "#0b3b75") +
            ggplot2::labs(title = "Y community matrix", x = "Species", y = "Sites"),
          "01_Y_heatmap.png", 8, 6)
prev <- data.frame(species = species, prevalence = colMeans(Y > 0, na.rm = TRUE))
save_plot(ggplot2::ggplot(prev, ggplot2::aes(stats::reorder(species, prevalence), prevalence)) +
            ggplot2::geom_col(fill = "#1f6f8b") + ggplot2::coord_flip() +
            ggplot2::labs(title = "Species prevalence", x = NULL, y = "Prevalence"),
          "02_species_prevalence.png", 6, 6)
rich <- data.frame(site = rownames(Y), richness = rowSums(Y > 0, na.rm = TRUE))
save_plot(ggplot2::ggplot(rich, ggplot2::aes(richness)) +
            ggplot2::geom_histogram(bins = 20, fill = "#75b843", color = "white") +
            ggplot2::labs(title = "Site richness", x = "Observed richness", y = "Sites"),
          "03_site_richness.png")
xcor <- stats::cor(XData, use = "pairwise.complete.obs")
xlong <- as.data.frame(as.table(xcor)); names(xlong) <- c("var1", "var2", "correlation")
save_plot(ggplot2::ggplot(xlong, ggplot2::aes(var1, var2, fill = correlation)) +
            ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(low = "#1f62d0", mid = "white", high = "#cf232b") +
            ggplot2::labs(title = "XData correlation", x = NULL, y = NULL),
          "04_X_correlation.png")
trait_long <- as.data.frame(as.table(as.matrix(traits)))
names(trait_long) <- c("species", "trait", "value")
save_plot(ggplot2::ggplot(trait_long, ggplot2::aes(trait, species, fill = value)) +
            ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(low = "#1f62d0", mid = "white", high = "#cf232b") +
            ggplot2::labs(title = "Trait matrix", x = NULL, y = NULL),
          "06_traits_heatmap.png", 6, 7)
save_plot(ggplot2::ggplot(as.data.frame(as.table(C)), ggplot2::aes(Var1, Var2, fill = Freq)) +
            ggplot2::geom_tile() + ggplot2::scale_fill_gradient(low = "white", high = "#0b3b75") +
            ggplot2::labs(title = "Phylogenetic correlation matrix", x = NULL, y = NULL),
          "09_phylo_correlation_matrix.png", 7, 6)
png(file.path(fig_dir, "08_phylo_tree.png"), width = 1800, height = 1400, res = 220)
plot(phy, cex = 0.7, main = "Time-calibrated phylogeny")
axisPhylo()
dev.off()

message("Phase 2: fitting HMSC MCMC on real data")
samples <- as.integer(arg_value("samples", if (quick) "50" else "1000"))
transient <- as.integer(arg_value("transient", if (quick) "50" else "1000"))
thin <- as.integer(arg_value("thin", if (quick) "1" else "10"))
nChains <- as.integer(arg_value("nChains", if (quick) "2" else "4"))
nParallel <- as.integer(arg_value("nParallel", "1"))

m <- Hmsc::Hmsc(
  Y = Y,
  XData = XData,
  XFormula = XFormula,
  TrData = traits,
  TrFormula = TrFormula,
  phyloTree = phy,
  distr = "probit",
  studyDesign = studyDesign,
  ranLevels = ranLevels
)

m <- Hmsc::sampleMcmc(
  m,
  samples = samples,
  transient = transient,
  thin = thin,
  nChains = nChains,
  nParallel = nParallel,
  verbose = verbose
)
saveRDS(m, file.path(rds_dir, "hmsc_model_sampled.rds"))
diag <- hee_hmsc_diagnostics(m)
saveRDS(diag, file.path(rds_dir, "hmsc_diagnostics.rds"))

beta <- diag$Beta$mean
gamma <- diag$Gamma$mean
omega <- diag$Omega$mean
rho <- diag$Rho
write_table(cbind(axis = rownames(beta), as.data.frame(beta, check.names = FALSE)), "hmsc_beta_posterior_mean.csv")
write_table(cbind(trait = rownames(gamma), as.data.frame(gamma, check.names = FALSE)), "hmsc_gamma_posterior_mean.csv")
write_table(cbind(row = rownames(omega), as.data.frame(omega, check.names = FALSE)), "hmsc_omega_posterior_mean.csv")
write_table(rho, "hmsc_rho_posterior.csv")
if (!is.null(diag$ESS)) write_table(diag$ESS, "hmsc_ess.csv")
if (!is.null(diag$Rhat)) write_table(diag$Rhat, "hmsc_rhat.csv")

beta_long <- as.data.frame(as.table(beta)); names(beta_long) <- c("axis", "species", "value")
save_plot(ggplot2::ggplot(beta_long, ggplot2::aes(species, axis, fill = value)) +
            ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(low = "#1f62d0", mid = "white", high = "#cf232b") +
            ggplot2::labs(title = "HMSC Beta posterior mean", x = "Species", y = "Predictor"),
          "13_hmsc_beta_heatmap.png", 8, 5)
gamma_long <- as.data.frame(as.table(gamma)); names(gamma_long) <- c("trait", "axis", "value")
save_plot(ggplot2::ggplot(gamma_long, ggplot2::aes(axis, trait, fill = value)) +
            ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(low = "#1f62d0", mid = "white", high = "#cf232b") +
            ggplot2::labs(title = "HMSC Gamma posterior mean", x = "Beta axis", y = "Trait"),
          "14_hmsc_gamma_heatmap.png", 7, 5)
if (nrow(rho) > 0) {
  save_plot(ggplot2::ggplot(rho, ggplot2::aes(parameter, mean, ymin = lwr, ymax = upr)) +
              ggplot2::geom_pointrange(color = "#7f3b8c") + ggplot2::coord_flip() +
              ggplot2::labs(title = "HMSC Rho posterior", x = NULL, y = "rho"),
            "15_hmsc_rho_posterior.png")
}
omega_long <- as.data.frame(as.table(omega)); names(omega_long) <- c("species1", "species2", "value")
save_plot(ggplot2::ggplot(omega_long, ggplot2::aes(species1, species2, fill = value)) +
            ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(low = "#1f62d0", mid = "white", high = "#cf232b") +
            ggplot2::labs(title = "Omega residual association", x = NULL, y = NULL),
          "16_hmsc_omega_heatmap.png", 7, 6)

message("Phase 2: modern prediction and validation")
pred_modern <- Hmsc::computePredictedValues(m)
fit <- Hmsc::evaluateModelFit(hM = m, predY = pred_modern)
write_table(as.data.frame(fit), "hmsc_model_fit_metrics.csv")
pred_mean <- if (length(dim(pred_modern)) == 3) apply(pred_modern, c(1, 2), mean, na.rm = TRUE) else pred_modern
obs_pred <- data.frame(observed = as.vector(Y), predicted = as.vector(pred_mean))
save_plot(ggplot2::ggplot(obs_pred, ggplot2::aes(observed, predicted)) +
            ggplot2::geom_jitter(width = 0.05, height = 0, alpha = 0.25, color = "#0b3b75") +
            ggplot2::geom_smooth(method = "loess", se = TRUE, color = "#cf232b") +
            ggplot2::labs(title = "Observed vs predicted occurrence", x = "Observed", y = "Predicted probability"),
          "21_hmsc_observed_vs_predicted.png")

message("Phase 4: requiring real BioGeoBEARS accessibility")
bgb_models_file <- file.path(data_root, "bgb_model_compare.csv")
bgb_access_file <- file.path(data_root, "bgb_accessibility.csv")
bgb_events_file <- file.path(data_root, "bgb_bsm_events.csv")
static_dispersal_file <- file.path(data_root, "static_dispersal.csv")
dynamic_dispersal_file <- file.path(data_root, "dynamic_dispersal.csv")
if (!file.exists(bgb_models_file) || !file.exists(bgb_access_file)) {
  stop("Real-data workflow requires precomputed BioGeoBEARS outputs or a custom runner. Missing optional real-data file(s): data_raw/bgb_model_compare.csv and/or data_raw/bgb_accessibility.csv. No mock BioGeoBEARS result is generated for Case 2.",
       call. = FALSE)
}
if (!file.exists(static_dispersal_file) || !file.exists(dynamic_dispersal_file)) {
  stop("Real-data workflow requires precomputed dispersal filters for M4/M5. Missing optional real-data file(s): data_raw/static_dispersal.csv and/or data_raw/dynamic_dispersal.csv. No constant or mock dispersal filter is generated for Case 2.",
       call. = FALSE)
}
bgb_models <- utils::read.csv(bgb_models_file, stringsAsFactors = FALSE)
bgb_cmp <- hee_bgb_compare(bgb_models, n = nrow(tip_ranges))
bgb_access <- hee_bgb_accessibility(
  utils::read.csv(bgb_access_file, stringsAsFactors = FALSE)
)
bgb_events <- if (file.exists(bgb_events_file)) utils::read.csv(bgb_events_file, stringsAsFactors = FALSE) else data.frame()
static_filter <- utils::read.csv(static_dispersal_file, stringsAsFactors = FALSE)
dynamic_filter <- utils::read.csv(dynamic_dispersal_file, stringsAsFactors = FALSE)
normalize_filter <- function(x, label, target = c("static_weight", "dynamic_weight")) {
  target <- match.arg(target)
  if (!target %in% names(x)) {
    candidates <- intersect(c(target, "dispersal_weight", "weight", "D_static", "D_dynamic"), names(x))
    if (length(candidates) == 0) {
      stop(label, " must contain ", target, ", dispersal_weight, weight, D_static, or D_dynamic.",
           call. = FALSE)
    }
    x[[target]] <- x[[candidates[1]]]
  }
  x
}
static_filter <- normalize_filter(static_filter, "static_dispersal.csv", "static_weight")
dynamic_filter <- normalize_filter(dynamic_filter, "dynamic_dispersal.csv", "dynamic_weight")
write_table(bgb_cmp, "bgb_model_compare.csv")
if (nrow(bgb_events) > 0) {
  bsm_metrics <- hee_bgb_bsm_metrics(bgb_events)
  write_table(bsm_metrics$region_metrics, "bgb_bsm_region_metrics.csv")
  write_table(bsm_metrics$dispersal_rates, "bgb_bsm_dispersal_rates.csv")
  write_table(bsm_metrics$source_sink_long, "bgb_bsm_source_sink_turnover_long.csv")
}
save_plot(ggplot2::ggplot(bgb_cmp, ggplot2::aes(stats::reorder(model, AICc_weight), AICc_weight, fill = model)) +
            ggplot2::geom_col(show.legend = FALSE) + ggplot2::coord_flip() +
            ggplot2::labs(title = "BioGeoBEARS model weights", x = NULL, y = "AICc weight"),
          "23_bgb_model_weights.png")

message("Phase 3/5: projecting five models through 540-0 Ma")
origin <- if ("origin_ma" %in% names(tip_ranges)) stats::setNames(tip_ranges$origin_ma, rownames(tip_ranges)) else NULL
if (is.null(origin)) {
  branching <- ape::branching.times(phy)
  origin <- stats::setNames(rep(max(times), length(species)), species)
  origin[names(origin)] <- pmin(max(times), max(branching, na.rm = TRUE))
}
E_phylo <- hee_phylo_time_mask(phy, times = times, species_origin = origin, species = species,
                               level = "species")
write_table(data.frame(species = rownames(E_phylo), E_phylo, check.names = FALSE), "phylo_time_mask.csv")

beta_species <- t(beta)
common_axes <- intersect(colnames(beta_species), c("(Intercept)", env_variables))
beta_species <- beta_species[, common_axes, drop = FALSE]
suit <- hee_project_hmsc_timecube(beta = beta_species, env_cube = env_cube, recipe = recipe,
                                  times = times, variables = env_variables, link = "probit",
                                  land_only = TRUE)
region_var <- intersect(c("region", "region_id", "biogeo_region", "realm",
                          "paleo_region"), region_cube$variables)[1]
if (is.na(region_var) || !nzchar(region_var)) {
  stop("region_cube must contain a region variable such as region, region_id, biogeo_region, realm, or paleo_region.",
       call. = FALSE)
}
region_grid <- hee_make_paleo_grid(region_cube, times = times,
                                   variables = region_var, land_only = FALSE)
if (!"cell_id" %in% names(region_grid)) {
  region_grid$cell_id <- paste0("lon", sprintf("%+.3f", region_grid$lon),
                                "_lat", sprintf("%+.3f", region_grid$lat))
}
names(region_grid)[names(region_grid) == region_var] <- "region"
region_grid$region <- as.character(region_grid$region)
region_lookup <- region_grid[, c("cell_id", "time_ma", "region"), drop = FALSE]
suit <- safe_left_join(suit, region_lookup, by = c("cell_id", "time_ma"),
                       label = "region_cube")
if (any(is.na(suit$region))) {
  stop("Some projected cells could not be matched to region_cube. Check that env, region, and landmask cubes share the same grid and time axis.",
       call. = FALSE)
}
landscape_state <- hee_make_paleo_grid(env_cube, times = times,
                                       variables = env_variables,
                                       land_only = TRUE)
landscape_state <- safe_left_join(landscape_state, region_lookup,
                                  by = c("cell_id", "time_ma"),
                                  label = "region_cube")
if (any(is.na(landscape_state$region))) {
  stop("Some landscape cells could not be matched to region_cube for geological-process diagnostics.",
       call. = FALSE)
}
landmask <- unique(suit[, intersect(c("cell_id", "time_ma", "land_mask_dem"), names(suit)), drop = FALSE])
if ("land_mask_dem" %in% names(landmask)) {
  names(landmask)[names(landmask) == "land_mask_dem"] <- "land_weight"
}
projection_defs <- hee_projection_formula(include_q = FALSE)
write_table(projection_defs, "projection_model_definitions.csv")
projection_example <- hee_projection_component_example(
  S_HMSC = 0.80, A_BGB = 0.25, E_phylo = 1, L_land = 1,
  D_static = 0.60, D_dynamic = 0.10, Q_extrap = 1,
  extrapolation_mode = "flag"
)
write_table(projection_example, "projection_numeric_example.csv")
save_plot(ggplot2::ggplot(projection_example,
                          ggplot2::aes(x = model_id, y = probability, fill = model_id)) +
            ggplot2::geom_col(show.legend = FALSE) + ggplot2::coord_flip() +
            ggplot2::labs(title = "Numeric example: layered deep-time projection",
                          x = NULL, y = "Final probability P"),
          "projection_numeric_example.png", 7, 4.5)
M1 <- hee_combine(suit, landmask = landmask)
M2 <- hee_combine(suit, phylo_mask = E_phylo, landmask = landmask)
M3 <- hee_combine(suit, phylo_mask = E_phylo, accessibility = bgb_access, landmask = landmask)
M4 <- hee_combine(suit, phylo_mask = E_phylo, accessibility = bgb_access,
                  landmask = landmask, static_filter = static_filter)
M5 <- hee_combine(suit, phylo_mask = E_phylo, accessibility = bgb_access,
                  landmask = landmask, dynamic_filter = dynamic_filter)
models <- list(M1_env_only = M1, M2_env_phylo = M2,
               M3_env_phylo_bgb_no_dispersal = M3,
               M4_env_phylo_bgb_static_dispersal = M4,
               M5_env_phylo_bgb_dynamic_dispersal = M5)
richness <- do.call(rbind, lapply(names(models), function(nm) {
  z <- hee_richness(models[[nm]])
  z$model_id <- nm
  z
}))
write_table(richness, "M1_M5_richness_through_time_cell.csv")
rtt <- stats::aggregate(expected_richness ~ model_id + time_ma, richness, mean, na.rm = TRUE)
write_table(rtt, "M1_M5_richness_through_time.csv")
save_plot(ggplot2::ggplot(rtt, ggplot2::aes(time_ma, expected_richness, color = model_id)) +
            ggplot2::geom_line(linewidth = 0.8) + ggplot2::scale_x_reverse() +
            ggplot2::labs(title = "M1-M5 richness through time", x = "Time (Ma)", y = "Mean expected richness"),
          "M1_M5_richness_through_time.png", 8, 5)
selected_deep_models <- c("M3_env_phylo_bgb_no_dispersal",
                          "M4_env_phylo_bgb_static_dispersal",
                          "M5_env_phylo_bgb_dynamic_dispersal")
model_short <- function(x) {
  out <- gsub("M3_env_phylo_bgb_no_dispersal", "M3 no dispersal", x, fixed = TRUE)
  out <- gsub("M4_env_phylo_bgb_static_dispersal", "M4 static dispersal", out, fixed = TRUE)
  out <- gsub("M5_env_phylo_bgb_dynamic_dispersal", "M5 dynamic dispersal", out, fixed = TRUE)
  out
}
m3_m5_rtt <- rtt[rtt$model_id %in% selected_deep_models, , drop = FALSE]
m3_m5_rtt$model_short <- model_short(m3_m5_rtt$model_id)
write_table(m3_m5_rtt, "M3_M5_richness_through_time.csv")
save_plot(ggplot2::ggplot(m3_m5_rtt, ggplot2::aes(time_ma, expected_richness, color = model_short)) +
            ggplot2::geom_line(linewidth = 0.8) + ggplot2::scale_x_reverse() +
            ggplot2::labs(title = "M3-M5 richness through deep time",
                          x = "Time (Ma)", y = "Mean expected richness", color = "Projection model"),
          "M3_M5_richness_through_time.png", 8, 5)
component_tables <- lapply(names(models), function(nm) {
  z <- models[[nm]]
  cols <- intersect(c("suitability", "accessibility", "phylo_existence",
                      "land_weight", "static_weight", "dynamic_weight",
                      "D_static", "D_dynamic", "Q_extrap", "probability"),
                    names(z))
  out <- stats::aggregate(z[, cols, drop = FALSE], z[, "time_ma", drop = FALSE], mean, na.rm = TRUE)
  out$model_id <- nm
  out
})
component_common <- Reduce(intersect, lapply(component_tables, names))
components <- do.call(rbind, lapply(component_tables, function(z) z[, component_common, drop = FALSE]))
write_table(components, "M1_M5_projection_component_diagnostics.csv")
m3_m5_comp <- components[components$model_id %in% selected_deep_models, , drop = FALSE]
m3_m5_comp$model_short <- model_short(m3_m5_comp$model_id)
component_cols <- setdiff(names(m3_m5_comp), c("time_ma", "model_id", "model_short"))
m3_m5_comp_long <- reshape(m3_m5_comp,
                           varying = component_cols,
                           v.names = "mean_value",
                           timevar = "component",
                           times = component_cols,
                           direction = "long")
m3_m5_comp_long <- m3_m5_comp_long[is.finite(m3_m5_comp_long$mean_value), , drop = FALSE]
write_table(m3_m5_comp_long, "M3_M5_component_diagnostics_through_time.csv")
save_plot(ggplot2::ggplot(m3_m5_comp_long, ggplot2::aes(time_ma, mean_value, color = component)) +
            ggplot2::geom_line(linewidth = 0.65) + ggplot2::scale_x_reverse() +
            ggplot2::facet_wrap(~ model_short, ncol = 1) +
            ggplot2::labs(title = "M3-M5 component diagnostics through deep time",
                          x = "Time (Ma)", y = "Mean component value", color = "Component"),
          "M3_M5_component_diagnostics_through_time.png", 8, 7)
m3_m5_rep <- richness[richness$model_id %in% selected_deep_models &
                        richness$time_ma %in% rep_times, , drop = FALSE]
m3_m5_rep$model_short <- model_short(m3_m5_rep$model_id)
write_table(m3_m5_rep, "M3_M5_representative_richness_maps.csv")
if (all(c("lon", "lat") %in% names(m3_m5_rep))) {
  save_plot(ggplot2::ggplot(m3_m5_rep, ggplot2::aes(lon, lat, color = expected_richness)) +
              ggplot2::geom_point(size = 0.2) + ggplot2::facet_grid(model_short ~ time_ma) +
              ggplot2::coord_equal() + ggplot2::scale_color_gradient(low = "white", high = "#cf232b") +
              ggplot2::labs(title = "M3-M5 representative richness maps",
                            x = "Longitude", y = "Latitude", color = "Richness"),
            "M3_M5_representative_richness_maps.png", 12, 7)
}
dynamic_diag <- components[components$model_id == "M5_env_phylo_bgb_dynamic_dispersal", , drop = FALSE]
write_table(dynamic_diag, "M5_dynamic_dispersal_diagnostics.csv")
dynamic_cols <- setdiff(names(dynamic_diag), c("time_ma", "model_id"))
dynamic_long <- reshape(dynamic_diag,
                        varying = dynamic_cols,
                        v.names = "mean_value",
                        timevar = "diagnostic",
                        times = dynamic_cols,
                        direction = "long")
dynamic_long <- dynamic_long[is.finite(dynamic_long$mean_value), , drop = FALSE]
save_plot(ggplot2::ggplot(dynamic_long, ggplot2::aes(time_ma, mean_value, color = diagnostic)) +
            ggplot2::geom_line(linewidth = 0.75) + ggplot2::scale_x_reverse() +
            ggplot2::labs(title = "M5 dynamic dispersal diagnostics through deep time",
                          x = "Time (Ma)", y = "Mean value", color = "Diagnostic"),
          "M5_dynamic_dispersal_diagnostics.png", 8, 5)

message("Phase 5: geological-process diagnostics")
landscape_state <- hee_land_age(landscape_state, unit_col = "cell_id",
                                land_col = "land_mask_dem")
landscape_events <- hee_landscape_events(landscape_state, unit_col = "cell_id",
                                         land_col = "land_mask_dem",
                                         area_col = "land_area_km2")
write_table(landscape_events, "geoprocess_landscape_events.csv")
landscape_summary <- stats::aggregate(cbind(geographic_existence,
                                            land_age_myr) ~ time_ma,
                                      landscape_state, mean, na.rm = TRUE)
write_table(landscape_summary, "geoprocess_landscape_summary.csv")
save_plot(ggplot2::ggplot(landscape_summary, ggplot2::aes(time_ma)) +
            ggplot2::geom_line(ggplot2::aes(y = geographic_existence,
                                            color = "land fraction"),
                               linewidth = 0.8) +
            ggplot2::geom_line(ggplot2::aes(y = scales::rescale(land_age_myr,
                                                                 to = c(0, 1),
                                                                 from = range(land_age_myr, na.rm = TRUE)),
                                            color = "mean land age (scaled)"),
                               linewidth = 0.8) +
            ggplot2::scale_x_reverse() +
            ggplot2::labs(title = "Geological arena through time",
                          x = "Time (Ma)", y = "Scaled value", color = NULL),
          "geoprocess_landscape_summary.png", 8, 5)
eco_opp <- hee_ecological_opportunity(landscape_state, env_cols = env_variables,
                                      region_col = "region",
                                      area_col = "land_area_km2")
eco_opp_summary <- stats::aggregate(cbind(ecological_opportunity,
                                          area_gain_index,
                                          heterogeneity_index,
                                          novelty_index,
                                          young_land_index) ~ time_ma,
                                    eco_opp, mean, na.rm = TRUE)
write_table(eco_opp_summary, "geoprocess_ecological_opportunity_summary.csv")
eco_opp_long <- stats::reshape(eco_opp_summary,
                               varying = c("ecological_opportunity",
                                           "area_gain_index",
                                           "heterogeneity_index",
                                           "novelty_index",
                                           "young_land_index"),
                               v.names = "value",
                               timevar = "component",
                               times = c("opportunity", "area gain",
                                         "heterogeneity", "novelty",
                                         "young land"),
                               direction = "long")
save_plot(ggplot2::ggplot(eco_opp_long,
                          ggplot2::aes(time_ma, value, color = component)) +
            ggplot2::geom_line(linewidth = 0.75) +
            ggplot2::scale_x_reverse() +
            ggplot2::labs(title = "Ecological opportunity components through time",
                          x = "Time (Ma)", y = "Index", color = NULL),
          "geoprocess_ecological_opportunity.png", 8, 5)
region_centroids <- stats::aggregate(cbind(lon, lat, land_area_km2) ~ region + time_ma,
                                     landscape_state, mean, na.rm = TRUE)
geo_distances <- do.call(rbind, lapply(split(region_centroids,
                                             region_centroids$time_ma),
                                       function(z) {
  pair <- expand.grid(from_region = z$region, to_region = z$region,
                      stringsAsFactors = FALSE)
  a <- z[match(pair$from_region, z$region), ]
  b <- z[match(pair$to_region, z$region), ]
  pair$time_ma <- z$time_ma[1]
  pair$paleodistance <- sqrt((a$lon - b$lon)^2 + (a$lat - b$lat)^2)
  md <- max(pair$paleodistance, na.rm = TRUE)
  pair$barrier_strength <- if (is.finite(md) && md > 0) {
    pmin(pair$paleodistance / md, 1)
  } else {
    0
  }
  pair
}))
geo_conn <- hee_build_connectivity_cube(geo_distances, distance_scale = 60)
geo_iso <- hee_isolation_history(geo_conn, threshold = 0.2)
geo_conn_summary <- stats::aggregate(cbind(structural_connectivity,
                                           connectivity,
                                           isolation,
                                           isolation_duration_ma) ~ time_ma,
                                     geo_iso, mean, na.rm = TRUE)
write_table(geo_conn_summary, "geoprocess_connectivity_isolation_summary.csv")
save_plot(ggplot2::ggplot(geo_conn_summary, ggplot2::aes(time_ma)) +
            ggplot2::geom_line(ggplot2::aes(y = connectivity,
                                            color = "connectivity"),
                               linewidth = 0.8) +
            ggplot2::geom_line(ggplot2::aes(y = isolation,
                                            color = "isolation"),
                               linewidth = 0.8) +
            ggplot2::scale_x_reverse() +
            ggplot2::labs(title = "Connectivity and isolation through time",
                          x = "Time (Ma)", y = "Index", color = NULL),
          "geoprocess_connectivity_isolation.png", 8, 5)
extrap <- hee_extrapolation_risk(XData_raw, env_cube, variables = env_variables,
                                 times = times, method = "range")
geo_process <- hee_geoprocess_diagnostics(
  landscape_state = landscape_state,
  projection = M5,
  connectivity_cube = geo_iso,
  accessibility = bgb_access,
  phylo_mask = E_phylo,
  extrapolation = extrap,
  tip_ranges = tip_ranges,
  bsm_events = bgb_events,
  traits = traits,
  times = times,
  env_cols = env_variables
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
if (!is.null(pair_flow)) write_table(pair_flow, "geoprocess_complete_pair_flow.csv")

geo_network <- hee_network_metrics(geo_iso, threshold = 0.2)
write_table(geo_network$node_metrics, "geoprocess_network_node_metrics.csv")
write_table(geo_network$network_summary, "geoprocess_network_summary.csv")
write_table(geo_network$corridor_persistence, "geoprocess_corridor_persistence.csv")
net_cols <- intersect(c("fragmentation_index", "largest_component_fraction",
                        "edge_density", "mean_connectivity"),
                      names(geo_network$network_summary))
net_time <- stats::aggregate(geo_network$network_summary[, net_cols, drop = FALSE],
                             geo_network$network_summary[, "time_ma", drop = FALSE],
                             mean, na.rm = TRUE)
net_long <- stats::reshape(net_time, varying = net_cols, v.names = "value",
                           timevar = "network_metric", times = net_cols,
                           direction = "long")
net_long <- net_long[is.finite(net_long$value), , drop = FALSE]
save_plot(ggplot2::ggplot(net_long,
                          ggplot2::aes(time_ma, value,
                                       color = network_metric)) +
            ggplot2::geom_line(linewidth = 0.75) +
            ggplot2::scale_x_reverse() +
            ggplot2::labs(title = "Dynamic connectivity network metrics through time",
                          x = "Time (Ma)", y = "Network index",
                          color = "Metric"),
          "geoprocess_network_fragmentation_through_time.png", 9, 5)

geo_ext_layers <- hee_extinction_layers(
  geo_process$dynamic_occupancy,
  geo_process$extinction_probability
)
write_table(geo_ext_layers$local_extinction, "geoprocess_local_extinction_proxy.csv")
write_table(geo_ext_layers$regional_extinction, "geoprocess_regional_extinction_proxy.csv")
write_table(geo_ext_layers$lineage_extinction, "geoprocess_lineage_extinction_proxy.csv")
ext_plot <- do.call(rbind, list(
  transform(stats::aggregate(local_extinction_risk ~ time_ma,
                             geo_ext_layers$local_extinction,
                             mean, na.rm = TRUE),
            scale = "local", value = local_extinction_risk)[, c("time_ma", "scale", "value")],
  transform(stats::aggregate(regional_extinction_proxy ~ time_ma,
                             geo_ext_layers$regional_extinction,
                             mean, na.rm = TRUE),
            scale = "regional", value = regional_extinction_proxy)[, c("time_ma", "scale", "value")],
  transform(stats::aggregate(lineage_extinction_proxy ~ time_ma,
                             geo_ext_layers$lineage_extinction,
                             mean, na.rm = TRUE),
            scale = "lineage", value = lineage_extinction_proxy)[, c("time_ma", "scale", "value")]
))
save_plot(ggplot2::ggplot(ext_plot,
                          ggplot2::aes(time_ma, value, color = scale)) +
            ggplot2::geom_line(linewidth = 0.75) +
            ggplot2::scale_x_reverse() +
            ggplot2::scale_y_continuous(limits = c(0, 1)) +
            ggplot2::labs(title = "Local, regional, and lineage extinction diagnostics",
                          x = "Time (Ma)", y = "Mean proxy",
                          color = "Scale"),
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
save_plot(ggplot2::ggplot(cmg_counts,
                          ggplot2::aes(time_ma, n_region_times,
                                       color = dominant_region_status)) +
            ggplot2::geom_line(linewidth = 0.75) +
            ggplot2::geom_point(size = 1.1) +
            ggplot2::scale_x_reverse() +
            ggplot2::labs(title = "Cradle, museum, grave, source and sink proxy states",
                          x = "Time (Ma)", y = "Region-time count",
                          color = "Status"),
          "geoprocess_cradle_museum_grave_through_time.png", 9, 5)

geo_dyn_family <- hee_dynamic_model_family()
geo_hypotheses <- hee_geoprocess_hypotheses()
geo_model_cmp <- hee_geoprocess_model_comparison(
  list(
    M1_env_only = M1,
    M2_env_phylo = M2,
    M3_env_phylo_bgb_no_dispersal = M3,
    M4_env_phylo_bgb_static_dispersal = M4,
    M5_env_phylo_bgb_dynamic_dispersal = M5
  ),
  richness = richness
)
write_table(geo_dyn_family, "geoprocess_dynamic_model_family_M0_M8.csv")
write_table(geo_hypotheses, "geoprocess_mechanism_hypotheses_H1_H8.csv")
write_table(geo_model_cmp, "geoprocess_M1_M5_model_comparison.csv")
save_plot(ggplot2::ggplot(geo_dyn_family,
                          ggplot2::aes(model_id, implementation_status,
                                       fill = implementation_status)) +
            ggplot2::geom_tile(color = "white", linewidth = 0.3) +
            ggplot2::geom_text(ggplot2::aes(label = model_name), size = 2.6) +
            ggplot2::labs(title = "Dynamic Earth-Biota model family M0-M8",
                          x = "Model", y = "Implementation status",
                          fill = "Status") +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30,
                                                               hjust = 1),
                           legend.position = "bottom"),
          "geoprocess_dynamic_model_family_M0_M8.png", 10, 5.5)
save_plot(ggplot2::ggplot(geo_hypotheses,
                          ggplot2::aes(hypothesis_id, hypothesis,
                                       fill = hypothesis_id)) +
            ggplot2::geom_tile(color = "white", linewidth = 0.3,
                               show.legend = FALSE) +
            ggplot2::geom_text(ggplot2::aes(label = hypothesis_id),
                               color = "white", fontface = "bold") +
            ggplot2::labs(title = "Geological-process mechanism hypotheses H1-H8",
                          x = "Hypothesis", y = "Mechanism") +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30,
                                                               hjust = 1)),
          "geoprocess_mechanism_hypotheses_H1_H8.png", 9.5, 5.5)
save_plot(ggplot2::ggplot(geo_model_cmp,
                          ggplot2::aes(model_id, mean_probability,
                                       fill = model_id)) +
            ggplot2::geom_col(show.legend = FALSE) +
            ggplot2::coord_flip() +
            ggplot2::labs(title = "M1-M5 geological-process model comparison",
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
process_long <- stats::reshape(process_summary,
                               varying = process_cols,
                               v.names = "mean_value",
                               timevar = "process",
                               times = process_cols,
                               direction = "long")
process_long <- process_long[is.finite(process_long$mean_value), , drop = FALSE]
save_plot(ggplot2::ggplot(process_long,
                          ggplot2::aes(time_ma, mean_value, color = process)) +
            ggplot2::geom_line(linewidth = 0.65) +
            ggplot2::scale_x_reverse() +
            ggplot2::labs(title = "Complete geological-process diagnostics through time",
                          x = "Time (Ma)", y = "Mean process value",
                          color = "Process"),
          "geoprocess_complete_process_summary.png", 9, 5.5)
cer <- Reduce(function(x, y) merge(x, y, by = "time_ma", all = TRUE,
                                   sort = FALSE), list(
  stats::aggregate(colonisation_probability ~ time_ma,
                   geo_process$colonisation_probability, mean, na.rm = TRUE),
  stats::aggregate(extinction_probability ~ time_ma,
                   geo_process$extinction_probability, mean, na.rm = TRUE),
  stats::aggregate(rescue_effect ~ time_ma,
                   geo_process$rescue_effect, mean, na.rm = TRUE)
))
cer_long <- stats::reshape(cer,
                           varying = c("colonisation_probability",
                                       "extinction_probability",
                                       "rescue_effect"),
                           v.names = "mean_value",
                           timevar = "process",
                           times = c("colonisation_probability",
                                     "extinction_probability",
                                     "rescue_effect"),
                           direction = "long")
cer_long <- cer_long[is.finite(cer_long$mean_value), , drop = FALSE]
save_plot(ggplot2::ggplot(cer_long,
                          ggplot2::aes(time_ma, mean_value, color = process)) +
            ggplot2::geom_line(linewidth = 0.8) +
            ggplot2::scale_x_reverse() +
            ggplot2::labs(title = "Colonisation, extinction and rescue diagnostics",
                          x = "Time (Ma)", y = "Mean probability/index",
                          color = NULL),
          "geoprocess_colonisation_extinction_rescue.png", 8, 5)
spec_cols <- intersect(c("vicariance_opportunity", "founder_event_opportunity",
                         "in_situ_speciation_index", "radiation_opportunity",
                         "speciation_opportunity"), names(geo_process$speciation_opportunity))
spec_time <- stats::aggregate(geo_process$speciation_opportunity[, spec_cols, drop = FALSE],
                              geo_process$speciation_opportunity[, "time_ma", drop = FALSE],
                              mean, na.rm = TRUE)
spec_long <- stats::reshape(spec_time,
                            varying = spec_cols,
                            v.names = "mean_value",
                            timevar = "speciation_component",
                            times = spec_cols,
                            direction = "long")
spec_long <- spec_long[is.finite(spec_long$mean_value), , drop = FALSE]
save_plot(ggplot2::ggplot(spec_long,
                          ggplot2::aes(time_ma, mean_value,
                                       color = speciation_component)) +
            ggplot2::geom_line(linewidth = 0.8) +
            ggplot2::scale_x_reverse() +
            ggplot2::labs(title = "Speciation opportunity components through time",
                          x = "Time (Ma)", y = "Opportunity index",
                          color = "Component"),
          "geoprocess_speciation_opportunity_through_time.png", 8.5, 5)
pool <- geo_process$species_pool
save_plot(ggplot2::ggplot(pool,
                          ggplot2::aes(time_ma, species_pool_size,
                                       color = region)) +
            ggplot2::geom_line(linewidth = 0.75) +
            ggplot2::scale_x_reverse() +
            ggplot2::labs(title = "Regional species-pool dynamics",
                          x = "Time (Ma)", y = "Expected species-pool size",
                          color = "Region"),
          "geoprocess_species_pool_through_time.png", 8, 5)
lim_region <- geo_process$limitation_summary
lim_region_plot <- lim_region[is.finite(lim_region$mechanism_fraction), , drop = FALSE]
save_plot(ggplot2::ggplot(lim_region_plot,
                          ggplot2::aes(time_ma, mechanism_fraction,
                                       color = dominant_community_mechanism)) +
            ggplot2::geom_line(linewidth = 0.75) +
            ggplot2::scale_x_reverse() +
            ggplot2::labs(title = "Region-level dominant limitation through time",
                          x = "Time (Ma)", y = "Weighted mechanism fraction",
                          color = "Mechanism"),
          "geoprocess_complete_limitation_region_through_time.png", 9, 5)

message("Phase 5b: method stress tests for geological-process functions")
method_geo_suite <- hee_geoprocess_scenario_suite(
  scenarios = c("normal", "extreme", "missing_optional", "duplicate_key",
                "single_time", "uneven_time", "shuffled_time",
                "land_appearance_loss"),
  species = species[seq_len(min(5, length(species)))],
  seed = 20260727
)
saveRDS(method_geo_suite,
        file.path(rds_dir, "method_geoprocess_scenario_suite.rds"),
        compress = "xz")
write_table(method_geo_suite$scenario_summary,
            "method_geoprocess_scenario_summary.csv")
write_table(method_geo_suite$validation_table,
            "method_geoprocess_scenario_validation.csv")
write_table(method_geo_suite$formula_catalog,
            "method_geoprocess_formula_catalog.csv")
write_table(method_geo_suite$interpretation_table,
            "method_geoprocess_interpretation_table.csv")
for (nm in names(method_geo_suite$input_tables)) {
  write_table(method_geo_suite$input_tables[[nm]],
              paste0("method_geoprocess_scenario_input_", nm, ".csv"))
}
for (nm in names(method_geo_suite$result_tables)) {
  write_table(method_geo_suite$result_tables[[nm]],
              paste0("method_geoprocess_scenario_result_", nm, ".csv"))
}
method_geo_plots <- plot_geoprocess_scenario_suite(method_geo_suite)
for (nm in names(method_geo_plots)) {
  save_plot(method_geo_plots[[nm]],
            paste0("method_geoprocess_scenario_", nm, ".png"),
            width = if (nm == "process_heatmap") 10 else 8,
            height = if (nm == "process_heatmap") 5.5 else 5)
}
method_geo_evidence <- method_geo_suite$interpretation_table
method_geo_evidence$input_table <- "method_geoprocess_scenario_input_landscape_state.csv; method_geoprocess_scenario_input_projection.csv"
method_geo_evidence$intermediate_table <- "method_geoprocess_scenario_result_region_state.csv"
method_geo_evidence$final_table <- ifelse(
  method_geo_evidence$module == "network_metrics",
  "geoprocess_network_summary.csv; geoprocess_corridor_persistence.csv",
  ifelse(method_geo_evidence$module == "extinction_layers",
         "geoprocess_local_extinction_proxy.csv; geoprocess_regional_extinction_proxy.csv; geoprocess_lineage_extinction_proxy.csv",
         ifelse(method_geo_evidence$module == "cradle_museum_grave",
                "geoprocess_cradle_museum_grave.csv",
                ifelse(method_geo_evidence$module == "dynamic_model_family",
                       "geoprocess_dynamic_model_family_M0_M8.csv",
                       ifelse(method_geo_evidence$module == "mechanism_hypotheses",
                              "geoprocess_mechanism_hypotheses_H1_H8.csv",
                              "method_geoprocess_scenario_result_process_summary.csv"))))
)
method_geo_evidence$figure <- paste0(
  "method_geoprocess_scenario_",
  c(
    arena = "landscape_events.png",
    land_age = "landscape_events.png",
    area_heterogeneity = "process_heatmap.png",
    connectivity = "network_fragmentation.png",
    network_metrics = "network_fragmentation.png",
    isolation = "network_fragmentation.png",
    colonisation = "process_heatmap.png",
    rescue_extinction = "extinction_layers.png",
    extinction_layers = "extinction_layers.png",
    species_pool = "process_heatmap.png",
    cradle_museum_grave = "cradle_museum_grave.png",
    speciation_opportunity = "process_heatmap.png",
    dynamic_model_family = "process_heatmap.png",
    mechanism_hypotheses = "process_heatmap.png",
    dynamic_occupancy = "dynamic_occupancy.png",
    richness_turnover_refugia = "process_heatmap.png",
    limitation = "process_heatmap.png",
    reliability = "process_heatmap.png",
    external_processes = "validation_status.png"
  )[match(method_geo_evidence$module,
          names(c(
            arena = "landscape_events.png",
            land_age = "landscape_events.png",
            area_heterogeneity = "process_heatmap.png",
            connectivity = "network_fragmentation.png",
            network_metrics = "network_fragmentation.png",
            isolation = "network_fragmentation.png",
            colonisation = "process_heatmap.png",
            rescue_extinction = "extinction_layers.png",
            extinction_layers = "extinction_layers.png",
            species_pool = "process_heatmap.png",
            cradle_museum_grave = "cradle_museum_grave.png",
            speciation_opportunity = "process_heatmap.png",
            dynamic_model_family = "process_heatmap.png",
            mechanism_hypotheses = "process_heatmap.png",
            dynamic_occupancy = "dynamic_occupancy.png",
            richness_turnover_refugia = "process_heatmap.png",
            limitation = "process_heatmap.png",
            reliability = "process_heatmap.png",
            external_processes = "validation_status.png"
          )))]
)
method_geo_evidence$automatic_validation <- "method_geoprocess_scenario_validation.csv"
method_geo_evidence$case_status <- "Method stress-test only; not a real-data geological estimate."
write_table(method_geo_evidence, "method_geoprocess_module_evidence_index.csv")

message("Phase 6: fossil validation and report")
if (!all(c("species", "time_ma") %in% names(fossils))) {
  stop("fossils.csv must contain species and time_ma columns.", call. = FALSE)
}
fv <- fossils
fv$note <- "Join fossil coordinates to projected cells before formal validation."
write_table(fv, "fossil_validation.csv")

inventory <- data.frame(
  path = list.files(output_root, recursive = TRUE, full.names = TRUE),
  stringsAsFactors = FALSE
)
inventory$type <- tools::file_ext(inventory$path)
write_table(inventory, "output_file_inventory.csv")

report_rmd <- file.path(rep_dir, "case02_realdata_full_540Ma_report.Rmd")
writeLines(c(
  "---",
  "title: 'HmscEcoEvo Case 2 Real-Data Full 540-0 Ma Workflow'",
  "output: word_document",
  "params:",
  "  output_root: '.'",
  "---",
  "",
  "```{r setup, include=FALSE}",
  "knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)",
  "fig_dir <- file.path(params$output_root, 'figures')",
  "tab_dir <- file.path(params$output_root, 'tables')",
  "```",
  "",
  "# Executive summary",
  "This real-data template used user-supplied community, environment, trait, phylogeny, range, plate, fossil, and BioGeoBEARS files. No simulated fallback was used.",
  "",
  "# HMSC inputs",
  "![Y matrix](`r file.path(fig_dir, '01_Y_heatmap.png')`)",
  "![X correlation](`r file.path(fig_dir, '04_X_correlation.png')`)",
  "![Trait heatmap](`r file.path(fig_dir, '06_traits_heatmap.png')`)",
  "![Phylogenetic correlation](`r file.path(fig_dir, '09_phylo_correlation_matrix.png')`)",
  "",
  "# HMSC posterior and validation",
  "![Beta](`r file.path(fig_dir, '13_hmsc_beta_heatmap.png')`)",
  "![Gamma](`r file.path(fig_dir, '14_hmsc_gamma_heatmap.png')`)",
  "![Omega](`r file.path(fig_dir, '16_hmsc_omega_heatmap.png')`)",
  "![Observed vs predicted](`r file.path(fig_dir, '21_hmsc_observed_vs_predicted.png')`)",
  "",
  "# BioGeoBEARS and projection",
  "![BGB weights](`r file.path(fig_dir, '23_bgb_model_weights.png')`)",
  "![M1-M5 richness](`r file.path(fig_dir, 'M1_M5_richness_through_time.png')`)",
  "",
  "# Geological-process diagnostics",
  "![Landscape arena](`r file.path(fig_dir, 'geoprocess_landscape_summary.png')`)",
  "![Ecological opportunity](`r file.path(fig_dir, 'geoprocess_ecological_opportunity.png')`)",
  "![Connectivity and isolation](`r file.path(fig_dir, 'geoprocess_connectivity_isolation.png')`)",
  "![Network fragmentation](`r file.path(fig_dir, 'geoprocess_network_fragmentation_through_time.png')`)",
  "![Complete process summary](`r file.path(fig_dir, 'geoprocess_complete_process_summary.png')`)",
  "![Colonisation extinction rescue](`r file.path(fig_dir, 'geoprocess_colonisation_extinction_rescue.png')`)",
  "![Extinction layers](`r file.path(fig_dir, 'geoprocess_extinction_layers_through_time.png')`)",
  "![Speciation opportunity](`r file.path(fig_dir, 'geoprocess_speciation_opportunity_through_time.png')`)",
  "![Species pool](`r file.path(fig_dir, 'geoprocess_species_pool_through_time.png')`)",
  "![Cradle museum grave](`r file.path(fig_dir, 'geoprocess_cradle_museum_grave_through_time.png')`)",
  "![Dynamic model family](`r file.path(fig_dir, 'geoprocess_dynamic_model_family_M0_M8.png')`)",
  "![Mechanism hypotheses](`r file.path(fig_dir, 'geoprocess_mechanism_hypotheses_H1_H8.png')`)",
  "![Limitation region](`r file.path(fig_dir, 'geoprocess_complete_limitation_region_through_time.png')`)",
  "",
  "# Geological-process method stress tests",
  "These stress tests are deterministic software and scientific-logic checks. They are not real geological estimates and are not used as fitted palaeo-results.",
  "```{r}",
  "knitr::kable(read.csv(file.path(tab_dir, 'method_geoprocess_scenario_summary.csv')))",
  "```",
  "![Scenario process heatmap](`r file.path(fig_dir, 'method_geoprocess_scenario_process_heatmap.png')`)",
  "![Scenario validation](`r file.path(fig_dir, 'method_geoprocess_scenario_validation_status.png')`)",
  "![Landscape events](`r file.path(fig_dir, 'method_geoprocess_scenario_landscape_events.png')`)",
  "![Dynamic occupancy](`r file.path(fig_dir, 'method_geoprocess_scenario_dynamic_occupancy.png')`)",
  "![Network fragmentation](`r file.path(fig_dir, 'method_geoprocess_scenario_network_fragmentation.png')`)",
  "![Extinction layers](`r file.path(fig_dir, 'method_geoprocess_scenario_extinction_layers.png')`)",
  "![Cradle museum grave](`r file.path(fig_dir, 'method_geoprocess_scenario_cradle_museum_grave.png')`)",
  "```{r}",
  "knitr::kable(head(read.csv(file.path(tab_dir, 'method_geoprocess_module_evidence_index.csv')), 20))",
  "```",
  "",
  "# Limitations",
  "BioGeoBEARS results must be supplied by the user or fitted through a custom runner. Plate reconstructions are based on user-supplied plate points unless rgplates/pygplates integration is added.",
  "",
  "# Session info",
  "```{r}",
  "sessionInfo()",
  "```"
), report_rmd)
report_docx <- hee_render_word_report(report_rmd,
                                       file.path(rep_dir, "hmscecoevo_realdata_full_540Ma_workflow.docx"),
                                       params = list(output_root = output_root))
report_docx <- hee_write_full_workflow_report(
  output_root,
  output_file = file.path(rep_dir, "hmscecoevo_realdata_full_540Ma_workflow.docx"),
  include_all_figures = TRUE,
  max_table_rows = 8
)

summary <- data.frame(output_root = output_root, report_docx = report_docx,
                      ordered_results = file.path(output_root, "ordered_results"),
                      n_figures_png = length(list.files(fig_dir, "\\.png$")),
                      n_tables = length(list.files(tab_dir, "\\.csv$")),
                      n_time_slices = length(times), n_species = length(species),
                      n_sites = nrow(Y), quick = quick,
                      projection_mode = projection_mode,
                      stringsAsFactors = FALSE)
write_table(summary, "case02_run_summary.csv")
message("Writing ordered output copies")
hee_write_ordered_outputs(output_root, include_rds = TRUE, overwrite = TRUE)
print(summary)
