#!/usr/bin/env Rscript

# Case 03 simulated demo for the rebuilt HmscEE nested framework.
# This script uses virtual data only when --analysis_mode=demo.

args_full <- commandArgs(FALSE)
cmd_file <- args_full[grep("^--file=", args_full)]
script_path <- if (length(cmd_file)) {
  normalizePath(sub("^--file=", "", cmd_file[[1]]), winslash = "/",
                mustWork = FALSE)
} else {
  normalizePath("scripts/case03_geoprocess_simulated_demo.R",
                winslash = "/", mustWork = FALSE)
}
pkg_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/",
                          mustWork = FALSE)
load_hmscecoevo <- function(script_path) {
  if (requireNamespace("HmscEcoEvo", quietly = TRUE)) {
    library(HmscEcoEvo)
    return(invisible(TRUE))
  }
  roots <- unique(normalizePath(c(
    file.path(dirname(script_path), ".."),
    file.path(dirname(script_path), "..", "..")
  ), winslash = "/", mustWork = FALSE))
  roots <- roots[file.exists(file.path(roots, "DESCRIPTION")) &
                   dir.exists(file.path(roots, "R"))]
  if (length(roots) && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(roots[[1]], quiet = TRUE)
    return(invisible(TRUE))
  }
  stop("HmscEcoEvo is not installed and a source package root could not be found.",
       call. = FALSE)
}
load_hmscecoevo(script_path)

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = "") {
  hit <- grep(paste0("^(--)?", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^(--)?", name, "="), "", hit[[1]])
}
as_bool <- function(x) tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")

analysis_mode <- tolower(arg_value("analysis_mode", "demo"))
if (!identical(analysis_mode, "demo")) {
  stop("This simulated demo only runs with --analysis_mode=demo. ",
       "Use case03_geoprocess_realdata_past_to_present.R for real data.",
       call. = FALSE)
}
quick <- as_bool(arg_value("quick", "TRUE"))
seed <- as.integer(arg_value("seed", "20260817"))
set.seed(seed)
output <- normalizePath(arg_value("output",
                                  if (quick) "outputs/case03_hmscee_nested_demo_quick" else
                                    "outputs/case03_hmscee_nested_demo_full"),
                        winslash = "/", mustWork = FALSE)
dirs <- file.path(output, c("tables", "metadata", "rds", "figures"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
tab_dir <- file.path(output, "tables")
fig_dir <- file.path(output, "figures")
meta_dir <- file.path(output, "metadata")
writeLines("demo", file.path(meta_dir, "analysis_mode.txt"))

times <- if (quick) c(100, 50, 0) else seq(540, 0, by = -20)
geography_scenarios <- c("g_landbridge", "g_no_bridge")
climate_scenarios <- c("c_cool", "c_warm")
lon <- seq(-40, 40, length.out = if (quick) 8 else 16)
lat <- seq(-20, 20, length.out = if (quick) 5 else 10)
earth <- expand.grid(geography_scenario = geography_scenarios,
                     climate_scenario = climate_scenarios,
                     lon = lon, lat = lat, time_ma = times,
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
earth$cell_id <- paste0("c", match(paste(earth$lon, earth$lat),
                                   unique(paste(earth$lon, earth$lat))))
earth$region <- ifelse(earth$lon < 0, "West", "East")
earth$land <- as.integer(!(earth$time_ma == max(times) & earth$lat > 10 &
                             earth$geography_scenario == "g_no_bridge"))
earth$bio1 <- as.numeric(scale(earth$lon / 40 + earth$time_ma / max(times) +
                                 ifelse(earth$climate_scenario == "c_warm",
                                        0.25, 0)))
earth$bio12 <- as.numeric(scale(earth$lat / 20 - earth$time_ma / max(times) -
                                  ifelse(earth$climate_scenario == "c_warm",
                                         0.15, 0)))
earth$climate_weight_given_geography <- ifelse(
  earth$climate_scenario == "c_cool", 0.55, 0.45
)
E <- hee_paleo_earth_state(
  earth, env_cols = c("bio1", "bio12"),
  geography_scenario_col = "geography_scenario",
  climate_scenario_col = "climate_scenario",
  climate_weight_col = "climate_weight_given_geography"
)

lineages <- c("lin_A", "lin_B")
ancestor_traits <- expand.grid(
  lineage = lineages, response_draw = c("s1", "s2"), time_ma = times,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
ancestor_traits$leaf_economics <- ifelse(ancestor_traits$lineage == "lin_A",
                                         1.0, -0.4) +
  ifelse(ancestor_traits$response_draw == "s2", 0.1, 0) -
  ancestor_traits$time_ma / max(times) * 0.2
ancestor_traits$wood_density <- ifelse(ancestor_traits$lineage == "lin_B",
                                       0.9, 0.2) +
  ancestor_traits$time_ma / max(times) * 0.15
Gamma <- matrix(
  c(0.9, 0.35,
    -0.45, 0.7),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("leaf_economics", "wood_density"),
                  c("bio1", "bio12"))
)
residual_response <- ancestor_traits[, c("lineage", "time_ma",
                                         "response_draw")]
residual_response$bio1 <- ifelse(residual_response$lineage == "lin_B",
                                 -0.15, 0.05)
residual_response$bio12 <- ifelse(residual_response$response_draw == "s2",
                                  0.08, 0)
responses <- hee_trait_mediated_ancestral_response(
  ancestor_traits, Gamma, residual_response = residual_response
)
S <- hee_lineage_suitability(E, responses, basis_cols = c("bio1", "bio12"))

region_history <- expand.grid(
  geography_scenario = geography_scenarios,
  lineage = lineages, region = c("West", "East"), time_ma = times,
  history_draw = c("h1", "h2"), KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
region_history$R_region <- with(region_history, as.integer(
    (lineage == "lin_A" & region == "West") |
    (lineage == "lin_A" & region == "East" & time_ma <=
       ifelse(geography_scenario == "g_landbridge", 50, 0)) |
    (lineage == "lin_B" & region == "East") |
    (lineage == "lin_B" & region == "West" & history_draw == "h2" &
       geography_scenario == "g_landbridge" & time_ma == 0)
))
region_history$geography_weight <- ifelse(
  region_history$geography_scenario == "g_landbridge", 0.65, 0.35
)
region_history$history_weight_given_geography <- ifelse(
  region_history$history_draw == "h1",
  ifelse(region_history$geography_scenario == "g_landbridge", 0.75, 0.45),
  ifelse(region_history$geography_scenario == "g_landbridge", 0.25, 0.55)
)
R <- hee_bsm_region_history(region_history)

cells <- unique(E[, c("geography_scenario", "cell_id", "lon", "lat", "region")])
paths0 <- merge(cells, cells, by = c("geography_scenario", "region"),
                suffixes = c("_from", "_to"))
paths0$least_cost_distance_km <- hee_great_circle_distance_km(
  paths0$lon_from, paths0$lat_from, paths0$lon_to, paths0$lat_to
)
paths0$least_cost_distance_km <- paths0$least_cost_distance_km *
  ifelse(paths0$geography_scenario == "g_no_bridge", 1.35, 1)
paths0 <- paths0[paths0$least_cost_distance_km <= if (quick) 1800 else 900, ]
paths <- do.call(rbind, lapply(tail(times, -1), function(tt) {
  z <- paths0
  z$time_ma <- tt
  z
}))
paths <- merge(expand.grid(lineage = lineages, row_id = seq_len(nrow(paths)),
                           KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE),
               cbind(row_id = seq_len(nrow(paths)), paths), by = "row_id")
names(paths)[names(paths) == "cell_id_from"] <- "from_cell_id"
names(paths)[names(paths) == "cell_id_to"] <- "to_cell_id"
names(paths)[names(paths) == "region"] <- "from_region"
paths$to_region <- paths$from_region
K <- hee_within_region_movement_kernel(
  paths[, c("geography_scenario", "lineage", "from_cell_id", "to_cell_id", "from_region",
            "to_region", "time_ma", "least_cost_distance_km")],
  delta_t = if (quick) 50 else 20,
  movement_intercept = -0.5,
  distance_decay = 0.001
)

fit <- hee_nested_region_cell_occupancy(S, R, E, movement_kernel = K,
                                        initial_rho = 0.35)
summary <- fit$summary
regional_pool <- hee_region_pool_from_cell_occupancy(summary)

extrapolation <- unique(E[, c("cell_id", "time_ma")])
extrapolation$extrapolation_score <- as.numeric(extrapolation$time_ma ==
                                                  max(extrapolation$time_ma))
summary_with_uncertainty <- hee_extrapolation_uncertainty_report(
  summary, extrapolation
)

utils::write.csv(hee_hmscee_process_catalog(),
                 file.path(tab_dir, "00.0_hmscee_five_process_catalog.csv"),
                 row.names = FALSE)
utils::write.csv(hee_hmscee_formula_catalog(),
                 file.path(tab_dir, "00.1_hmscee_nested_formula_catalog.csv"),
                 row.names = FALSE)
utils::write.csv(E, file.path(tab_dir, "01.1_paleo_earth_X_H_state.csv"),
                 row.names = FALSE)
utils::write.csv(ancestor_traits,
                 file.path(tab_dir, "01.2_ancestral_traits_T.csv"),
                 row.names = FALSE)
utils::write.csv(responses,
                 file.path(tab_dir, "01.3_trait_mediated_response_beta.csv"),
                 row.names = FALSE)
utils::write.csv(R, file.path(tab_dir, "02.1_bsm_region_history_R.csv"),
                 row.names = FALSE)
utils::write.csv(K, file.path(tab_dir, "03.1_within_region_movement_kernel_K.csv"),
                 row.names = FALSE)
utils::write.csv(fit$draws, file.path(tab_dir, "04.1_cell_occupancy_draws_q.csv"),
                 row.names = FALSE)
utils::write.csv(summary_with_uncertainty,
                 file.path(tab_dir, "04.2_cell_probability_summary_p.csv"),
                 row.names = FALSE)
utils::write.csv(regional_pool,
                 file.path(tab_dir, "05.1_region_pool_derived_from_cells.csv"),
                 row.names = FALSE)
saveRDS(fit, file.path(output, "rds", "case03_hmscee_nested_demo_fit.rds"))

metadata <- data.frame(
  item = c("analysis_mode", "seed", "time_axis", "n_lineages", "n_cells",
           "framework", "biological_process_engines",
           "external_layers_not_processes"),
  value = c("demo", seed, paste(times, collapse = ";"), length(lineages),
            length(unique(E$cell_id)),
            "Five-process HmscEE with BioGeoBEARS/BSM region history + within-region dynamic occupancy",
            "P1 abiotic filtering; P2 dispersal-colonisation; P3 biotic interactions; P4 niche/trait evolution; P5 speciation-lineage extinction",
            "Climate, geology, geography, habitat, dated trees, BioGeoBEARS, fossils, extrapolation and phyloregion diagnostics"),
  stringsAsFactors = FALSE
)
utils::write.csv(metadata, file.path(meta_dir, "case03_hmscee_demo_metadata.csv"),
                 row.names = FALSE)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  p <- ggplot2::ggplot(summary_with_uncertainty,
                       ggplot2::aes(lon, lat, fill = probability_mean)) +
    ggplot2::geom_tile() +
    ggplot2::facet_grid(lineage ~ time_ma) +
    ggplot2::coord_equal() +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), na.value = "grey90") +
    ggplot2::labs(
      title = "Case 03 nested HmscEE demo",
      subtitle = "p = E[R_region x q]; extrapolation is reported separately",
      x = "Longitude", y = "Latitude", fill = "p"
    ) +
    ggplot2::theme_minimal(base_size = 9)
  ggplot2::ggsave(file.path(fig_dir, "case03_nested_cell_probability_maps.png"),
                  p, width = 10, height = 6, dpi = 160)
}

message("Case 03 nested demo written to: ", output)
