#!/usr/bin/env Rscript

# Reweight complete Case05 history shards by a declared modern endpoint. This
# script does not change an individual forward trajectory; it writes an
# auditable weight table for case04_merge_posterior_draw_shards.R.

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x)) y else x

parse_args <- function(x) {
  out <- list()
  for (arg in x) {
    if (!grepl("^--", arg)) next
    key_value <- sub("^--", "", arg)
    pos <- regexpr("=", key_value, fixed = TRUE)
    if (pos < 0L) out[[key_value]] <- TRUE else {
      out[[substr(key_value, 1L, pos - 1L)]] <-
        substr(key_value, pos + 1L, nchar(key_value))
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
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(x, path)
  } else utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
  invisible(path)
}
coord_key <- function(lon, lat) {
  paste(formatC(as.numeric(lon), format = "f", digits = 8),
        formatC(as.numeric(lat), format = "f", digits = 8), sep = "\r")
}

cfg <- parse_args(args)
pkg_root <- normalizePath(cfg$pkg_root %||%
  "C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo", winslash = "/", mustWork = TRUE)
if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(pkg_root, quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(HmscEcoEvo))
}
shard_root <- normalizePath(cfg$shard_root %||% "", winslash = "/", mustWork = TRUE)
output <- safe_dir(cfg$output %||% dirname(shard_root))
data_role <- cfg$data_role %||% "training_diagnostic"
valid_roles <- c("independent_validation", "spatial_holdout", "training_diagnostic")
if (!data_role %in% valid_roles) {
  stop("data_role must be one of: ", paste(valid_roles, collapse = ", "), ".",
       call. = FALSE)
}
plant_input_dir <- cfg$plant_input_dir %||%
  "C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo_API8_20260818/data_external/plant_200_multifamily_extant_20260820/prepared_inputs/case03_data_raw/Plant200_multifamily_extant_1deg_traits"
comm_path <- normalizePath(cfg$comm %||% file.path(plant_input_dir, "comm.csv"),
                           winslash = "/", mustWork = TRUE)
sites_path <- normalizePath(cfg$sites %||% file.path(plant_input_dir, "sites.csv"),
                            winslash = "/", mustWork = TRUE)
shards <- list.dirs(shard_root, recursive = FALSE, full.names = TRUE)
shards <- shards[file.exists(file.path(shards, "20_modern_endpoint_validation"))]
if (!length(shards)) {
  stop("No completed shards with modern endpoint payloads were found.", call. = FALSE)
}
shards <- normalizePath(shards, winslash = "/", mustWork = TRUE)
endpoint_paths <- vapply(shards, function(path) {
  files <- list.files(file.path(path, "20_modern_endpoint_validation"),
                      pattern = "^modern_tip_endpoint_0Ma_.*\\.rds$",
                      full.names = TRUE)
  if (length(files) != 1L) {
    stop("Shard must contain exactly one modern endpoint payload: ", path,
         call. = FALSE)
  }
  normalizePath(files[[1L]], winslash = "/", mustWork = TRUE)
}, character(1))
payloads <- lapply(endpoint_paths, readRDS)
first <- payloads[[1L]]
needed <- c("cell_id", "lon", "lat", "species",
            "dynamic_occupancy_probability")
if (length(setdiff(needed, names(first)))) {
  stop("The first endpoint payload has an invalid schema.", call. = FALSE)
}
for (payload in payloads[-1L]) {
  if (!identical(as.character(payload$cell_id), as.character(first$cell_id)) ||
      !identical(as.character(payload$species), as.character(first$species)) ||
      !identical(dim(payload$dynamic_occupancy_probability),
                 dim(first$dynamic_occupancy_probability))) {
    stop("Endpoint payload grid, species, or matrix dimensions differ among shards.",
         call. = FALSE)
  }
}
comm <- utils::read.csv(comm_path, check.names = FALSE, stringsAsFactors = FALSE)
sites <- utils::read.csv(sites_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("site_id", "lon", "lat") %in% names(sites)) ||
    !"site_id" %in% names(comm)) {
  stop("comm and sites must both have site_id; sites also needs lon and lat.",
       call. = FALSE)
}
if (anyDuplicated(comm$site_id) || anyDuplicated(sites$site_id) ||
    !setequal(as.character(comm$site_id), as.character(sites$site_id))) {
  stop("comm.csv and sites.csv must have the same unique site_id values.",
       call. = FALSE)
}
sites <- sites[match(comm$site_id, sites$site_id), c("site_id", "lon", "lat"),
               drop = FALSE]
species <- as.character(first$species)
missing_species <- setdiff(species, names(comm))
if (length(missing_species)) {
  stop("Endpoint species are absent from comm.csv: ",
       paste(missing_species, collapse = ", "), call. = FALSE)
}
endpoint_key <- coord_key(first$lon, first$lat)
site_key <- coord_key(sites$lon, sites$lat)
site_index <- match(site_key, endpoint_key)
if (anyNA(site_index)) {
  stop("At least one modern site has no exact 0 Ma endpoint grid match.",
       call. = FALSE)
}
observed <- as.matrix(comm[, species, drop = FALSE])
storage.mode(observed) <- "double"
rownames(observed) <- comm$site_id
particles <- setNames(lapply(payloads, function(payload) {
  x <- payload$dynamic_occupancy_probability[site_index, species, drop = FALSE]
  rownames(x) <- comm$site_id
  colnames(x) <- species
  x
}), shards)

power_supplied <- !is.null(cfg$likelihood_power)
likelihood_power <- suppressWarnings(as.numeric(cfg$likelihood_power %||%
  if (data_role == "training_diagnostic") 1 / nrow(observed) else 1))
if (!is.finite(likelihood_power) || likelihood_power < 0) {
  stop("likelihood_power must be finite and non-negative.", call. = FALSE)
}
weights <- hee_endpoint_particle_weights(
  particles = particles, observed = observed, data_role = data_role,
  likelihood_power = likelihood_power
)
weights$shard <- shards[match(weights$particle_id, names(particles))]
weights$endpoint_payload <- endpoint_paths[match(weights$particle_id, names(particles))]
weights$likelihood_power_source <- if (power_supplied) {
  "user_supplied"
} else if (data_role == "training_diagnostic") {
  "default_per_site_tempering_for_nonindependent_training_diagnostic"
} else "default_full_independent_endpoint_likelihood"
weights <- weights[, c("shard", "particle_id", "weight", "log_endpoint_score",
                       "likelihood_power", "likelihood_power_source",
                       "ensemble_effective_sample_size", "data_role",
                       "scientific_role", "endpoint_payload")]
audit_dir <- safe_dir(file.path(output, "20_modern_endpoint_validation"))
write_csv(weights, file.path(audit_dir, "complete_history_endpoint_weights.csv"))
write_csv(data.frame(
  site_id = comm$site_id, lon = sites$lon, lat = sites$lat,
  endpoint_cell_id = first$cell_id[site_index],
  endpoint_lon = first$lon[site_index], endpoint_lat = first$lat[site_index],
  stringsAsFactors = FALSE
), file.path(audit_dir, "complete_history_endpoint_site_match.csv"))
write_csv(data.frame(
  shard_root = shard_root, n_particles = length(particles),
  n_sites = nrow(observed), n_species = ncol(observed), data_role = data_role,
  likelihood_power = likelihood_power,
  scientific_boundary = if (data_role == "training_diagnostic") {
    "Weights are a non-independent endpoint diagnostic. Do not call them an empirical posterior."
  } else {
    "Weights condition complete trajectories on the declared endpoint observations."
  },
  stringsAsFactors = FALSE
), file.path(audit_dir, "complete_history_endpoint_weight_inputs.csv"))
cat("Case05 complete-history endpoint weights written: ", audit_dir, "\n", sep = "")
