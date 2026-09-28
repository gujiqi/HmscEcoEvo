#!/usr/bin/env Rscript

# Merge independently checkpointed Case04 response-draw shards.  Each shard
# has run the full time-ordered, complete 1 degree land grid for one genuine
# HMSC Beta posterior draw.  This script averages only aligned cell-time
# summaries; it never resamples cells or treats the 40 QA points as inputs.

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

read_csv <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    as.data.frame(data.table::fread(path, data.table = FALSE))
  } else utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_csv <- function(x, path) {
  safe_dir(dirname(path))
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(x, path)
  } else utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
  invisible(path)
}

cfg <- parse_args(args)
shard_root <- cfg$shard_root %||% ""
shard_paths_arg <- cfg$shard_paths %||% ""
output <- cfg$output %||% ""
if ((!dir.exists(shard_root) && !nzchar(shard_paths_arg)) || !nzchar(output)) {
  stop("Use --shard_root=<directory> or --shard_paths=<path1|path2>, plus --output=<scenario directory>.",
       call. = FALSE)
}
output <- safe_dir(output)

shards <- if (nzchar(shard_paths_arg)) {
  trimws(strsplit(shard_paths_arg, "|", fixed = TRUE)[[1L]])
} else {
  list.dirs(shard_root, recursive = FALSE, full.names = TRUE)
}
shards <- shards[nzchar(shards) & file.exists(file.path(shards, "16_state_space_inference"))]
if (!length(shards)) stop("No completed Case04 draw shards found.", call. = FALSE)
shards <- normalizePath(shards, winslash = "/", mustWork = TRUE)

# A supplied table can reweight whole, already completed history particles
# (for example alternative root patches) by a declared endpoint likelihood.
# Without one, the historical equal-weight aggregation remains unchanged.
weight_file <- cfg$shard_weights %||% ""
shard_weights <- rep(1 / length(shards), length(shards))
weight_role <- "equal_weight_complete_history_particles"
if (nzchar(weight_file)) {
  weight_file <- normalizePath(weight_file, winslash = "/", mustWork = TRUE)
  supplied_weights <- read_csv(weight_file)
  required_weight_cols <- c("shard", "weight")
  missing_weight_cols <- setdiff(required_weight_cols, names(supplied_weights))
  if (length(missing_weight_cols)) {
    stop("shard_weights is missing required columns: ",
         paste(missing_weight_cols, collapse = ", "), call. = FALSE)
  }
  supplied_weights$shard <- normalizePath(as.character(supplied_weights$shard),
                                           winslash = "/", mustWork = TRUE)
  if (anyDuplicated(supplied_weights$shard)) {
    stop("shard_weights contains duplicate shard paths.", call. = FALSE)
  }
  hit_weight <- match(shards, supplied_weights$shard)
  if (anyNA(hit_weight)) {
    stop("shard_weights does not contain every completed shard.", call. = FALSE)
  }
  shard_weights <- suppressWarnings(as.numeric(supplied_weights$weight[hit_weight]))
  if (any(!is.finite(shard_weights)) || any(shard_weights < 0) ||
      sum(shard_weights) <= 0) {
    stop("shard_weights must be finite, non-negative, and have positive total weight.",
         call. = FALSE)
  }
  shard_weights <- shard_weights / sum(shard_weights)
  weight_role <- if ("scientific_role" %in% names(supplied_weights)) {
    paste(unique(as.character(supplied_weights$scientific_role[hit_weight])),
          collapse = ";")
  } else "externally_supplied_shard_weights"
}

weighted_rows <- function(values, weights) {
  valid <- is.finite(values)
  numerator <- as.vector(replace(values, !valid, 0) %*% weights)
  denominator <- as.vector(valid %*% weights)
  out <- numerator / denominator
  out[denominator <= 0] <- NA_real_
  out
}

metric_rel <- function(path) {
  metric_dir <- file.path(path, "16_state_space_inference")
  csv <- list.files(metric_dir, pattern = "^cell_metrics_.*\\.csv$",
                    full.names = FALSE)
  rds <- list.files(metric_dir, pattern = "^cell_metrics_.*\\.rds$",
                    full.names = FALSE)
  if (length(csv) && length(rds)) {
    stop("A shard mixes CSV and RDS metric storage: ", path, call. = FALSE)
  }
  c(csv, rds)
}
read_metric <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) readRDS(path) else read_csv(path)
}
metric_names <- metric_rel(shards[[1L]])
if (!length(metric_names) || any(vapply(shards, function(x) {
  !setequal(metric_rel(x), metric_names)
}, logical(1)))) {
  stop("Every shard must contain the same per-time full-grid metric files.",
       call. = FALSE)
}
metric_output_names <- sub("\\.rds$", ".csv", metric_names,
                           ignore.case = TRUE)
if (anyDuplicated(metric_output_names)) {
  stop("Metric storage names cannot be converted unambiguously to final CSV names.",
       call. = FALSE)
}

fixed_columns <- c("cell_id", "time_ma", "lon", "lat", "H_state", "land_area_km2")
merge_one_time <- function(file_name) {
  tabs <- lapply(shards, function(x) {
    z <- read_metric(file.path(x, "16_state_space_inference", file_name))
    required <- c(fixed_columns, "mean_occupancy_probability")
    missing <- setdiff(required, names(z))
    if (length(missing)) {
      stop("Shard table missing columns: ", paste(missing, collapse = ", "),
           call. = FALSE)
    }
    z <- z[order(z$cell_id, z$time_ma), , drop = FALSE]
    rownames(z) <- NULL
    z
  })
  keys <- vapply(tabs, function(z) paste(z$cell_id, z$time_ma, sep = "\r"),
                 character(nrow(tabs[[1L]])))
  if (is.null(dim(keys))) keys <- matrix(keys, ncol = 1L)
  if (any(vapply(seq_len(ncol(keys)), function(i) !identical(keys[, i], keys[, 1L]),
                 logical(1)))) {
    stop("Cell-time keys differ among shards for ", file_name, ".", call. = FALSE)
  }
  out <- tabs[[1L]][, intersect(fixed_columns, names(tabs[[1L]])), drop = FALSE]
  numeric_cols <- names(tabs[[1L]])[vapply(tabs[[1L]], is.numeric, logical(1))]
  metrics <- setdiff(numeric_cols, fixed_columns)
  for (nm in metrics) {
    values <- do.call(cbind, lapply(tabs, `[[`, nm))
    out[[nm]] <- weighted_rows(values, shard_weights)
  }
  out$analysis_mode <- "global_no_bgb_complete_history_particle_aggregate"
  out$interpretation <-
    "weighted_mean_across_complete_1deg_history_particles;sampled_surviving_lineage_not_total_flora"
  out
}

metric_dir <- safe_dir(file.path(output, "16_state_space_inference"))
all_tables <- lapply(metric_names, merge_one_time)
names(all_tables) <- metric_names
for (i in seq_along(metric_names)) {
  write_csv(all_tables[[metric_names[[i]]]],
            file.path(metric_dir, metric_output_names[[i]]))
}
all_metrics <- do.call(rbind, all_tables)
write_csv(all_metrics, file.path(metric_dir, "cell_time_all_metrics.csv"))

audit_dir <- safe_dir(file.path(output, "01_input_audit"))
first_points <- file.path(shards[[1L]], "01_input_audit", "case04_diagnostic_points_by_time.csv")
if (file.exists(first_points)) {
  points <- read_csv(first_points)
  write_csv(points, file.path(audit_dir, "case04_diagnostic_points_by_time.csv"))
  # Cell identity and time are the stable keys. Numeric cell-area fields are
  # derived independently in audit and metric readers, so using them as exact
  # floating-point join keys can drop valid high-latitude diagnostic cells.
  point_metrics <- merge(points, all_metrics,
                         by = c("cell_id", "time_ma"),
                         all.x = TRUE, sort = FALSE)
  if (nrow(point_metrics) != nrow(points) ||
      any(!is.finite(point_metrics$mean_occupancy_probability))) {
    stop("Merged posterior diagnostic points did not join to complete-grid results.",
         call. = FALSE)
  }
  point_metrics$trace_role <-
    "quality_assurance_points_only_not_a_model_or_map_subset"
  write_csv(point_metrics,
            file.path(metric_dir, "case04_diagnostic_point_process_trace.csv"))
  coverage <- aggregate(cell_id ~ time_ma, point_metrics,
                        function(x) length(unique(x)))
  names(coverage)[2L] <- "n_diagnostic_points_joined"
  coverage$expected_diagnostic_points <-
    pmin(40L, as.integer(table(points$time_ma)[as.character(coverage$time_ma)]))
  coverage$trace_join_complete <- coverage$n_diagnostic_points_joined ==
    coverage$expected_diagnostic_points
  write_csv(coverage,
            file.path(metric_dir, "case04_diagnostic_point_trace_coverage.csv"))
}

metric_cols <- setdiff(names(all_metrics)[vapply(all_metrics, is.numeric, logical(1))],
                       fixed_columns)
time_summary <- stats::aggregate(all_metrics[, metric_cols, drop = FALSE],
                                 all_metrics[, "time_ma", drop = FALSE],
                                 function(z) mean(z, na.rm = TRUE))
names(time_summary)[-1L] <- paste0("global_mean_", names(time_summary)[-1L])
report_dir <- safe_dir(file.path(output, "22_report"))
write_csv(time_summary, file.path(report_dir, "case04_final_time_summary.csv"))

first_config_path <- file.path(shards[[1L]], "00_config", "case04_final_config.csv")
if (!file.exists(first_config_path)) stop("First shard is missing Case04 config.", call. = FALSE)
config <- read_csv(first_config_path)
set_config <- function(key, value) {
  where <- match(key, config$parameter)
  if (is.na(where)) {
    config <<- rbind(config, data.frame(parameter = key, value = as.character(value)))
  } else config$value[[where]] <<- as.character(value)
}
shard_configs <- lapply(shards, function(x) {
  read_csv(file.path(x, "00_config", "case04_final_config.csv"))
})
draw_ids <- vapply(shard_configs, function(x) {
  v <- x$value[match("response_draw_ids_requested", x$parameter)]
  if (length(v)) as.character(v) else NA_character_
}, character(1))
set_config("n_draws", length(shards))
set_config("response_draw_ids_requested", paste(draw_ids, collapse = ","))
set_config("posterior_aggregation", weight_role)
set_config("palaeo_input_mode", "streaming_complete_1deg_time_slices")
config_dir <- safe_dir(file.path(output, "00_config"))
write_csv(config, file.path(config_dir, "case04_final_config.csv"))

manifest <- data.frame(
  shard = shards,
  response_draw = draw_ids,
  aggregation_weight = shard_weights,
  aggregation_weight_role = weight_role,
  n_cell_time_rows = vapply(shards, function(x) {
    sum(vapply(metric_names, function(nm) {
      nrow(read_metric(file.path(x, "16_state_space_inference", nm)))
    }, numeric(1)))
  }, numeric(1)),
  stringsAsFactors = FALSE
)
write_csv(manifest, file.path(config_dir, "posterior_draw_shard_manifest.csv"))

# Preserve the draw-level P2 kernel audit used by Case05 to distinguish the
# finite-propagule particle scheme from the two deterministic movement schemes.
movement_paths <- file.path(shards, "11_dispersal", "movement_matrix_summary.csv")
if (all(file.exists(movement_paths))) {
  movement_summary <- do.call(rbind, lapply(seq_along(movement_paths), function(i) {
    x <- read_csv(movement_paths[[i]])
    x$response_draw <- draw_ids[[i]]
    x
  }))
  dispersal_dir <- safe_dir(file.path(output, "11_dispersal"))
  write_csv(movement_summary, file.path(dispersal_dir, "movement_matrix_summary.csv"))
}

# Preserve the full 0 Ma tip-by-cell state required for a real endpoint audit.
# The historical run stores one compressed matrix per posterior response draw;
# merge them sequentially so the entire draw ensemble is never held in memory.
endpoint_paths <- unlist(lapply(shards, function(path) {
  list.files(file.path(path, "20_modern_endpoint_validation"),
             pattern = "^modern_tip_endpoint_0Ma_.*\\.rds$", full.names = TRUE)
}), use.names = FALSE)
if (length(endpoint_paths)) {
  if (length(endpoint_paths) != length(shards)) {
    stop("Each posterior shard must contain exactly one 0 Ma tip endpoint payload.",
         call. = FALSE)
  }
  endpoint_first <- readRDS(endpoint_paths[[1L]])
  required_endpoint <- c(
    "response_draw", "time_ma", "cell_id", "lon", "lat", "species",
    "dynamic_occupancy_probability", "hmsc_fixed_effect_nowcast_probability"
  )
  missing_endpoint <- setdiff(required_endpoint, names(endpoint_first))
  if (length(missing_endpoint)) {
    stop("Modern endpoint payload is missing: ",
         paste(missing_endpoint, collapse = ", "), call. = FALSE)
  }
  dynamic_sum <- endpoint_first$dynamic_occupancy_probability * shard_weights[[1L]]
  nowcast_sum <- endpoint_first$hmsc_fixed_effect_nowcast_probability * shard_weights[[1L]]
  has_forward_raw <- "forward_raw_occupancy_probability" %in% names(endpoint_first)
  forward_raw_sum <- if (has_forward_raw) {
    endpoint_first$forward_raw_occupancy_probability * shard_weights[[1L]]
  } else NULL
  anchor_mode <- endpoint_first$terminal_anchor_mode %||% "none"
  anchor_weight <- endpoint_first$terminal_anchor_weight %||% 0
  anchor_role <- endpoint_first$terminal_anchor_role %||% "not_recorded"
  if (!all(dim(dynamic_sum) == c(length(endpoint_first$cell_id),
                                 length(endpoint_first$species))) ||
      !identical(dim(dynamic_sum), dim(nowcast_sum))) {
    stop("First modern endpoint payload has incompatible matrix dimensions.",
         call. = FALSE)
  }
  if (length(endpoint_paths) > 1L) for (i in 2:length(endpoint_paths)) {
    endpoint <- readRDS(endpoint_paths[[i]])
    same_grid <- identical(as.character(endpoint$cell_id),
                           as.character(endpoint_first$cell_id))
    same_species <- identical(as.character(endpoint$species),
                              as.character(endpoint_first$species))
    same_shape <- identical(dim(endpoint$dynamic_occupancy_probability),
                            dim(dynamic_sum)) &&
      identical(dim(endpoint$hmsc_fixed_effect_nowcast_probability), dim(nowcast_sum))
    if (!same_grid || !same_species || !same_shape) {
      stop("Modern endpoint payloads have non-identical cell or species axes.",
           call. = FALSE)
    }
    dynamic_sum <- dynamic_sum +
      endpoint$dynamic_occupancy_probability * shard_weights[[i]]
    nowcast_sum <- nowcast_sum +
      endpoint$hmsc_fixed_effect_nowcast_probability * shard_weights[[i]]
    endpoint_has_raw <- "forward_raw_occupancy_probability" %in% names(endpoint)
    if (!identical(endpoint_has_raw, has_forward_raw)) {
      stop("Modern endpoint payloads disagree on forward_raw_occupancy_probability.",
           call. = FALSE)
    }
    if (has_forward_raw) {
      if (!identical(dim(endpoint$forward_raw_occupancy_probability), dim(forward_raw_sum))) {
        stop("Modern endpoint forward_raw_occupancy_probability has incompatible dimensions.",
             call. = FALSE)
      }
      forward_raw_sum <- forward_raw_sum +
        endpoint$forward_raw_occupancy_probability * shard_weights[[i]]
    }
    if (!identical(endpoint$terminal_anchor_mode %||% "none", anchor_mode) ||
        !isTRUE(all.equal(endpoint$terminal_anchor_weight %||% 0, anchor_weight))) {
      stop("Modern endpoint payloads disagree on terminal-anchor settings.",
           call. = FALSE)
    }
  }
  endpoint_dir <- safe_dir(file.path(output, "20_modern_endpoint_validation"))
  endpoint_mean <- list(
    time_ma = 0,
    cell_id = endpoint_first$cell_id,
    lon = endpoint_first$lon,
    lat = endpoint_first$lat,
    species = endpoint_first$species,
    dynamic_occupancy_probability = dynamic_sum,
    hmsc_fixed_effect_nowcast_probability = nowcast_sum,
    forward_raw_occupancy_probability = forward_raw_sum,
    n_response_draws = length(endpoint_paths),
    dynamic_definition = paste(
      "weighted mean of global-grid dynamic-occupancy complete-history particles;",
      weight_role
    ),
    nowcast_definition = endpoint_first$nowcast_definition,
    terminal_anchor_mode = anchor_mode,
    terminal_anchor_weight = anchor_weight,
    terminal_anchor_role = anchor_role
  )
  saveRDS(endpoint_mean,
          file.path(endpoint_dir, "modern_tip_endpoint_0Ma_posterior_mean.rds"),
          compress = "gzip")
  write_csv(data.frame(
    response_draw = vapply(endpoint_paths, function(path) readRDS(path)$response_draw,
                           character(1)),
    aggregation_weight = shard_weights,
    aggregation_weight_role = weight_role,
    endpoint_path = normalizePath(endpoint_paths, winslash = "/", mustWork = TRUE),
    stringsAsFactors = FALSE
  ), file.path(endpoint_dir, "modern_tip_endpoint_draw_manifest.csv"))
}

validation <- data.frame(
  check = c("complete_draw_shards", "identical_cell_time_keys",
            "complete_1deg_grid_preserved", "diagnostic_points_not_model_subset"),
  status = c("PASS", "PASS", "PASS", "PASS"),
  detail = c(
    paste0(length(shards), " independently completed HMSC response draw shard(s) merged."),
    "Every merged time table has identical cell_id x time_ma keys across shards.",
    "No metric-table row was spatially sampled during aggregation.",
    "The 40 points per time are rebuilt only as a QA trace from full-grid summaries."
  ), stringsAsFactors = FALSE
)
engine_validation_paths <- file.path(shards, "00_config", "case04_final_validation_checks.csv")
if (any(!file.exists(engine_validation_paths))) {
  stop("Every completed posterior shard must include engine validation checks.", call. = FALSE)
}
engine_validation <- read_csv(engine_validation_paths[[1L]])
if (!all(c("check", "status", "detail") %in% names(engine_validation))) {
  stop("Shard engine validation table has an unexpected schema.", call. = FALSE)
}
engine_signature <- lapply(engine_validation_paths, function(path) {
  x <- read_csv(path)
  x <- x[order(x$check), c("check", "status"), drop = FALSE]
  paste(x$check, x$status, collapse = "|")
})
if (!all(vapply(engine_signature, identical, logical(1), engine_signature[[1L]]))) {
  stop("Posterior shards disagree on their engine validation status.", call. = FALSE)
}
validation$validation_source <- "posterior_shard_aggregation"
engine_validation$validation_source <- "per_draw_engine"
validation <- rbind(engine_validation[, names(validation), drop = FALSE], validation)
write_csv(validation, file.path(config_dir, "case04_final_validation_checks.csv"))

writeLines(c(
  "# Case04 posterior-draw shard aggregate",
  "",
  paste0("Merged full-grid HMSC posterior draw shards: ", length(shards)),
  "Every source shard used all valid 1 degree palaeo-land cells at each time slice.",
  paste0("Complete-history particle aggregation: ", weight_role, "."),
  "The 40 per-time points are diagnostic traces only and never enter dynamic calculation or map construction."
), file.path(report_dir, "case04_posterior_draw_shard_aggregate.md"))

cat("Merged ", length(shards), " complete-grid posterior draw shards into ", output, "\n", sep = "")
