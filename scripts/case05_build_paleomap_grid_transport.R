#!/usr/bin/env Rscript

# Build explicit 1 degree palaeogeographic cell-carriage tables for Case05.
#
# The PALEOMAP input is a published gridded plate reconstruction with one H3
# track per row and reconstructed longitude/latitude fields from 0 to 540 Ma.
# It is used here only to carry an existing grid occupancy to the next Earth
# slice. The map is target-centred (semi-Lagrangian): every younger land cell
# is followed backward along its nearest plate track to the older slice. This
# prevents a coarser H3 carrier grid from collapsing two 1 degree source cells
# into one target and leaving an artificial chequerboard of uncarried land.
# It does not create biological dispersal edges and it does not replace a
# polygon-overlap transport product when one becomes available.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || identical(x, "") || is.na(x)) y else x
}

parse_args <- function(args) {
  out <- list()
  for (a in args) {
    if (!grepl("^--", a)) next
    z <- sub("^--", "", a)
    eq <- regexpr("=", z, fixed = TRUE)
    if (eq < 0L) out[[z]] <- TRUE else {
      out[[substr(z, 1L, eq - 1L)]] <- substr(z, eq + 1L, nchar(z))
    }
  }
  out
}

as_num <- function(x, default = NA_real_) {
  if (is.null(x) || identical(x, "")) return(default)
  suppressWarnings(as.numeric(x)[1L])
}

as_int <- function(x, default = NA_integer_) {
  z <- as_num(x, default)
  if (!is.finite(z)) default else as.integer(z)
}

safe_dir <- function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE)
  normalizePath(x, winslash = "/", mustWork = FALSE)
}

write_csv <- function(x, path) {
  safe_dir(dirname(path))
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(x, path)
  } else {
    utils::write.csv(x, path, row.names = FALSE)
  }
  invisible(path)
}

ma_slug <- function(x) {
  gsub("\\.", "p", formatC(as.numeric(x), format = "f", digits = 4,
                               drop0trailing = TRUE))
}

unit_xyz <- function(lon, lat) {
  rad <- pi / 180
  lon <- as.numeric(lon) * rad
  lat <- as.numeric(lat) * rad
  cbind(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
}

chord_to_km <- function(x) {
  2 * 6371.0088 * asin(pmin(1, pmax(0, as.numeric(x) / 2)))
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
pkg_root <- normalizePath(
  args$pkg_root %||% "C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo",
  winslash = "/", mustWork = FALSE
)
slice_index_path <- args$palaeo_slice_index %||% file.path(
  pkg_root, "derived_inputs", "case04_plant200_palaeo_1deg_slices_20260919",
  "case04_palaeo_earth_slice_index.csv"
)
plate_tracks_path <- args$plate_tracks %||% file.path(
  "C:/Users/Google/Documents/HMSC-HIST/data/raw/plate_reconstruction",
  "palaeoverse_gridded_recon_zenodo", "PALEOMAP.RDS"
)
output <- safe_dir(args$output %||% file.path(
  pkg_root, "derived_inputs", "case05_paleomap_h3_grid_transport_1deg_v1"
))
start_ma <- as_num(args$start_ma, 325)
max_match_distance_km <- as_num(args$max_match_distance_km, 350)
time_selection <- tolower(args$time_selection %||% "all")
n_time_slices <- as_int(args$n_time_slices, NA_integer_)
if (!time_selection %in% c("all", "evenly_spaced")) {
  stop("time_selection must be 'all' or 'evenly_spaced'.", call. = FALSE)
}
if (identical(time_selection, "evenly_spaced") &&
    (!is.finite(n_time_slices) || n_time_slices < 2L)) {
  stop("n_time_slices must be at least two for evenly_spaced selection.",
       call. = FALSE)
}
if (!is.finite(start_ma) || start_ma <= 0 ||
    !is.finite(max_match_distance_km) || max_match_distance_km <= 0) {
  stop("start_ma and max_match_distance_km must be positive finite values.",
       call. = FALSE)
}
if (!requireNamespace("RANN", quietly = TRUE)) {
  stop("RANN is required to construct plate-track nearest-neighbour maps.",
       call. = FALSE)
}
if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("pkgload is required to load the current HmscEcoEvo source.", call. = FALSE)
}
pkgload::load_all(pkg_root, quiet = TRUE)

if (!file.exists(slice_index_path) || !file.exists(plate_tracks_path)) {
  stop("Missing palaeo slice index or PALEOMAP track RDS.", call. = FALSE)
}
index <- utils::read.csv(slice_index_path, stringsAsFactors = FALSE,
                         check.names = FALSE)
if (!all(c("time_ma", "file") %in% names(index))) {
  stop("palaeo_slice_index must contain time_ma and file.", call. = FALSE)
}
index$time_ma <- suppressWarnings(as.numeric(index$time_ma))
if (any(!is.finite(index$time_ma)) || anyDuplicated(index$time_ma) ||
    any(!file.exists(index$file))) {
  stop("palaeo_slice_index has invalid times or missing files.", call. = FALSE)
}
times <- sort(index$time_ma[index$time_ma <= start_ma + 1e-8], decreasing = TRUE)
if (identical(time_selection, "evenly_spaced") && length(times) > n_time_slices) {
  times <- times[unique(round(seq.int(1L, length(times), length.out = n_time_slices)))]
}
if (length(times) < 2L) stop("Fewer than two eligible palaeo time slices.", call. = FALSE)

tracks <- readRDS(plate_tracks_path)
if (!is.data.frame(tracks) || !all(c("h3", "lng", "lat") %in% names(tracks))) {
  stop("PALEOMAP tracks must be a data frame with h3, lng, and lat columns.",
       call. = FALSE)
}
track_at_time <- function(tm) {
  suffix <- if (abs(tm) < 1e-8) "" else paste0("_", as.integer(round(tm)))
  lon_col <- paste0("lng", suffix)
  lat_col <- paste0("lat", suffix)
  if (!all(c(lon_col, lat_col) %in% names(tracks))) {
    stop("PALEOMAP track table has no coordinates for ", tm, " Ma.", call. = FALSE)
  }
  z <- data.frame(
    track_id = as.character(tracks$h3),
    lon = suppressWarnings(as.numeric(tracks[[lon_col]])),
    lat = suppressWarnings(as.numeric(tracks[[lat_col]])),
    stringsAsFactors = FALSE
  )
  z[is.finite(z$lon) & is.finite(z$lat), , drop = FALSE]
}
read_slice <- function(tm) {
  p <- index$file[[match(tm, index$time_ma)]]
  z <- as.data.frame(readRDS(p))
  need <- c("cell_id", "time_ma", "lon", "lat", "H_state", "land_area_km2")
  if (!all(need %in% names(z))) stop("Earth slice lacks required columns at ", tm, " Ma.")
  z <- z[is.finite(z$lon) & is.finite(z$lat) &
           is.finite(z$H_state) & z$H_state > 0, need, drop = FALSE]
  if (anyDuplicated(z$cell_id)) stop("Earth slice has duplicate cell IDs at ", tm, " Ma.")
  z$land_area_km2 <- suppressWarnings(as.numeric(z$land_area_km2))
  z$land_area_km2[!is.finite(z$land_area_km2) | z$land_area_km2 <= 0] <- 1
  z
}
nearest <- function(query_lon, query_lat, ref_lon, ref_lat) {
  ans <- RANN::nn2(data = unit_xyz(ref_lon, ref_lat),
                   query = unit_xyz(query_lon, query_lat), k = 1L)
  list(index = as.integer(ans$nn.idx[, 1L]), distance_km = chord_to_km(ans$nn.dists[, 1L]))
}

tables_dir <- safe_dir(file.path(output, "transport_tables"))
audit_dir <- safe_dir(file.path(output, "audit"))
rows <- vector("list", length(times) - 1L)
for (ii in seq_len(length(times) - 1L)) {
  from_ma <- times[[ii]]
  to_ma <- times[[ii + 1L]]
  source <- read_slice(from_ma)
  target <- read_slice(to_ma)
  tr_from <- track_at_time(from_ma)
  tr_to <- track_at_time(to_ma)
  # Track IDs must retain one-to-one identity across the two time points.
  tr_to <- tr_to[match(tr_from$track_id, tr_to$track_id), , drop = FALSE]
  if (anyNA(tr_to$track_id)) stop("PALEOMAP track IDs are not stable across times.")
  # Semi-Lagrangian transport: start from each younger *target* land cell,
  # associate it with the nearest reconstructed H3 carrier at the younger
  # time, then look up where that same carrier was at the source time. A
  # target with no sufficiently close old land source is newly available land
  # and deliberately receives zero inherited occupancy.
  hit_track <- nearest(target$lon, target$lat, tr_to$lon, tr_to$lat)
  source_lon <- tr_from$lon[hit_track$index]
  source_lat <- tr_from$lat[hit_track$index]
  hit_source <- nearest(source_lon, source_lat, source$lon, source$lat)
  keep <- hit_track$distance_km <= max_match_distance_km &
    hit_source$distance_km <= max_match_distance_km
  tab <- data.frame(
    time_from_ma = from_ma,
    time_to_ma = to_ma,
    source_cell_id = as.character(source$cell_id[hit_source$index[keep]]),
    target_cell_id = as.character(target$cell_id[keep]),
    track_id = tr_from$track_id[hit_track$index[keep]],
    source_area_km2 = source$land_area_km2[hit_source$index[keep]],
    target_area_km2 = target$land_area_km2[keep],
    target_to_track_km = hit_track$distance_km[keep],
    track_to_source_km = hit_source$distance_km[keep],
    stringsAsFactors = FALSE
  )
  if (!nrow(tab)) stop("No target cells matched a valid old land source at ", from_ma, " -> ", to_ma, " Ma.")
  # Each target obtains one backward-tracked source state. This is an
  # interpolation of a probability field, not an expected-area sum, so the
  # appropriate target mixture weight is one rather than a source-area flux.
  tab$target_weight <- 1
  check <- hee_plate_grid_transport(
    tab, source_cells = source$cell_id, target_cells = target$cell_id
  )
  file <- file.path(tables_dir, paste0("paleomap_h3_transport_",
                                       ma_slug(from_ma), "_to_",
                                       ma_slug(to_ma), "Ma.rds"))
  saveRDS(check, file, compress = "gzip")
  map_source <- unique(tab$source_cell_id)
  map_target <- unique(tab$target_cell_id)
  rows[[ii]] <- data.frame(
    time_from_ma = from_ma,
    time_to_ma = to_ma,
    file = normalizePath(file, winslash = "/", mustWork = FALSE),
    n_source_land_cells = nrow(source),
    n_target_land_cells = nrow(target),
    n_mapped_source_cells = length(map_source),
    n_mapped_target_cells = length(map_target),
    source_coverage_fraction = length(map_source) / nrow(source),
    target_coverage_fraction = length(map_target) / nrow(target),
    target_to_track_km_q95 = as.numeric(stats::quantile(tab$target_to_track_km, .95)),
    track_to_source_km_q95 = as.numeric(stats::quantile(tab$track_to_source_km, .95)),
    target_to_track_km_max = max(tab$target_to_track_km),
    track_to_source_km_max = max(tab$track_to_source_km),
    transport_type = "PALEOMAP_H3_target_centred_nearest_track_carriage",
    scientific_boundary = paste(
      "External PALEOMAP target-centred track-based grid carriage;",
      "nearest-track approximation, not polygon-overlap or an organismal dispersal kernel."
    ),
    stringsAsFactors = FALSE
  )
  message("Built transport ", from_ma, " -> ", to_ma, " Ma: ",
          nrow(tab), " links; target coverage ",
          formatC(rows[[ii]]$target_coverage_fraction, digits = 3, format = "f"))
}
out_index <- do.call(rbind, rows)
write_csv(out_index, file.path(output, "paleomap_h3_transport_index.csv"))
write_csv(data.frame(
  parameter = c("plate_tracks", "palaeo_slice_index", "time_selection",
                "n_time_slices", "start_ma", "max_match_distance_km",
                "transport_formula", "scientific_boundary"),
  value = c(normalizePath(plate_tracks_path, winslash = "/"),
            normalizePath(slice_index_path, winslash = "/"),
            time_selection, length(times), start_ma, max_match_distance_km,
            "q_target=q_nearest_backward_tracked_source; one target-centred source per target grid cell",
            "Palaeogeographic carriage uses external PALEOMAP H3 reconstructed tracks; it is target-centred nearest-track interpolation, not active dispersal or polygon overlap."),
  stringsAsFactors = FALSE
), file.path(audit_dir, "transport_metadata.csv"))
write_csv(out_index[, c("time_from_ma", "time_to_ma", "n_source_land_cells",
                        "n_target_land_cells", "source_coverage_fraction",
                        "target_coverage_fraction", "target_to_track_km_q95",
                        "track_to_source_km_q95")],
          file.path(audit_dir, "transport_coverage_by_interval.csv"))
writeLines(c(
  "# Case05 PALEOMAP H3 grid transport",
  "",
  "Each interval table carries a source time-slice occupancy field to the next younger 1 degree palaeo-land slice before biological movement.",
  "The external PALEOMAP H3 reconstruction provides the track positions. Each younger land cell is followed backward to one older source cell along its nearest H3 track.",
  "This target-centred method avoids artificial target-grid gaps when H3 is coarser than the 1 degree grid. It is a documented nearest-track approximation, not a GPlates polygon-overlap transport matrix. It must not be described as organismal dispersal or as exact land-area overlap."
), file.path(output, "README.md"), useBytes = TRUE)
cat("DONE\nOutput:", output, "\n")
