#' Build an effort-aware target-group HMSC training design
#'
#' Aggregates focal-taxon occurrence records and a modern terrestrial grid to a
#' regular geographic grid.  It is designed for data such as GBIF occurrences
#' where non-records are **not** confirmed absences.  The returned response
#' matrix is consequently labelled a target-group presence-background design:
#' a zero means that a focal species was not recorded in a grid cell with at
#' least the declared target-group recording intensity.  It must not be
#' interpreted as a detection-corrected absence.
#'
#' The target-group record count is returned as an observation covariate.  It
#' may be used in a modern HMSC sensitivity model, but it is not an ecological
#' predictor and must not be projected into palaeoenvironmental time slices.
#'
#' @param sites Data frame identifying the original response rows.  It needs
#'   `site_id_col`, `lon_col`, and `lat_col`.
#' @param Y Numeric/logical site-by-species occurrence matrix.  Its row names
#'   or its `site_id_col` attribute must match `sites`.
#' @param grid Complete modern land grid, including coordinates.  All its cells
#'   are retained in `grid_cells`, even when not eligible for model fitting.
#' @param records Optional focal-taxon occurrence table.  It needs longitude,
#'   latitude and a species column.  When omitted, target-group effort is
#'   approximated by focal-species incidences in `Y` and labelled accordingly.
#' @param grid_deg Positive regular grid resolution in degrees.
#' @param min_target_group_records Minimum number of focal-taxon records needed
#'   for a cell to enter the target-group presence-background training design.
#' @param site_id_col,lon_col,lat_col Column names for site/grid coordinates.
#' @param record_species_col,record_lon_col,record_lat_col Record columns.
#'
#' @return A list with `Y`, eligible `sites`, complete `grid_cells`, a
#'   `cell_map`, and `metadata`.
#' @export
#'
#' @examples
#' sites <- data.frame(site_id = c("a", "b"), lon = c(1, 7), lat = c(1, 1))
#' y <- matrix(c(1, 0, 0, 1), 2, 2,
#'             dimnames = list(c("a", "b"), c("sp1", "sp2")))
#' grid <- data.frame(site_id = c("a", "b", "c"), lon = c(1, 7, 11),
#'                    lat = c(1, 1, 1))
#' hee_hmsc_target_group_design(sites, y, grid, grid_deg = 4,
#'                              min_target_group_records = 1)
hee_hmsc_target_group_design <- function(
    sites,
    Y,
    grid,
    records = NULL,
    grid_deg = 4,
    min_target_group_records = 5L,
    site_id_col = "site_id",
    lon_col = "lon",
    lat_col = "lat",
    record_species_col = "selected_species_id",
    record_lon_col = "lon",
    record_lat_col = "lat") {
  sites <- as.data.frame(sites, stringsAsFactors = FALSE)
  grid <- as.data.frame(grid, stringsAsFactors = FALSE)
  .require_cols(sites, c(site_id_col, lon_col, lat_col), "sites")
  .require_cols(grid, c(site_id_col, lon_col, lat_col), "grid")
  y <- as.matrix(Y)
  storage.mode(y) <- "double"
  if (!is.numeric(y) || !nrow(y) || !ncol(y) ||
      any(!is.finite(y)) || any(!(y %in% c(0, 1)))) {
    stop("Y must be a non-empty finite binary site-by-species matrix.",
         call. = FALSE)
  }
  site_ids <- as.character(sites[[site_id_col]])
  if (any(!nzchar(site_ids)) || anyDuplicated(site_ids)) {
    stop("sites must have unique non-blank site IDs.", call. = FALSE)
  }
  if (!is.null(rownames(y)) && all(nzchar(rownames(y)))) {
    missing_y <- setdiff(site_ids, rownames(y))
    if (length(missing_y)) {
      stop("Y row names are missing site IDs: ",
           paste(utils::head(missing_y, 5L), collapse = ", "), call. = FALSE)
    }
    y <- y[site_ids, , drop = FALSE]
  } else if (nrow(y) != nrow(sites)) {
    stop("Y must have sites rows or be named by site_id.", call. = FALSE)
  }
  if (is.null(colnames(y)) || any(!nzchar(colnames(y))) ||
      anyDuplicated(colnames(y))) {
    stop("Y must have unique non-blank species column names.", call. = FALSE)
  }
  grid_deg <- suppressWarnings(as.numeric(grid_deg)[1L])
  min_target_group_records <- suppressWarnings(
    as.integer(min_target_group_records)[1L]
  )
  if (!is.finite(grid_deg) || grid_deg <= 0 || grid_deg > 180 ||
      !is.finite(min_target_group_records) || min_target_group_records < 1L) {
    stop("grid_deg must be in (0, 180] and min_target_group_records >= 1.",
         call. = FALSE)
  }
  if (abs(360 / grid_deg - round(360 / grid_deg)) > 1e-8 ||
      abs(180 / grid_deg - round(180 / grid_deg)) > 1e-8) {
    stop("grid_deg must divide both 360 and 180 exactly.", call. = FALSE)
  }
  n_lon <- as.integer(round(360 / grid_deg))
  n_lat <- as.integer(round(180 / grid_deg))
  bin_coordinates <- function(lon, lat, label) {
    lon <- suppressWarnings(as.numeric(lon))
    lat <- suppressWarnings(as.numeric(lat))
    if (any(!is.finite(lon) | !is.finite(lat) | lon < -180 | lon > 180 |
            lat < -90 | lat > 90)) {
      stop(label, " has invalid longitude/latitude values.", call. = FALSE)
    }
    lon_bin <- pmin(n_lon - 1L, pmax(0L, floor((lon + 180) / grid_deg)))
    lat_bin <- pmin(n_lat - 1L, pmax(0L, floor((lat + 90) / grid_deg)))
    data.frame(
      grid_lon_index = as.integer(lon_bin),
      grid_lat_index = as.integer(lat_bin),
      grid_cell_id = sprintf("grid%03d_%03d", lon_bin + 1L, lat_bin + 1L),
      grid_lon = -180 + (lon_bin + 0.5) * grid_deg,
      grid_lat = -90 + (lat_bin + 0.5) * grid_deg,
      stringsAsFactors = FALSE
    )
  }
  site_bin <- bin_coordinates(sites[[lon_col]], sites[[lat_col]], "sites")
  grid_bin <- bin_coordinates(grid[[lon_col]], grid[[lat_col]], "grid")
  grid_cells <- grid_bin[!duplicated(grid_bin$grid_cell_id), , drop = FALSE]
  grid_count <- tabulate(match(grid_cells$grid_cell_id, grid_bin$grid_cell_id),
                         nbins = nrow(grid_cells))
  grid_cells$original_grid_cells <- as.integer(grid_count)
  rownames(grid_cells) <- NULL

  aggregate_binary <- function(values, keys, levels) {
    out <- matrix(0, nrow = length(levels), ncol = ncol(values),
                  dimnames = list(levels, colnames(values)))
    for (i in seq_along(levels)) {
      hit <- which(keys == levels[[i]])
      if (length(hit)) out[i, ] <- as.numeric(colSums(values[hit, , drop = FALSE]) > 0)
    }
    out
  }
  cell_levels <- as.character(grid_cells$grid_cell_id)
  y_grid <- aggregate_binary(y, site_bin$grid_cell_id, cell_levels)

  if (is.null(records)) {
    effort <- rowSums(y_grid)
    effort_source <- "focal_species_incidence_count_fallback"
  } else {
    records <- as.data.frame(records, stringsAsFactors = FALSE)
    .require_cols(records, c(record_species_col, record_lon_col, record_lat_col),
                  "records")
    record_species <- as.character(records[[record_species_col]])
    keep <- !is.na(record_species) & record_species %in% colnames(y_grid)
    if (!any(keep)) {
      stop("records contains no species that match Y column names.", call. = FALSE)
    }
    record_bin <- bin_coordinates(records[[record_lon_col]][keep],
                                  records[[record_lat_col]][keep], "records")
    effort <- tabulate(match(record_bin$grid_cell_id, cell_levels),
                       nbins = length(cell_levels))
    effort_source <- "focal_taxon_occurrence_record_count"
  }
  grid_cells$target_group_record_count <- as.integer(effort)
  grid_cells$target_group_species_count <- rowSums(y_grid)
  grid_cells$log_target_group_record_count <- log1p(grid_cells$target_group_record_count)
  grid_cells$eligible_target_group_background <-
    grid_cells$target_group_record_count >= min_target_group_records
  eligible <- which(grid_cells$eligible_target_group_background)
  if (length(eligible) < 2L || any(colSums(y_grid[eligible, , drop = FALSE]) == 0)) {
    stop("The selected effort threshold leaves too few training cells or omits at least one focal species.",
         call. = FALSE)
  }
  site_to_grid <- data.frame(
    original_site_id = site_ids,
    grid_cell_id = site_bin$grid_cell_id,
    stringsAsFactors = FALSE
  )
  result_sites <- grid_cells[eligible, c("grid_cell_id", "grid_lon", "grid_lat",
                                         "target_group_record_count",
                                         "target_group_species_count",
                                         "log_target_group_record_count"),
                            drop = FALSE]
  names(result_sites)[names(result_sites) == "grid_cell_id"] <- site_id_col
  names(result_sites)[names(result_sites) == "grid_lon"] <- lon_col
  names(result_sites)[names(result_sites) == "grid_lat"] <- lat_col
  list(
    Y = y_grid[eligible, , drop = FALSE],
    sites = result_sites,
    grid_cells = grid_cells,
    cell_map = site_to_grid,
    metadata = data.frame(
      field = c("observation_design", "zero_semantics", "target_group_effort_source",
                "grid_deg", "min_target_group_records", "n_complete_grid_cells",
                "n_training_cells", "n_species"),
      value = c(
        "target_group_presence_background",
        "focal_species_not_recorded_in_target_group_effort_qualified_cell_not_confirmed_absence",
        effort_source, grid_deg, min_target_group_records, nrow(grid_cells),
        length(eligible), ncol(y_grid)
      ),
      stringsAsFactors = FALSE
    )
  )
}
