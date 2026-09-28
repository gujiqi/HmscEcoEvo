#' Write a cell-time palaeo grid as independently readable time slices
#'
#' Converts a complete palaeoenvironment or palaeogeography table into one RDS
#' file per `time_ma`.  The helper is intended for large 1-degree time cubes:
#' downstream dynamic runs can read a complete valid grid for one time slice,
#' update occupancy, write its outputs, and release the slice instead of
#' holding the entire 540--0 Ma cube in every worker process.
#'
#' This is an I/O operation only. It never samples cells, changes CRS, rescales
#' environmental variables, or alters the values in `x`.
#'
#' @param x A data frame with a time column and one row per cell-time record.
#' @param output_dir Output directory for RDS slices and their CSV index.
#' @param time_col Name of the numeric time column, normally `"time_ma"`.
#' @param keep_cols Optional columns to retain. `time_col` is always retained.
#' @param key_cols Columns that must be unique within each time slice. Use
#'   `"cell_id"` when a fixed grid identifier is available.
#' @param prefix Filename prefix for time-slice RDS files.
#' @param compression RDS compression passed to [saveRDS()].
#' @param overwrite Logical; allow replacement of existing slice files.
#' @return A data frame containing `time_ma`, `file`, `n_rows`, `n_columns`,
#'   and `retained_columns`. The same index is written as
#'   `<prefix>_index.csv`.
#' @export
#' @examples
#' cube <- data.frame(
#'   cell_id = rep(c("c1", "c2"), 2), time_ma = rep(c(10, 0), each = 2),
#'   temperature = c(20, 21, 18, 19)
#' )
#' out <- tempfile("hee_slices_")
#' hee_write_timecube_slices(cube, out, keep_cols = names(cube))
hee_write_timecube_slices <- function(x,
                                      output_dir,
                                      time_col = "time_ma",
                                      keep_cols = NULL,
                                      key_cols = "cell_id",
                                      prefix = "palaeo_grid",
                                      compression = "gzip",
                                      overwrite = FALSE) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!time_col %in% names(x)) {
    stop("x is missing time_col: ", time_col, call. = FALSE)
  }
  x[[time_col]] <- suppressWarnings(as.numeric(x[[time_col]]))
  if (any(!is.finite(x[[time_col]]))) {
    stop("x contains non-finite values in ", time_col, ".", call. = FALSE)
  }
  if (is.null(keep_cols)) keep_cols <- names(x)
  keep_cols <- unique(c(time_col, as.character(keep_cols)))
  missing_keep <- setdiff(keep_cols, names(x))
  if (length(missing_keep)) {
    stop("keep_cols are missing from x: ", paste(missing_keep, collapse = ", "),
         call. = FALSE)
  }
  key_cols <- as.character(key_cols)
  missing_key <- setdiff(key_cols, keep_cols)
  if (length(missing_key)) {
    stop("key_cols must be retained in keep_cols: ",
         paste(missing_key, collapse = ", "), call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  times <- sort(unique(x[[time_col]]), decreasing = TRUE)
  slug <- function(z) gsub("\\.", "p", formatC(z, format = "f", digits = 3,
                                                    drop0trailing = TRUE))
  index <- vector("list", length(times))
  for (i in seq_along(times)) {
    tm <- times[[i]]
    slice <- x[x[[time_col]] == tm, keep_cols, drop = FALSE]
    if (anyDuplicated(slice[, key_cols, drop = FALSE])) {
      stop("Duplicate key(s) in time slice ", tm, " for ",
           paste(key_cols, collapse = ", "), ".", call. = FALSE)
    }
    path <- file.path(output_dir, paste0(prefix, "_", slug(tm), "Ma.rds"))
    if (file.exists(path) && !isTRUE(overwrite)) {
      stop("Slice already exists: ", path,
           ". Set overwrite = TRUE to replace it.", call. = FALSE)
    }
    saveRDS(slice, path, compress = compression)
    index[[i]] <- data.frame(
      time_ma = tm,
      file = normalizePath(path, winslash = "/", mustWork = TRUE),
      n_rows = nrow(slice),
      n_columns = ncol(slice),
      retained_columns = paste(keep_cols, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  index <- do.call(rbind, index)
  utils::write.csv(index, file.path(output_dir, paste0(prefix, "_index.csv")),
                   row.names = FALSE, fileEncoding = "UTF-8")
  index
}
