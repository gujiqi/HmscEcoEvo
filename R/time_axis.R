#' Read the full time axis from a deep-time environmental cube
#'
#' @param env_cube A time-cube object, compact RDS path, or object accepted by
#'   [hee_load_timecube()].
#' @param decreasing Return older-to-younger ages when `TRUE`.
#' @return Numeric ages in Ma.
#' @export
hee_time_axis_from_env_cube <- function(env_cube, decreasing = TRUE) {
  cube <- .hee_as_timecube(env_cube)
  times <- unique(as.numeric(cube$times))
  sort(times, decreasing = decreasing)
}

#' Validate a full deep-time time axis
#'
#' @param times Numeric ages in Ma.
#' @param quick Logical; skip strict full-mode checks when `TRUE`.
#' @param min_time_slices Minimum number of unique time slices required in full mode.
#' @param required_range Numeric vector of required endpoint ages.
#' @param tolerance Numeric tolerance for endpoint matching.
#' @return Invisibly returns sorted unique `times`.
#' @export
hee_assert_full_time_axis <- function(times,
                                      quick = FALSE,
                                      min_time_slices = 3,
                                      required_range = c(540, 0),
                                      tolerance = 1e-8) {
  times <- sort(unique(as.numeric(times)), decreasing = TRUE)
  if (!quick) {
    if (length(times) < min_time_slices) {
      stop("Full mode requires at least ", min_time_slices,
           " time slices; got ", length(times),
           ". Use quick = TRUE for toy/debug runs.", call. = FALSE)
    }
    missing_endpoints <- required_range[!vapply(required_range, function(z) {
      any(abs(times - z) <= tolerance)
    }, logical(1))]
    if (length(missing_endpoints) > 0) {
      stop("Full mode time axis is missing required endpoint age(s): ",
           paste(missing_endpoints, collapse = ", "), " Ma.", call. = FALSE)
    }
  }
  invisible(times)
}
