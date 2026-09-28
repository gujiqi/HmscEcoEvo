#' Build a master deep-time axis with exact lineage and observation events
#'
#' Environmental rasters are often available only at coarse time slices, while
#' speciation nodes and fossil age boundaries occur between those slices.  This
#' helper keeps the two roles separate: it inserts exact biological and
#' observational events into the process time axis, and records the enclosing
#' palaeoenvironment interval.  A runner must declare how it represents Earth
#' covariates within that interval; this function does not silently interpolate
#' land masks or plate geometry.
#'
#' Larger `time_ma` values are older.  If a dated tree is supplied, node ages
#' are inferred from its branch lengths and only times from the root to the
#' present are returned.
#'
#' @param env_times Numeric palaeoenvironment slice ages in Ma.
#' @param tree Optional ultrametric `ape::phylo` tree with branch lengths in
#'   Ma.
#' @param fossil_times Optional numeric fossil ages or age boundaries in Ma.
#' @param additional_times Optional numeric geological or user-defined event
#'   ages in Ma.
#' @param tolerance Numeric matching tolerance in Ma.
#' @return A data frame with `time_ma`, event flags and the enclosing older and
#'   younger environmental slices.
#' @export
hee_master_time_axis <- function(env_times,
                                 tree = NULL,
                                 fossil_times = NULL,
                                 additional_times = NULL,
                                 tolerance = 1e-8) {
  env <- sort(unique(suppressWarnings(as.numeric(env_times))), decreasing = TRUE)
  if (length(env) < 2L || any(!is.finite(env))) {
    stop("env_times must contain at least two finite distinct ages in Ma.",
         call. = FALSE)
  }
  tolerance <- suppressWarnings(as.numeric(tolerance)[1L])
  if (!is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be a non-negative finite number.", call. = FALSE)
  }
  tree_times <- numeric()
  root_age <- max(env)
  if (!is.null(tree)) {
    if (!requireNamespace("ape", quietly = TRUE) || !inherits(tree, "phylo") ||
        is.null(tree$edge.length)) {
      stop("tree must be a dated ape::phylo object with branch lengths in Ma.",
           call. = FALSE)
    }
    if (!isTRUE(ape::is.ultrametric(tree))) {
      stop("tree must be ultrametric for a common Ma time axis.", call. = FALSE)
    }
    tree_times <- as.numeric(ape::branching.times(tree))
    root_age <- max(tree_times)
  }
  fossil_times <- suppressWarnings(as.numeric(fossil_times))
  additional_times <- suppressWarnings(as.numeric(additional_times))
  extra <- c(tree_times, fossil_times[is.finite(fossil_times)],
             additional_times[is.finite(additional_times)])
  all_times <- sort(unique(c(env, extra)), decreasing = TRUE)
  all_times <- all_times[all_times <= root_age + tolerance & all_times >= -tolerance]
  all_times[abs(all_times) <= tolerance] <- 0
  if (!any(abs(all_times - root_age) <= tolerance)) all_times <-
    sort(c(all_times, root_age), decreasing = TRUE)

  has_time <- function(values, value) any(abs(values - value) <= tolerance)
  bracket <- lapply(all_times, function(tm) {
    exact <- which(abs(env - tm) <= tolerance)
    if (length(exact)) return(c(env[[exact[[1L]]]], env[[exact[[1L]]]], TRUE))
    older <- env[env > tm]
    younger <- env[env < tm]
    if (!length(older) || !length(younger)) return(c(NA_real_, NA_real_, FALSE))
    c(min(older), max(younger), FALSE)
  })
  bracket <- do.call(rbind, bracket)
  out <- data.frame(
    time_ma = all_times,
    is_environment_slice = vapply(all_times, function(tm) has_time(env, tm), logical(1)),
    is_tree_node = vapply(all_times, function(tm) has_time(tree_times, tm), logical(1)),
    is_fossil_event = vapply(all_times, function(tm) has_time(fossil_times, tm), logical(1)),
    is_additional_event = vapply(all_times, function(tm) has_time(additional_times, tm), logical(1)),
    earth_slice_older_ma = as.numeric(bracket[, 1L]),
    earth_slice_younger_ma = as.numeric(bracket[, 2L]),
    earth_slice_exact = as.logical(bracket[, 3L]),
    stringsAsFactors = FALSE
  )
  out$earth_representation <- ifelse(
    out$earth_slice_exact, "exact_environment_slice",
    "runner_must_declare_piecewise_or_interpolated_earth_state"
  )
  out
}
