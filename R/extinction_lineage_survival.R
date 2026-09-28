#' Condition extant-tree lineages on survival somewhere in the available world
#'
#' An extant-only dated tree establishes that every active branch leading to a
#' sampled tip survived until its next observed node (or to the present). A
#' cell-level mean-field occupancy update, however, can let the probability of
#' *all* cells occupied by such a branch decay towards zero. This helper
#' conditions each selected lineage's cell probabilities on the event that at
#' least one available cell remains occupied:
#'
#' \deqn{P(Z_c=1 \mid \sum_u Z_u > 0) =
#' P(Z_c=1) / [1 - \prod_u (1-P(Z_u=1))].}
#'
#' It preserves the relative spatial weights produced by environmental
#' filtering, dispersal, colonisation and local persistence, but prevents an
#' impossible global loss of a branch known from the dated extant tree to have
#' descendants. It is **not** an estimate of historical global extinction and
#' it must not be applied to an extinct-tip or fossilized-birth-death analysis.
#'
#' This is a **conditional marginal post-processing operation**, not a forward
#' dynamic-occupancy transition.  In particular, applying it after every
#' numerical time step changes the transition kernel by reinflating any branch
#' whose total occupancy becomes small.  Use it only after posterior particles
#' or complete trajectories have been weighted by an explicitly defined
#' survival/endpoint event.  A forward process scenario must retain the
#' unconditioned state and record a zero-survival particle as incompatible at
#' its node or endpoint, rather than silently renormalising it.
#'
#' If a selected lineage has exactly zero support across all available cells,
#' the current Earth/process scenario is incompatible with the required tree
#' survival condition and the function stops rather than inventing a range.
#'
#' @param occupancy Numeric cell-by-lineage matrix of finite probabilities in
#'   `[0, 1]`.
#' @param condition_lineage Logical vector, or lineage names, selecting columns
#'   to condition. The default conditions every column.
#' @param zero_tolerance Non-negative tolerance below which a lineage has no
#'   numerical support and is treated as incompatible.
#' @return A list with conditioned `occupancy` and an `audit` table containing
#'   pre-conditioning survival probabilities and expected occupancy masses.
#' @export
#'
#' @examples
#' q <- cbind(a = c(0.1, 0.2), b = c(0.5, 0))
#' hee_extinction_condition_extant_lineage_survival(q)$audit
hee_extinction_condition_extant_lineage_survival <- function(
    occupancy,
    condition_lineage = NULL,
    zero_tolerance = 1e-300) {
  if (is.data.frame(occupancy)) occupancy <- as.matrix(occupancy)
  if (!is.matrix(occupancy) || !is.numeric(occupancy) || !nrow(occupancy) || !ncol(occupancy)) {
    stop("occupancy must be a non-empty numeric cell-by-lineage matrix.", call. = FALSE)
  }
  if (any(!is.finite(occupancy)) || any(occupancy < 0 | occupancy > 1)) {
    stop("occupancy must contain finite probabilities in [0, 1].", call. = FALSE)
  }
  if (length(zero_tolerance) != 1L || !is.finite(zero_tolerance) || zero_tolerance < 0) {
    stop("zero_tolerance must be one finite non-negative number.", call. = FALSE)
  }
  if (is.null(condition_lineage)) {
    selected <- rep(TRUE, ncol(occupancy))
  } else if (is.character(condition_lineage)) {
    if (is.null(colnames(occupancy))) {
      stop("Named condition_lineage requires occupancy column names.", call. = FALSE)
    }
    unknown <- setdiff(condition_lineage, colnames(occupancy))
    if (length(unknown)) {
      stop("condition_lineage names are not occupancy columns: ",
           paste(utils::head(unknown, 10L), collapse = ", "), call. = FALSE)
    }
    selected <- colnames(occupancy) %in% condition_lineage
  } else if (is.logical(condition_lineage) && length(condition_lineage) == ncol(occupancy) && !anyNA(condition_lineage)) {
    selected <- condition_lineage
  } else {
    stop("condition_lineage must be NULL, column names, or a non-missing logical vector with one entry per lineage.", call. = FALSE)
  }

  out <- occupancy
  lineage <- if (is.null(colnames(out))) paste0("lineage_", seq_len(ncol(out))) else colnames(out)
  audit <- data.frame(
    lineage = lineage,
    conditioned = selected,
    expected_occupancy_mass_before = colSums(out),
    survival_probability_before = NA_real_,
    expected_occupancy_mass_after = NA_real_,
    stringsAsFactors = FALSE
  )
  for (j in which(selected)) {
    p <- out[, j]
    # log1p retains useful precision for very small p, unlike prod(1 - p).
    log_none <- sum(log1p(-p))
    p_any <- -expm1(log_none)
    audit$survival_probability_before[[j]] <- p_any
    if (!is.finite(p_any) || p_any <= zero_tolerance) {
      stop("Lineage '", lineage[[j]], "' has no surviving cell support. ",
           "The Earth/process scenario is incompatible with the required extant-tree branch survival.",
           call. = FALSE)
    }
    out[, j] <- pmin(1, p / p_any)
  }
  audit$expected_occupancy_mass_after <- colSums(out)
  audit$survival_probability_before[!selected] <- 1 - exp(colSums(log1p(-out[, !selected, drop = FALSE])))
  list(occupancy = out, audit = audit)
}
