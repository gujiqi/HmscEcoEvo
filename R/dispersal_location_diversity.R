#' Calibrate terminal range-support scales for location-density maps
#'
#' Converts a normalized lineage-location density into a terminal-calibrated
#' **range-support** surface. It is useful when an explicit phylogeographic
#' diffusion model returns one normalized location density per lineage, while
#' present-day occurrence grids provide a defensible reference for each tip's
#' spatial support size. This is not a demographic occupancy model.
#'
#' For lineage \eqn{l}, location density \eqn{pi_lc}, and scale
#' \eqn{kappa_l}, the support transform is
#' \deqn{a_lc = min(1, kappa_l pi_lc).}
#' `kappa_l` is chosen so that the sum of \eqn{a_lc} approaches a declared
#' fraction of the terminal observed grid count. The fraction is below one so
#' that a finite calibration exists even when a terminal density has support
#' only in its observed cells.
#'
#' The resulting quantities are explicitly called *range support*, not
#' occupancy probabilities, species richness, colonisation, persistence, or
#' extinction. Use an independently calibrated dynamic range/occupancy model
#' before making those claims.
#'
#' @param location_density Numeric state-by-lineage matrix. Each column must
#'   be finite, non-negative, and have positive mass; it is internally
#'   normalized to one.
#' @param terminal_support_cells Named positive numeric vector with one
#'   observed terminal grid count per lineage.
#' @param target_fraction A number in `(0, 1)`. The calibration target is this
#'   fraction times each observed terminal grid count.
#' @param tolerance Numerical tolerance for the monotone root solve.
#' @return A data frame with lineage-specific support scales and achieved
#'   terminal support totals.
#' @export
hee_dispersal_range_support_calibrate <- function(location_density,
                                                   terminal_support_cells,
                                                   target_fraction = 0.995,
                                                   tolerance = 1e-8) {
  density <- as.matrix(location_density)
  storage.mode(density) <- "double"
  if (is.null(colnames(density)) || anyDuplicated(colnames(density)) ||
      !ncol(density) || !nrow(density)) {
    stop("location_density must be a non-empty matrix with unique lineage columns.",
         call. = FALSE)
  }
  total <- colSums(density)
  if (any(!is.finite(density) | density < 0) || any(!is.finite(total) | total <= 0)) {
    stop("location_density must be finite, non-negative, and positive by lineage.",
         call. = FALSE)
  }
  density <- sweep(density, 2L, total, "/")
  if (!is.numeric(target_fraction) || length(target_fraction) != 1L ||
      !is.finite(target_fraction) || target_fraction <= 0 || target_fraction >= 1) {
    stop("target_fraction must be one finite number strictly between zero and one.",
         call. = FALSE)
  }
  support <- suppressWarnings(as.numeric(terminal_support_cells))
  names(support) <- names(terminal_support_cells)
  if (is.null(names(support)) || anyDuplicated(names(support)) ||
      !all(colnames(density) %in% names(support))) {
    stop("terminal_support_cells must be uniquely named for every lineage column.",
         call. = FALSE)
  }
  support <- support[colnames(density)]
  nonzero <- colSums(density > 0)
  if (any(!is.finite(support) | support <= 0 | support > nonzero)) {
    stop("terminal_support_cells must be positive and no larger than each density support.",
         call. = FALSE)
  }
  target <- target_fraction * support
  solve_one <- function(probability, target_value) {
    objective <- function(scale) sum(pmin(1, scale * probability)) - target_value
    upper <- 1
    while (objective(upper) < 0) upper <- upper * 2
    stats::uniroot(objective, lower = 0, upper = upper, tol = tolerance)$root
  }
  scale <- vapply(seq_len(ncol(density)), function(i) {
    solve_one(density[, i], target[[i]])
  }, numeric(1))
  achieved <- vapply(seq_len(ncol(density)), function(i) {
    sum(pmin(1, scale[[i]] * density[, i]))
  }, numeric(1))
  data.frame(
    lineage = colnames(density),
    observed_terminal_support_cells = support,
    target_fraction = target_fraction,
    target_terminal_support_cells = target,
    support_scale = scale,
    achieved_terminal_support_cells = achieved,
    stringsAsFactors = FALSE
  )
}

#' Summarise dispersal-conditioned lineage-support diversity
#'
#' Applies terminal-calibrated range-support scales to a matrix of explicit
#' lineage-location densities and returns non-demographic diversity summaries.
#' These summaries answer how many modelled lineage-support fields overlap in a
#' cell under a chosen diffusion model. They do **not** estimate ecological
#' occupancy, historical species richness, colonisation, persistence, local
#' extinction, or census diversity.
#'
#' The support intensity is \eqn{I_c = sum_l a_lc}. The Shannon mixture is
#' \eqn{H_c = -sum_l w_lc log(w_lc)}, where
#' \eqn{w_lc = a_lc / I_c}; `effective_lineages` is \eqn{exp(H_c)}.
#'
#' @param location_density Numeric state-by-lineage location-density matrix.
#' @param support_scale Named positive vector calibrated by
#'   [hee_dispersal_range_support_calibrate()].
#' @param return_support Should the state-by-lineage support matrix be returned?
#' @return A list with `summary` and optionally `support`.
#' @export
hee_dispersal_range_support_diversity <- function(location_density,
                                                   support_scale,
                                                   return_support = FALSE) {
  density <- as.matrix(location_density)
  storage.mode(density) <- "double"
  if (is.null(rownames(density)) || is.null(colnames(density)) ||
      anyDuplicated(rownames(density)) || anyDuplicated(colnames(density))) {
    stop("location_density must have unique named state rows and lineage columns.",
         call. = FALSE)
  }
  total <- colSums(density)
  if (any(!is.finite(density) | density < 0) || any(!is.finite(total) | total <= 0)) {
    stop("location_density must be finite, non-negative, and positive by lineage.",
         call. = FALSE)
  }
  if (is.null(names(support_scale)) || !all(colnames(density) %in% names(support_scale))) {
    stop("support_scale must be a named vector covering every lineage column.",
         call. = FALSE)
  }
  scale <- suppressWarnings(as.numeric(support_scale[colnames(density)]))
  if (any(!is.finite(scale) | scale <= 0)) {
    stop("support_scale must be finite and strictly positive.", call. = FALSE)
  }
  density <- sweep(density, 2L, total, "/")
  support <- sweep(density, 2L, scale, "*")
  support[] <- pmin(1, support)
  intensity <- rowSums(support)
  mixture <- support / ifelse(intensity > 0, intensity, 1)
  mixture[intensity <= 0, ] <- 0
  log_mixture <- matrix(0, nrow = nrow(mixture), ncol = ncol(mixture))
  positive <- mixture > 0
  log_mixture[positive] <- mixture[positive] * log(mixture[positive])
  shannon <- -rowSums(log_mixture)
  simpson <- 1 - rowSums(mixture^2)
  effective <- exp(shannon)
  effective[intensity <= 0] <- NA_real_
  summary <- data.frame(
    state = rownames(density),
    terminal_calibrated_lineage_support_intensity = intensity,
    support_mixture_shannon = shannon,
    support_mixture_simpson = simpson,
    support_mixture_effective_lineages = effective,
    stringsAsFactors = FALSE
  )
  output <- list(
    summary = summary,
    interpretation = paste(
      "Dispersal-conditioned, terminal-calibrated lineage-support diversity.",
      "It is not posterior occupancy, historical species richness, or ecological alpha diversity."
    )
  )
  if (isTRUE(return_support)) output$support <- support
  output
}
