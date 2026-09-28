#' Anchor a forward occupancy prediction to an endpoint nowcast
#'
#' Form a transparent terminal model ensemble:
#' 
#' `q_anchor = (1 - w) q_forward + w q_endpoint`.
#'
#' This helper is deliberately **not** a historical process equation and does
#' not make an endpoint model independent evidence. It is useful when a
#' forward deep-time scenario must be reconciled with an explicitly declared
#' independent or held-out endpoint nowcast. If the endpoint nowcast is made
#' from the same records used to fit HMSC, use it only as a training diagnostic
#' or modular sensitivity layer; do not describe the anchored result as an
#' independently validated palaeodistribution reconstruction.
#'
#' @param forward_occupancy Numeric probability vector or cell-by-lineage
#'   matrix from a forward dynamic model.
#' @param endpoint_nowcast Numeric vector or matrix with identical dimensions,
#'   containing a 0 Ma endpoint nowcast.
#' @param anchor_weight Weight in `[0, 1]` assigned to the endpoint nowcast.
#'   `0` returns the forward scenario and `1` returns the endpoint nowcast.
#' @return A numeric vector or matrix of anchored probabilities in `[0, 1]`.
#' @export
#'
#' @examples
#' forward <- matrix(c(0.1, 0.8), ncol = 1)
#' nowcast <- matrix(c(0.3, 0.6), ncol = 1)
#' hee_endpoint_nowcast_anchor(forward, nowcast, anchor_weight = 0.5)
hee_endpoint_nowcast_anchor <- function(forward_occupancy,
                                        endpoint_nowcast,
                                        anchor_weight = 0.5) {
  weight <- suppressWarnings(as.numeric(anchor_weight)[1L])
  if (!is.finite(weight) || weight < 0 || weight > 1) {
    stop("anchor_weight must be one finite value in [0, 1].", call. = FALSE)
  }
  forward <- as.matrix(forward_occupancy)
  endpoint <- as.matrix(endpoint_nowcast)
  if (!identical(dim(forward), dim(endpoint))) {
    stop("forward_occupancy and endpoint_nowcast must have identical dimensions.",
         call. = FALSE)
  }
  storage.mode(forward) <- "double"
  storage.mode(endpoint) <- "double"
  valid <- is.finite(forward) & is.finite(endpoint) &
    forward >= 0 & forward <= 1 & endpoint >= 0 & endpoint <= 1
  if (!all(valid)) {
    stop("forward_occupancy and endpoint_nowcast must be finite probabilities in [0, 1].",
         call. = FALSE)
  }
  out <- (1 - weight) * forward + weight * endpoint
  dimnames(out) <- dimnames(forward)
  if (is.null(dim(forward_occupancy))) as.numeric(out) else out
}

#' Choose the minimum endpoint-anchor weight that meets a prevalence tolerance
#'
#' This helper selects the smallest weight for
#' [hee_endpoint_nowcast_anchor()] that brings the mean probability of a
#' forward endpoint into a declared tolerance interval relative to an endpoint
#' nowcast. It is intended for an explicitly declared endpoint *calibration*
#' subset, preferably spatially held out from the HMSC fit. When the same
#' records were used to fit HMSC, it is only a training diagnostic and terminal
#' model-ensemble setting, not independent validation or particle smoothing.
#'
#' The selected weight is analytic because the anchored mean is a convex
#' combination of the forward and endpoint means. It does not alter any
#' historical state before the endpoint.
#'
#' @param forward_occupancy Numeric probability vector or matrix from a
#'   forward model.
#' @param endpoint_nowcast Numeric probability vector or matrix of identical
#'   dimensions containing the endpoint nowcast.
#' @param prevalence_ratio_tolerance Numeric length-two vector giving lower and
#'   upper acceptable ratios of anchored to nowcast mean probability. It must
#'   contain `1`, for example `c(0.9, 1.1)`.
#' @return A one-row data frame containing forward and nowcast prevalence,
#'   their ratio, the selected anchor weight, anchored prevalence and its ratio,
#'   and whether the tolerance was met.
#' @export
#'
#' @examples
#' forward <- matrix(c(0.9, 0.7, 0.5, 0.3), ncol = 1)
#' nowcast <- matrix(c(0.2, 0.2, 0.1, 0.1), ncol = 1)
#' hee_endpoint_anchor_weight(forward, nowcast, c(0.9, 1.1))
hee_endpoint_anchor_weight <- function(
    forward_occupancy,
    endpoint_nowcast,
    prevalence_ratio_tolerance = c(0.9, 1.1)) {
  forward <- as.matrix(forward_occupancy)
  endpoint <- as.matrix(endpoint_nowcast)
  if (!identical(dim(forward), dim(endpoint))) {
    stop("forward_occupancy and endpoint_nowcast must have identical dimensions.",
         call. = FALSE)
  }
  storage.mode(forward) <- "double"
  storage.mode(endpoint) <- "double"
  valid <- is.finite(forward) & is.finite(endpoint) &
    forward >= 0 & forward <= 1 & endpoint >= 0 & endpoint <= 1
  if (!all(valid)) {
    stop("forward_occupancy and endpoint_nowcast must be finite probabilities in [0, 1].",
         call. = FALSE)
  }
  tol <- suppressWarnings(as.numeric(prevalence_ratio_tolerance))
  if (length(tol) != 2L || any(!is.finite(tol)) || tol[1L] <= 0 ||
      tol[1L] > 1 || tol[2L] < 1) {
    stop("prevalence_ratio_tolerance must be finite, positive and contain 1.",
         call. = FALSE)
  }
  lower <- min(tol)
  upper <- max(tol)
  forward_mean <- mean(forward)
  endpoint_mean <- mean(endpoint)
  if (!is.finite(endpoint_mean) || endpoint_mean <= 0) {
    stop("endpoint_nowcast must have a positive mean probability for prevalence calibration.",
         call. = FALSE)
  }
  raw_ratio <- forward_mean / endpoint_mean
  weight <- 0
  if (raw_ratio > upper) {
    weight <- (raw_ratio - upper) / (raw_ratio - 1)
  } else if (raw_ratio < lower) {
    weight <- (lower - raw_ratio) / (1 - raw_ratio)
  }
  weight <- min(1, max(0, weight))
  anchored_mean <- (1 - weight) * forward_mean + weight * endpoint_mean
  anchored_ratio <- anchored_mean / endpoint_mean
  data.frame(
    forward_prevalence = forward_mean,
    endpoint_prevalence = endpoint_mean,
    forward_to_endpoint_ratio = raw_ratio,
    tolerance_lower = lower,
    tolerance_upper = upper,
    selected_anchor_weight = weight,
    anchored_prevalence = anchored_mean,
    anchored_to_endpoint_ratio = anchored_ratio,
    tolerance_met = anchored_ratio >= lower - 1e-12 &&
      anchored_ratio <= upper + 1e-12,
    stringsAsFactors = FALSE
  )
}
