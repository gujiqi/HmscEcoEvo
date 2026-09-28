#' Derive lineage-specific environmental-response reference thresholds
#'
#' HMSC response coefficients without their prevalence intercept describe a
#' response shape, but `g^{-1}(X beta)` has an arbitrary 0.5 baseline and must
#' not be summed as occurrence-based richness.  This helper derives a
#' lineage-specific reference on the **linear environmental-response** scale
#' from modern occupied survey cells.  The reference can then be reconstructed
#' along a dated tree just like any other continuous response parameter and
#' used by `hee_environmental_filtering_relative_support()`.
#'
#' The result is deliberately a calibration of a relative environmental
#' response, not an ancestral HMSC intercept, an occurrence probability, or a
#' replacement for independent endpoint validation.
#'
#' @param Y Site-by-lineage occurrence matrix with non-negative values; values
#'   greater than zero are treated as observed occurrences.
#' @param eta_tip Site-by-lineage matrix of modern `X beta` values calculated
#'   with the same locked environmental recipe as the palaeo projection.  It
#'   must have the same dimensions and lineage names as `Y`.
#' @param presence_quantile Quantile of the observed-presence response values
#'   used as the reference threshold.  A value of `0.1` defines an inclusive
#'   environmental support boundary while avoiding a single extreme record.
#' @param method Reference-calibration rule. `"presence_quantile"` uses the
#'   declared quantile among observed presences. `"prevalence_matched"` uses
#'   the quantile of all training-cell responses that has the same expected
#'   prevalence as the supplied observation column. The latter is useful for
#'   an explicitly labelled internal modern calibration when a beta-only
#'   environmental response must be compared with an occurrence-derived map.
#'   It is not independent endpoint validation and treats unrecorded cells in
#'   the same way as the supplied `Y` matrix.
#' @return A data frame with `lineage`, `eta_reference`, observed-presence
#'   counts, and transparent threshold metadata.
#' @export
#'
#' @examples
#' Y <- matrix(c(1, 0, 1, 0), ncol = 2,
#'             dimnames = list(NULL, c("a", "b")))
#' eta <- matrix(c(1, -1, 2, 0), ncol = 2,
#'               dimnames = list(NULL, c("a", "b")))
#' hee_environmental_filtering_reference_threshold(Y, eta)
hee_environmental_filtering_reference_threshold <- function(
    Y,
    eta_tip,
    presence_quantile = 0.1,
    method = c("presence_quantile", "prevalence_matched")) {
  y <- as.matrix(Y)
  eta <- as.matrix(eta_tip)
  storage.mode(y) <- "double"
  storage.mode(eta) <- "double"
  if (!identical(dim(y), dim(eta))) {
    stop("Y and eta_tip must have identical site-by-lineage dimensions.",
         call. = FALSE)
  }
  if (is.null(colnames(y)) || is.null(colnames(eta)) ||
      !identical(colnames(y), colnames(eta))) {
    stop("Y and eta_tip must have identical named lineage columns.",
         call. = FALSE)
  }
  presence_quantile <- suppressWarnings(as.numeric(presence_quantile)[1L])
  if (!is.finite(presence_quantile) || presence_quantile < 0 ||
      presence_quantile > 1) {
    stop("presence_quantile must be finite and within [0, 1].", call. = FALSE)
  }
  method <- match.arg(method)
  rows <- lapply(seq_len(ncol(y)), function(j) {
    occupied <- is.finite(y[, j]) & y[, j] > 0 & is.finite(eta[, j])
    available <- is.finite(y[, j]) & is.finite(eta[, j])
    prevalence <- if (any(available)) mean(y[available, j] > 0) else NA_real_
    value <- switch(
      method,
      presence_quantile = if (any(occupied)) {
        as.numeric(stats::quantile(eta[occupied, j],
                                   probs = presence_quantile,
                                   names = FALSE, type = 8))
      } else NA_real_,
      prevalence_matched = if (is.finite(prevalence) && prevalence > 0) {
        as.numeric(stats::quantile(eta[available, j],
                                   probs = 1 - prevalence,
                                   names = FALSE, type = 8))
      } else NA_real_
    )
    data.frame(
      lineage = colnames(y)[[j]],
      eta_reference = value,
      n_observed_presence = sum(occupied),
      training_prevalence = prevalence,
      presence_quantile = presence_quantile,
      reference_method = method,
      reference_source = if (method == "presence_quantile") {
        "modern_observed_presence_X_beta_quantile"
      } else "modern_training_prevalence_matched_X_beta_quantile",
      interpretation = paste(
        "Relative environmental-response reference; not an ancestral HMSC",
        "intercept or an independently validated occurrence-probability threshold."
      ),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Convert an environmental linear predictor to calibrated relative support
#'
#' This function centres `eta = X beta` on a lineage-specific modern-response
#' reference and maps the contrast to a bounded **relative support index**:
#' 
#' \deqn{S^*_{l,c,t}=\operatorname{logit}^{-1}
#' \left\{(\eta_{l,c,t}-\theta_{l,t})/\sigma_l\right\}.}
#'
#' `S*` is appropriate for thresholded environmental-support maps and for an
#' explicitly declared establishment covariate.  It is not a historical
#' occurrence probability: accessibility, local persistence, observation, and
#' endpoint information remain separate model components.
#'
#' @param eta Cell-by-lineage matrix of environmental linear predictors.
#' @param reference_eta Named numeric vector aligned to `colnames(eta)`, or a
#'   one-row-per-lineage data frame containing `lineage` and `eta_reference`.
#' @param scale Positive scalar or named lineage vector setting the response
#'   contrast scale.  It is a declared calibration width, not a new ecological
#'   process.
#' @return A list with `relative_support`, `binary_environmental_support`,
#'   aligned `eta_reference`, and `scale`.
#' @export
#'
#' @examples
#' eta <- matrix(c(-1, 1, 0, 2), nrow = 2,
#'               dimnames = list(NULL, c("a", "b")))
#' hee_environmental_filtering_relative_support(
#'   eta, c(a = 0, b = 1), scale = 1
#' )
hee_environmental_filtering_relative_support <- function(
    eta,
    reference_eta,
    scale = 1) {
  x <- as.matrix(eta)
  storage.mode(x) <- "double"
  if (is.null(colnames(x)) || any(!nzchar(colnames(x))) ||
      anyDuplicated(colnames(x))) {
    stop("eta must have unique non-blank lineage column names.", call. = FALSE)
  }
  ref <- if (is.data.frame(reference_eta)) {
    .require_cols(reference_eta, c("lineage", "eta_reference"),
                  "reference_eta")
    if (anyDuplicated(reference_eta$lineage)) {
      stop("reference_eta must be unique by lineage.", call. = FALSE)
    }
    stats::setNames(as.numeric(reference_eta$eta_reference),
                    as.character(reference_eta$lineage))
  } else {
    reference_eta
  }
  if (is.null(names(ref))) {
    stop("reference_eta must be named by lineage.", call. = FALSE)
  }
  ref <- suppressWarnings(as.numeric(ref[colnames(x)]))
  if (any(!is.finite(ref))) {
    stop("reference_eta is missing or non-finite for one or more eta lineages.",
         call. = FALSE)
  }
  if (length(scale) == 1L) {
    width <- rep(suppressWarnings(as.numeric(scale)[1L]), ncol(x))
  } else {
    if (is.null(names(scale))) {
      stop("A vector scale must be named by lineage.", call. = FALSE)
    }
    width <- suppressWarnings(as.numeric(scale[colnames(x)]))
  }
  if (any(!is.finite(width) | width <= 0)) {
    stop("scale must be finite and strictly positive for every lineage.",
         call. = FALSE)
  }
  contrast <- sweep(x, 2L, ref, "-")
  contrast <- sweep(contrast, 2L, width, "/")
  support <- stats::plogis(contrast)
  support[!is.finite(x)] <- NA_real_
  list(
    relative_support = support,
    binary_environmental_support = ifelse(is.finite(support) & support >= 0.5,
                                          1L, 0L),
    eta_reference = stats::setNames(ref, colnames(x)),
    scale = stats::setNames(width, colnames(x)),
    interpretation = paste(
      "Reference-calibrated relative environmental support, not historical",
      "occupancy or occurrence probability."
    )
  )
}

#' Standardise an intercept-free ancestral response over available habitat
#'
#' Direct reconstruction of HMSC `Beta` coefficients deliberately excludes the
#' modern species intercept because it can encode prevalence, sampling and
#' range size. Its raw linear predictor is therefore not an absolute
#' historical occurrence probability. This helper converts its *shape* into a
#' robust relative environmental gradient within the currently available
#' habitat at one time slice:
#'
#' \deqn{\widetilde{\eta}_{lct}=\operatorname{clip}_{[-b,b]}
#' \left\{\frac{\eta_{lct}-\operatorname{median}_{u\in H_t}(\eta_{lut})}
#' {\operatorname{IQR}_{u\in H_t}(\eta_{lut})/1.349}\right\}.}
#'
#' Use this gradient in conditional colonisation and persistence scenarios.
#' Keep a prevalence-calibrated reference support layer separate for reporting
#' and no-analog diagnosis; do not turn that reference into an extra survival
#' penalty.
#'
#' @param eta Numeric cell-by-lineage intercept-free response matrix.
#' @param habitat_state Numeric or logical vector with one entry per cell;
#'   positive/`TRUE` entries define available habitat.
#' @param clip Positive finite maximum absolute standardised value.
#' @param minimum_scale Positive fallback threshold for nearly constant
#'   response vectors.
#' @return A list containing `process_eta` and a per-lineage `audit` table.
#' @export
#'
#' @examples
#' eta <- cbind(a = c(-1, 0, 1), b = c(4, 4, 4))
#' hee_environmental_filtering_relative_gradient(eta, c(1, 1, 0))$process_eta
hee_environmental_filtering_relative_gradient <- function(
    eta,
    habitat_state,
    clip = 4,
    minimum_scale = 1e-8) {
  if (is.data.frame(eta)) eta <- as.matrix(eta)
  if (!is.matrix(eta) || !is.numeric(eta) || !nrow(eta) || !ncol(eta)) {
    stop("eta must be a non-empty numeric cell-by-lineage matrix.", call. = FALSE)
  }
  if (length(habitat_state) != nrow(eta)) {
    stop("habitat_state must have one entry per eta row.", call. = FALSE)
  }
  if (length(clip) != 1L || !is.finite(clip) || clip <= 0 ||
      length(minimum_scale) != 1L || !is.finite(minimum_scale) ||
      minimum_scale <= 0) {
    stop("clip and minimum_scale must be positive finite numbers.", call. = FALSE)
  }
  active <- !is.na(habitat_state) & as.numeric(habitat_state) > 0
  if (!any(active)) stop("habitat_state contains no available cells.", call. = FALSE)
  out <- eta
  lineage <- if (is.null(colnames(eta))) paste0("lineage_", seq_len(ncol(eta))) else colnames(eta)
  audit <- data.frame(
    lineage = lineage,
    eta_median_available_habitat = NA_real_,
    eta_robust_scale_available_habitat = NA_real_,
    n_finite_available_habitat = 0L,
    stringsAsFactors = FALSE
  )
  for (j in seq_len(ncol(eta))) {
    value <- eta[active, j]
    value <- value[is.finite(value)]
    if (!length(value)) {
      out[, j] <- NA_real_
      next
    }
    centre <- stats::median(value)
    spread <- stats::IQR(value) / 1.349
    if (!is.finite(spread) || spread < minimum_scale) spread <- stats::sd(value)
    if (!is.finite(spread) || spread < minimum_scale) spread <- 1
    out[, j] <- pmax(-clip, pmin(clip, (eta[, j] - centre) / spread))
    out[!is.finite(out[, j]), j] <- NA_real_
    audit$eta_median_available_habitat[[j]] <- centre
    audit$eta_robust_scale_available_habitat[[j]] <- spread
    audit$n_finite_available_habitat[[j]] <- length(value)
  }
  colnames(out) <- colnames(eta)
  rownames(out) <- rownames(eta)
  list(process_eta = out, audit = audit)
}
