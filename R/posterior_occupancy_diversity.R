#' Summarize incidence diversity from binary occupancy posterior draws
#'
#' Calculates diversity only after each posterior draw or particle supplies a
#' binary cell-by-lineage occupancy state. This distinction matters because
#' the diversity of a mean occupancy field is generally not the posterior mean
#' of diversity: `D(E[Z]) != E[D(Z)]`. For every draw, lineages present in a
#' cell have equal incidence weight, so the Shannon and Gini--Simpson outputs
#' are incidence-based lineage diversity, not abundance diversity.
#'
#' This function intentionally refuses continuous probability matrices. Use
#' [hee_occupancy_weighted_diversity()] only for a clearly labelled mean-field
#' composition diagnostic; use this function for maps described as posterior
#' expected lineage richness or posterior incidence diversity.
#'
#' @param occupancy_draws A non-empty list of binary (`0`/`1`)
#'   cell-by-lineage matrices with identical dimensions and aligned dimnames.
#' @param draw_weights Optional non-negative particle weights. Uniform weights
#'   are used when omitted.
#' @param probabilities Numeric posterior summary probabilities in `[0, 1]`.
#'
#' @return A data frame with posterior means and weighted quantiles of expected
#'   sampled-lineage richness, incidence Shannon entropy, and incidence
#'   Gini--Simpson diversity for each cell.
#' @export
#'
#' @examples
#' draws <- list(
#'   d1 = rbind(c(1, 0), c(1, 1)),
#'   d2 = rbind(c(1, 1), c(0, 1))
#' )
#' hee_posterior_occupancy_diversity(draws)
hee_posterior_occupancy_diversity <- function(
    occupancy_draws,
    draw_weights = NULL,
    probabilities = c(0.025, 0.5, 0.975)) {
  if (!is.list(occupancy_draws) || !length(occupancy_draws)) {
    stop("occupancy_draws must be a non-empty list of binary matrices.",
         call. = FALSE)
  }
  probabilities <- as.numeric(probabilities)
  if (!length(probabilities) || any(!is.finite(probabilities)) ||
      any(probabilities < 0 | probabilities > 1)) {
    stop("probabilities must be finite values in [0, 1].", call. = FALSE)
  }
  draws <- lapply(seq_along(occupancy_draws), function(i) {
    x <- occupancy_draws[[i]]
    if (is.data.frame(x)) x <- as.matrix(x)
    if (!is.matrix(x) || !is.numeric(x) || !nrow(x) || !ncol(x) ||
        any(!is.finite(x)) || any(!(x %in% c(0, 1)))) {
      stop("Every occupancy draw must be a non-empty finite binary cell-by-lineage matrix; draw ",
           i, " is invalid.", call. = FALSE)
    }
    storage.mode(x) <- "double"
    x
  })
  template <- draws[[1L]]
  aligned <- vapply(draws, function(x) {
    identical(dim(x), dim(template)) &&
      (is.null(rownames(template)) || identical(rownames(x), rownames(template))) &&
      (is.null(colnames(template)) || identical(colnames(x), colnames(template)))
  }, logical(1))
  if (!all(aligned)) {
    stop("All occupancy draws must have identical dimensions and aligned dimnames.",
         call. = FALSE)
  }
  n_draws <- length(draws)
  if (is.null(draw_weights)) {
    draw_weights <- rep(1 / n_draws, n_draws)
  } else {
    draw_weights <- as.numeric(draw_weights)
    if (length(draw_weights) != n_draws || any(!is.finite(draw_weights)) ||
        any(draw_weights < 0) || sum(draw_weights) <= 0) {
      stop("draw_weights must be non-negative, finite, match occupancy_draws, and sum above zero.",
           call. = FALSE)
    }
    draw_weights <- draw_weights / sum(draw_weights)
  }

  n_cell <- nrow(template)
  richness <- shannon <- simpson <- matrix(0, nrow = n_cell, ncol = n_draws)
  for (d in seq_len(n_draws)) {
    r <- rowSums(draws[[d]])
    richness[, d] <- r
    has_lineage <- r > 0
    shannon[has_lineage, d] <- log(r[has_lineage])
    simpson[has_lineage, d] <- 1 - 1 / r[has_lineage]
  }
  weighted_quantile <- function(x, w, probs) {
    ord <- order(x)
    x <- x[ord]; w <- w[ord]
    cw <- cumsum(w) / sum(w)
    vapply(probs, function(p) x[[which(cw >= p)[[1L]]]], numeric(1))
  }
  summarize <- function(x, prefix) {
    q <- t(apply(x, 1L, weighted_quantile, w = draw_weights,
                 probs = probabilities))
    colnames(q) <- paste0(prefix, "_q", formatC(probabilities * 1000,
                                                  width = 3, flag = "0",
                                                  format = "f", digits = 0))
    cbind(
      setNames(data.frame(as.numeric(x %*% draw_weights)), paste0(prefix, "_mean")),
      as.data.frame(q)
    )
  }
  out <- cbind(
    summarize(richness, "expected_sampled_lineage_richness"),
    summarize(shannon, "posterior_incidence_shannon_entropy"),
    summarize(simpson, "posterior_incidence_gini_simpson")
  )
  out$n_posterior_draws <- n_draws
  out$posterior_effective_sample_size <- 1 / sum(draw_weights^2)
  out$scientific_boundary <-
    "binary_occupancy_draw_incidence_diversity_not_abundance_diversity_or_total_historical_flora"
  if (!is.null(rownames(template))) out$cell_id <- rownames(template)
  out
}
