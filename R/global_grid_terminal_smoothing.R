#' Convert a probabilistic terminal map into a soft occupancy likelihood
#'
#' A modern HMSC nowcast or an independent terminal range model can be used as
#' a *soft observation* of a binary occupancy state.  This helper returns the
#' likelihood of state zero and state one after optional power tempering.  It
#' does not alter an environmental response coefficient, and it must not be
#' described as an ancestral niche intercept.
#'
#' @param support Numeric probability matrix in `[0, 1]`, with rows as cells
#'   and columns as lineages.
#' @param power Positive likelihood power. `1` retains the supplied soft
#'   observation; values below one temper a non-independent or less reliable
#'   terminal source.
#' @param epsilon Positive clipping constant that prevents exact zero terminal
#'   likelihoods from annihilating all historical paths.
#' @return A list containing `state0_likelihood`, `state1_likelihood`, and the
#'   validated `support`.
#' @export
#'
#' @examples
#' hee_terminal_soft_likelihood(matrix(c(0.1, 0.9), ncol = 1))
hee_terminal_soft_likelihood <- function(support,
                                         power = 1,
                                         epsilon = 1e-8) {
  p <- as.matrix(support)
  storage.mode(p) <- "double"
  if (any(!is.finite(p)) || any(p < 0 | p > 1)) {
    stop("support must be finite and within [0, 1].", call. = FALSE)
  }
  power <- suppressWarnings(as.numeric(power)[1L])
  epsilon <- suppressWarnings(as.numeric(epsilon)[1L])
  if (!is.finite(power) || power <= 0 || !is.finite(epsilon) ||
      epsilon <= 0 || epsilon >= 0.5) {
    stop("power must be positive and epsilon must be in (0, 0.5).",
         call. = FALSE)
  }
  p <- pmin(pmax(p, epsilon), 1 - epsilon)
  list(
    state0_likelihood = (1 - p)^power,
    state1_likelihood = p^power,
    support = p,
    power = power,
    epsilon = epsilon
  )
}

#' Backward CTMC message for a binary occupancy transition
#'
#' Given a likelihood message at the end of an interval and the interval's
#' binary CTMC transition probabilities, calculate the likelihood message at
#' its start.  `transition_p01` and `transition_p11` may be cell-by-lineage
#' matrices, so this function remains compatible with lineage-specific
#' movement, colonisation, and persistence scenarios.
#'
#' @param transition_p01 Probability of zero-to-one transition over the
#'   interval.
#' @param transition_p11 Probability of one-to-one transition over the
#'   interval.
#' @param next_state0_likelihood Likelihood of the later state being zero.
#' @param next_state1_likelihood Likelihood of the later state being one.
#' @return A list with start-of-interval likelihood matrices for state zero and
#'   state one. Each cell-lineage pair is normalized by its larger likelihood
#'   to avoid underflow without changing posterior odds.
#' @export
#'
#' @examples
#' hee_ctmc_backward_message(matrix(0.2, 1, 1), matrix(0.8, 1, 1),
#'   matrix(0.9, 1, 1), matrix(0.1, 1, 1))
hee_ctmc_backward_message <- function(transition_p01,
                                      transition_p11,
                                      next_state0_likelihood,
                                      next_state1_likelihood) {
  p01 <- as.matrix(transition_p01)
  p11 <- as.matrix(transition_p11)
  b0_next <- as.matrix(next_state0_likelihood)
  b1_next <- as.matrix(next_state1_likelihood)
  if (!identical(dim(p01), dim(p11)) || !identical(dim(p01), dim(b0_next)) ||
      !identical(dim(p01), dim(b1_next))) {
    stop("All transition and likelihood matrices must have identical dimensions.",
         call. = FALSE)
  }
  for (x in list(p01, p11, b0_next, b1_next)) {
    if (any(!is.finite(x))) {
      stop("Transition and likelihood matrices must be finite.", call. = FALSE)
    }
  }
  if (any(p01 < 0 | p01 > 1) || any(p11 < 0 | p11 > 1) ||
      any(b0_next < 0) || any(b1_next < 0)) {
    stop("Transitions must be in [0, 1] and likelihoods must be non-negative.",
         call. = FALSE)
  }
  b0 <- (1 - p01) * b0_next + p01 * b1_next
  b1 <- (1 - p11) * b0_next + p11 * b1_next
  scale <- pmax(b0, b1, .Machine$double.xmin)
  list(
    state0_likelihood = b0 / scale,
    state1_likelihood = b1 / scale
  )
}

#' Compose two binary CTMC transition intervals
#'
#' This utility composes cell-by-lineage transition probabilities from two
#' consecutive intervals. It permits an analysis runner to retain one compact
#' `P01`/`P11` pair per palaeoenvironment interval while integrating many
#' internal numerical steps.
#'
#' @param first_p01,first_p11 Zero-to-one and one-to-one probabilities for the
#'   first interval.
#' @param second_p01,second_p11 Corresponding probabilities for the immediately
#'   following interval.
#' @return A list with composed `transition_p01` and `transition_p11` matrices.
#' @export
#'
#' @examples
#' hee_ctmc_compose_transition(
#'   matrix(0.2, 1, 1), matrix(0.8, 1, 1),
#'   matrix(0.1, 1, 1), matrix(0.9, 1, 1)
#' )
hee_ctmc_compose_transition <- function(first_p01,
                                        first_p11,
                                        second_p01,
                                        second_p11) {
  p01_a <- as.matrix(first_p01)
  p11_a <- as.matrix(first_p11)
  p01_b <- as.matrix(second_p01)
  p11_b <- as.matrix(second_p11)
  if (!identical(dim(p01_a), dim(p11_a)) ||
      !identical(dim(p01_a), dim(p01_b)) ||
      !identical(dim(p01_a), dim(p11_b))) {
    stop("All transition matrices must have identical dimensions.", call. = FALSE)
  }
  for (x in list(p01_a, p11_a, p01_b, p11_b)) {
    if (any(!is.finite(x)) || any(x < 0 | x > 1)) {
      stop("All transition matrices must be finite and within [0, 1].",
           call. = FALSE)
    }
  }
  out_p01 <- (1 - p01_a) * p01_b + p01_a * p11_b
  out_p11 <- (1 - p11_a) * p01_b + p11_a * p11_b
  list(
    transition_p01 = pmin(pmax(out_p01, 0), 1),
    transition_p11 = pmin(pmax(out_p11, 0), 1)
  )
}

#' Apply a binary occupancy likelihood message to a forward probability
#'
#' @param occupancy Forward occupancy probability matrix.
#' @param state0_likelihood Backward likelihood for state zero.
#' @param state1_likelihood Backward likelihood for state one.
#' @return A cell-by-lineage smoothed probability matrix.
#' @export
#'
#' @examples
#' hee_ctmc_apply_message(matrix(0.4, 1, 1), matrix(0.8, 1, 1),
#'   matrix(0.2, 1, 1))
hee_ctmc_apply_message <- function(occupancy,
                                   state0_likelihood,
                                   state1_likelihood) {
  q <- as.matrix(occupancy)
  b0 <- as.matrix(state0_likelihood)
  b1 <- as.matrix(state1_likelihood)
  if (!identical(dim(q), dim(b0)) || !identical(dim(q), dim(b1))) {
    stop("occupancy and likelihood matrices must have identical dimensions.",
         call. = FALSE)
  }
  if (any(!is.finite(q)) || any(q < 0 | q > 1) || any(!is.finite(b0)) ||
      any(!is.finite(b1)) || any(b0 < 0) || any(b1 < 0)) {
    stop("occupancy must be in [0, 1] and likelihoods must be finite non-negative values.",
         call. = FALSE)
  }
  numerator <- q * b1
  denominator <- numerator + (1 - q) * b0
  out <- numerator / pmax(denominator, .Machine$double.xmin)
  out[denominator <= 0] <- q[denominator <= 0]
  dimnames(out) <- dimnames(q)
  pmin(pmax(out, 0), 1)
}
