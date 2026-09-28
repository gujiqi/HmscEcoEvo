.hee_response_interval_variance <- function(start_depth, end_depth, model, parameter) {
  delta <- max(end_depth - start_depth, 0)
  if (model == "OU" && parameter > 0) {
    return(-expm1(-2 * parameter * delta) / (2 * parameter))
  }
  if (model == "EB" && abs(parameter) > 0) {
    return(exp(parameter * start_depth) * expm1(parameter * delta) / parameter)
  }
  delta
}

.hee_response_transition <- function(delta, model, parameter) {
  if (model == "OU") exp(-parameter * delta) else 1
}

.hee_response_covariance <- function(geom, model, parameter = 0) {
  edge <- geom$tree$edge
  n <- length(geom$labels)
  depth <- geom$root_age - geom$node_age
  if (model == "OU") {
    covariance <- matrix(0, n, n)
    visited <- geom$root
    for (ei in seq_len(nrow(edge))) {
      parent <- edge[ei, 1L]
      child <- edge[ei, 2L]
      a <- .hee_response_transition(depth[child] - depth[parent], model, parameter)
      covariance[child, visited] <- a * covariance[parent, visited]
      covariance[visited, child] <- covariance[child, visited]
      covariance[child, child] <- a^2 * covariance[parent, parent] +
        .hee_response_interval_variance(depth[parent], depth[child], model,
                                        parameter)
      visited <- c(visited, child)
    }
    return(covariance)
  }
  paths <- matrix(0, nrow = n, ncol = nrow(edge))
  variances <- numeric(nrow(edge))
  for (ei in seq_len(nrow(edge))) {
    parent <- edge[ei, 1L]
    child <- edge[ei, 2L]
    paths[child, ] <- paths[parent, ]
    paths[child, ei] <- 1
    variances[ei] <- .hee_response_interval_variance(
      depth[parent], depth[child], model, parameter)
  }
  tcrossprod(sweep(paths, 2L, sqrt(variances), `*`))
}

.hee_fit_response_process <- function(y, geom, model) {
  n <- geom$n_tip
  ones <- rep(1, n)
  evaluate <- function(parameter) {
    covariance <- .hee_response_covariance(geom, model, parameter)
    tip_cov <- covariance[seq_len(n), seq_len(n), drop = FALSE]
    ch <- tryCatch(chol(tip_cov), error = function(e) NULL)
    if (is.null(ch)) return(NULL)
    solve_cov <- function(z) backsolve(ch, forwardsolve(t(ch), z))
    inv_one <- solve_cov(ones)
    denominator <- sum(inv_one)
    if (!is.finite(denominator) || denominator <= 0) return(NULL)
    mu <- sum(inv_one * y) / denominator
    residual <- y - mu
    inv_residual <- solve_cov(residual)
    rss <- max(sum(residual * inv_residual), 0)
    sigma2_ml <- rss / n
    log_likelihood <- if (sigma2_ml > 0) {
      -0.5 * (n * (log(2 * pi * sigma2_ml) + 1) +
                2 * sum(log(diag(ch))))
    } else NA_real_
    list(covariance = covariance, inverse_residual = inv_residual,
         mu = mu, sigma2_ml = sigma2_ml,
         sigma2_bridge = rss / (n - 1L), log_likelihood = log_likelihood)
  }

  constant_tips <- diff(range(y)) <= 1e-10 * max(1, max(abs(y)))
  boundary_status <- "none"
  if (model == "BM" || constant_tips) {
    parameter <- 0
    if (model != "BM" && constant_tips) {
      boundary_status <- "constant_tip_values"
    }
  } else {
    root_age <- geom$root_age
    if (!is.finite(root_age) || root_age <= 0) {
      stop("A dated tree with positive root age is required for OU/EB.",
           call. = FALSE)
    }
    bounds <- if (model == "OU") c(1e-6, 20) else c(-20, -1e-6)
    objective <- function(scaled_parameter) {
      fit <- evaluate(scaled_parameter / root_age)
      if (is.null(fit) || !is.finite(fit$log_likelihood)) return(1e100)
      -fit$log_likelihood
    }
    opt <- stats::optimize(objective, interval = bounds, tol = 1e-5)
    candidates <- c(bounds[1L], opt$minimum, bounds[2L])
    scores <- vapply(candidates, objective, numeric(1))
    if (all(scores >= 1e100)) {
      stop("Could not fit the ", model, " response process to tip Beta values.",
           call. = FALSE)
    }
    parameter <- candidates[which.min(scores)] / root_age
    if (which.min(scores) == 1L ||
        abs(parameter * root_age - bounds[1L]) < 1e-4) {
      boundary_status <- "lower_parameter_bound"
    } else if (which.min(scores) == length(candidates) ||
               abs(parameter * root_age - bounds[2L]) < 1e-4) {
      boundary_status <- "upper_parameter_bound"
    }
  }
  fit <- evaluate(parameter)
  if (is.null(fit)) {
    stop("The tree has singular ", model,
         " tip covariance; check zero-length terminal branches.", call. = FALSE)
  }
  k <- if (model == "BM") 2L else 3L
  fit$model <- model
  fit$parameter <- parameter
  fit$boundary_status <- boundary_status
  fit$aicc <- if (is.finite(fit$log_likelihood) && n > k + 1L) {
    2 * k - 2 * fit$log_likelihood + 2 * k * (k + 1) / (n - k - 1)
  } else Inf
  fit
}

.hee_response_branch_bridge <- function(parent_value, child_value,
                                        parent_depth, child_depth,
                                        sample_depths, fit, stochastic) {
  out <- numeric(length(sample_depths))
  previous_value <- parent_value
  previous_depth <- parent_depth
  for (i in seq_along(sample_depths)) {
    depth <- sample_depths[i]
    a <- .hee_response_transition(depth - previous_depth,
                                  fit$model, fit$parameter)
    b <- .hee_response_transition(child_depth - depth,
                                  fit$model, fit$parameter)
    v1 <- fit$sigma2_bridge * .hee_response_interval_variance(
      previous_depth, depth, fit$model, fit$parameter)
    v2 <- fit$sigma2_bridge * .hee_response_interval_variance(
      depth, child_depth, fit$model, fit$parameter)
    base <- fit$mu + a * (previous_value - fit$mu)
    expected_end <- fit$mu + b * (base - fit$mu)
    denominator <- b^2 * v1 + v2
    gain <- if (denominator > 0) v1 * b / denominator else 0
    conditional_mean <- base + gain * (child_value - expected_end)
    conditional_variance <- if (denominator > 0) v1 * v2 / denominator else 0
    value <- if (stochastic && conditional_variance > 0) {
      stats::rnorm(1L, conditional_mean, sqrt(conditional_variance))
    } else conditional_mean
    out[i] <- value
    previous_value <- value
    previous_depth <- depth
  }
  out
}
