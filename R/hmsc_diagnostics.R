#' Extract compact HMSC posterior diagnostics
#'
#' @param hM A fitted `Hmsc` model with posterior samples.
#' @param beta_limit Maximum number of Beta/Gamma parameters to include in
#'   MCMC diagnostic summaries.
#' @return A list containing posterior summaries, model fit, coda diagnostics,
#'   and variance partitioning when available.
#' @export
hee_hmsc_diagnostics <- function(hM, beta_limit = 200) {
  .require_pkg("Hmsc", "HMSC posterior diagnostics")
  .require_pkg("coda", "MCMC diagnostics")
  post_est <- function(name) {
    tryCatch(Hmsc::getPostEstimate(hM, name), error = function(e) NULL)
  }
  beta <- post_est("Beta")
  gamma <- post_est("Gamma")
  omega <- post_est("Omega")
  rho_draws <- .hee_hmsc_draw_table(hM, "rho")
  rho <- if (nrow(rho_draws) > 0) {
    do.call(rbind, lapply(split(rho_draws$value, rho_draws$parameter), function(z) {
      data.frame(mean = mean(z, na.rm = TRUE),
                 lwr = as.numeric(stats::quantile(z, 0.025, na.rm = TRUE)),
                 upr = as.numeric(stats::quantile(z, 0.975, na.rm = TRUE)))
    }))
  } else data.frame()
  if (nrow(rho) > 0) {
    rho <- data.frame(parameter = rownames(rho),
                      mean = rho$mean,
                      lwr = rho$lwr,
                      upr = rho$upr,
                      row.names = NULL)
  }
  pred <- tryCatch(Hmsc::computePredictedValues(hM), error = function(e) NULL)
  fit <- if (!is.null(pred)) tryCatch(Hmsc::evaluateModelFit(hM = hM, predY = pred), error = function(e) NULL) else NULL
  vp <- tryCatch(Hmsc::computeVariancePartitioning(hM), error = function(e) NULL)
  beta_draws <- .hee_hmsc_draw_table(hM, "Beta", limit = beta_limit)
  gamma_draws <- .hee_hmsc_draw_table(hM, "Gamma", limit = beta_limit)
  diag_draws <- rbind(beta_draws, gamma_draws, rho_draws)
  mcmc <- .hee_coda_diagnostics(diag_draws)
  list(beta = beta, gamma = gamma, omega = omega, rho = rho,
       predicted = pred, model_fit = fit, variance_partitioning = vp,
       mcmc_draws = diag_draws, ess = mcmc$ess, rhat = mcmc$rhat)
}

.hee_hmsc_draw_table <- function(hM, par, limit = Inf) {
  if (is.null(hM$postList) || length(hM$postList) == 0) return(data.frame())
  rows <- list()
  k <- 1L
  for (ch in seq_along(hM$postList)) {
    chain <- hM$postList[[ch]]
    for (it in seq_along(chain)) {
      obj <- chain[[it]][[par]]
      if (is.null(obj)) next
      vals <- as.numeric(obj)
      nms <- names(vals)
      if (is.null(nms)) nms <- paste0(par, "_", seq_along(vals))
      keep <- seq_len(min(length(vals), limit))
      rows[[k]] <- data.frame(chain = ch, iteration = it,
                              parameter = paste0(par, ":", nms[keep]),
                              value = vals[keep], stringsAsFactors = FALSE)
      k <- k + 1L
    }
  }
  if (length(rows) == 0) data.frame() else do.call(rbind, rows)
}

.hee_coda_diagnostics <- function(draws) {
  if (nrow(draws) == 0) return(list(ess = data.frame(), rhat = data.frame()))
  ess_rows <- list()
  rhat_rows <- list()
  for (p in unique(draws$parameter)) {
    d <- draws[draws$parameter == p, , drop = FALSE]
    chains <- split(d$value, d$chain)
    chains <- chains[vapply(chains, length, integer(1)) > 1]
    if (length(chains) == 0) next
    ml <- coda::mcmc.list(lapply(chains, coda::mcmc))
    ess_rows[[p]] <- data.frame(parameter = p,
                                ESS = as.numeric(coda::effectiveSize(ml)[1]))
    rh <- tryCatch(coda::gelman.diag(ml, autoburnin = FALSE)$psrf[1, 1],
                   error = function(e) NA_real_)
    rhat_rows[[p]] <- data.frame(parameter = p, Rhat = rh)
  }
  list(ess = if (length(ess_rows) == 0) data.frame() else do.call(rbind, ess_rows),
       rhat = if (length(rhat_rows) == 0) data.frame() else do.call(rbind, rhat_rows))
}
