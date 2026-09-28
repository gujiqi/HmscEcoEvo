#' Diagnose calculated historical indices
#'
#' Checks missing values, zero variance, correlations among history variables,
#' and overlap with environmental variables.
#'
#' @param hist An `hmsc_history_indices` object.
#' @param env Optional site-level environmental data.
#' @param cor_threshold Correlation threshold for warnings.
#' @return A list of diagnostic tables and messages.
#' @export
diagnose_history_indices <- function(hist, env = NULL, cor_threshold = 0.7) {
  if (!inherits(hist, "hmsc_history_indices")) stop("hist must be an hmsc_history_indices object.", call. = FALSE)
  X <- hist$XData_history
  messages <- hist$diagnostics$messages %||% character()
  warnings <- hist$diagnostics$warnings %||% character()

  missing_summary <- data.frame(
    variable = names(X),
    n_missing = vapply(X, function(z) sum(is.na(z)), integer(1)),
    stringsAsFactors = FALSE
  )

  zero_variance <- character()
  for (nm in .numeric_cols(X)) {
    z <- X[[nm]]
    if (all(is.na(z)) || stats::sd(z, na.rm = TRUE) == 0) zero_variance <- c(zero_variance, nm)
  }
  if (length(zero_variance) > 0) warnings <- c(warnings, paste("Zero-variance historical variable(s):", paste(zero_variance, collapse = ", ")))

  high_cor <- data.frame(var1 = character(), var2 = character(), cor = numeric(), stringsAsFactors = FALSE)
  num <- .numeric_cols(X)
  if (length(num) >= 2) {
    cm <- suppressWarnings(stats::cor(X[, num, drop = FALSE], use = "pairwise.complete.obs", method = "spearman"))
    idx <- which(abs(cm) >= cor_threshold & upper.tri(cm), arr.ind = TRUE)
    if (nrow(idx) > 0) {
      high_cor <- data.frame(
        var1 = rownames(cm)[idx[, 1]],
        var2 = colnames(cm)[idx[, 2]],
        cor = cm[idx],
        stringsAsFactors = FALSE
      )
      warnings <- c(warnings, paste0(nrow(high_cor), " pair(s) of historical variables exceed |r| >= ", cor_threshold, "."))
    }
  }

  env_overlap <- data.frame(history_variable = character(), env_variable = character(), cor = numeric(), stringsAsFactors = FALSE)
  if (!is.null(env)) {
    env <- .align_rows(env, rownames(X), "env")
    hnum <- .numeric_cols(X)
    enum <- .numeric_cols(env)
    rows <- list()
    for (h in hnum) for (e in enum) {
      cc <- suppressWarnings(stats::cor(X[[h]], env[[e]], use = "pairwise.complete.obs", method = "spearman"))
      if (!is.na(cc) && abs(cc) >= cor_threshold) {
        rows[[length(rows) + 1]] <- data.frame(history_variable = h, env_variable = e, cor = cc, stringsAsFactors = FALSE)
      }
    }
    env_overlap <- .safe_bind_rows(rows)
    if (is.null(env_overlap)) env_overlap <- data.frame(history_variable = character(), env_variable = character(), cor = numeric(), stringsAsFactors = FALSE)
    if (nrow(env_overlap) > 0) warnings <- c(warnings, paste0(nrow(env_overlap), " history-environment correlation(s) exceed |r| >= ", cor_threshold, "."))
  }

  list(
    messages = unique(messages),
    warnings = unique(warnings),
    missing_summary = missing_summary,
    zero_variance = zero_variance,
    high_correlations = high_cor,
    environment_overlap = env_overlap
  )
}

#' History index dictionary
#'
#' @return A data.frame describing main variables.
#' @export
history_dictionary <- function() {
  data.frame(
    variable = c(
      "historical_occupancy_balance", "colonization_age", "lineage_retention_index",
      "historical_gain_fraction", "historical_loss_fraction", "range_loss_fraction",
      "historical_turnover_predictor", "historical_gain_loss_balance",
      "insitu_speciation_prop", "prob_weighted_source_diversity",
      "dominant_source_region", "trait_conservatism", "trait_dispersal_coupling"
    ),
    group = c(
      "occupancy_dynamics", "occupancy_dynamics", "occupancy_dynamics",
      "occupancy_dynamics", "occupancy_dynamics", "occupancy_dynamics",
      "occupancy_dynamics", "occupancy_dynamics",
      "speciation", "dispersal_sources",
      "random_history", "trait_history", "trait_history"
    ),
    role = c(
      "XData", "XData", "XData", "XData", "XData", "XData", "XData",
      "XData", "XData", "XData", "random", "TrData", "TrData"
    ),
    description = c(
      "Composite predictor combining historical gain fraction, lineage retention, and loss fraction.",
      "Mean entry time of present species into the site region.",
      "Probability-weighted retention of lineages in the site region.",
      "Regional share of historical gain weight across all regions.",
      "Regional share of historical loss weight across all regions.",
      "Local fraction of gain/loss weight represented by losses.",
      "Regional share of total gain plus loss event weight.",
      "Local balance of gain and loss weight, ranging from loss-dominated to gain-dominated.",
      "Proportion of relevant speciation events occurring in the site region.",
      "Shannon diversity of historical source regions, weighted by probabilities.",
      "Dominant historical source region of the site assemblage.",
      "Species-level trait conservatism supplied by user.",
      "Coupling between dispersal trait and historical source-region breadth."
    ),
    stringsAsFactors = FALSE
  )
}
