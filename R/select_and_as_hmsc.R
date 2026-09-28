#' Select representative historical variables
#'
#' Selects a small representative set of historical variables for HMSC.
#'
#' @param hist An `hmsc_history_indices` object.
#' @param strategy One of "minimal", "balanced", or "full".
#' @param component One of "XData", "TrData", or "random".
#' @param drop_zero_variance Logical; drop constant/all-NA selected variables.
#' @return A data.frame of selected variables.
#' @export
select_history_variables <- function(hist,
                                     strategy = c("minimal", "balanced", "full"),
                                     component = c("XData", "TrData", "random"),
                                     drop_zero_variance = FALSE) {
  if (!inherits(hist, "hmsc_history_indices")) stop("hist must be an hmsc_history_indices object.", call. = FALSE)
  strategy <- match.arg(strategy)
  component <- match.arg(component)

  dat <- switch(component,
                XData = hist$XData_history,
                TrData = hist$TrData_history,
                random = hist$random_history)
  if (is.null(dat) || ncol(as.data.frame(dat)) == 0) return(NULL)

  available <- names(dat)
  if (drop_zero_variance) {
    available <- intersect(.nonzero_variance_names(dat), available)
  }
  wanted <- switch(paste(component, strategy, sep = ":"),
    "XData:minimal" = c(
      .first_available(available, c("historical_occupancy_balance", "colonization_age", "historical_turnover_predictor", "historical_gain_fraction", "lineage_retention_index")),
      .first_available(available, c("insitu_speciation_prop", "insitu_exsitu_balance", "insitu_diversification_strength")),
      .first_available(available, c("prob_weighted_source_diversity", "dispersal_source_diversity", "dominant_source_prop"))
    ),
    "XData:balanced" = c("colonization_age", "lineage_retention_index", "historical_loss_fraction",
                         "historical_turnover_predictor", "insitu_speciation_prop",
                         "prob_weighted_source_diversity"),
    "XData:full" = available,
    "TrData:minimal" = c("trait_conservatism", "trait_dispersal_coupling"),
    "TrData:balanced" = c("trait_conservatism", "trait_stasis_time", "trait_dispersal_coupling", "trait_insitu_coupling"),
    "TrData:full" = available,
    "random:minimal" = c("dominant_source_region"),
    "random:balanced" = c("dominant_source_region", "evoregion"),
    "random:full" = available
  )
  wanted <- unique(na.omit(wanted))
  wanted <- intersect(wanted, available)
  if (length(wanted) == 0) return(NULL)
  dat[, wanted, drop = FALSE]
}

.first_available <- function(available, priority) {
  hit <- intersect(priority, available)
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

#' Build XData for HMSC
#'
#' Combines environmental variables with selected site-level historical indices.
#'
#' @param hist An `hmsc_history_indices` object.
#' @param env Optional site-level environmental data.
#' @param strategy Variable selection strategy.
#' @param drop_zero_variance Logical; drop constant/all-NA historical variables
#'   before combining with `env`.
#' @return A site-level data.frame.
#' @export
as_hmsc_xdata <- function(hist, env = NULL, strategy = c("minimal", "balanced", "full"),
                          drop_zero_variance = TRUE) {
  strategy <- match.arg(strategy)
  Xhist <- select_history_variables(hist, strategy = strategy, component = "XData",
                                    drop_zero_variance = drop_zero_variance)
  sites <- rownames(hist$XData_history)
  if (!is.null(env)) {
    env <- .align_rows(env, sites, "env")
    if (!is.null(Xhist)) return(cbind(env, Xhist))
    return(env)
  }
  Xhist
}

#' Build TrData for HMSC
#'
#' Combines trait data with selected species-level historical trait indices.
#'
#' @param hist An `hmsc_history_indices` object.
#' @param traits Optional species trait data.
#' @param strategy Variable selection strategy.
#' @return A species-level data.frame or NULL.
#' @export
as_hmsc_trdata <- function(hist, traits = NULL, strategy = c("minimal", "balanced", "full")) {
  strategy <- match.arg(strategy)
  Trhist <- select_history_variables(hist, strategy = strategy, component = "TrData")
  species <- colnames(hist$project$comm)
  traits <- .align_traits(traits, species, "traits")
  if (!is.null(traits) && !is.null(Trhist)) {
    Trhist <- Trhist[rownames(traits), , drop = FALSE]
    # Avoid duplicated trait columns when a variable is present both in ordinary
    # traits and in trait-history summaries. The ordinary trait table is kept as
    # the primary source; the duplicated history column is dropped.
    Trhist <- Trhist[, setdiff(colnames(Trhist), colnames(traits)), drop = FALSE]
    if (ncol(Trhist) == 0) return(traits)
    return(cbind(traits, Trhist))
  }
  if (!is.null(traits)) return(traits)
  Trhist
}

#' Build historical random-effect table
#'
#' Returns selected categorical historical variables that may be used in HMSC
#' `studyDesign` and `ranLevels`.
#'
#' @param hist An `hmsc_history_indices` object.
#' @param strategy Variable selection strategy.
#' @return A site-level data.frame or NULL.
#' @export
as_hmsc_random <- function(hist, strategy = c("minimal", "balanced", "full")) {
  strategy <- match.arg(strategy)
  select_history_variables(hist, strategy = strategy, component = "random")
}

.nonzero_variance_names <- function(x) {
  if (is.null(x) || ncol(as.data.frame(x)) == 0) return(character())
  x <- as.data.frame(x)
  names(x)[vapply(x, function(z) {
    z <- z[!is.na(z)]
    if (length(z) == 0) return(FALSE)
    if (is.numeric(z) || is.integer(z)) {
      s <- stats::sd(z)
      return(!is.na(s) && s > 0)
    }
    length(unique(z)) > 1
  }, logical(1))]
}
