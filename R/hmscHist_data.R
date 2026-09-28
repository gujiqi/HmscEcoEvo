#' Create a historical-data project object
#'
#' Store the basic community, environment, phylogeny, coordinate, and trait data
#' used by the historical-predictor workflow. This function does not calculate
#' historical indices and does not run HMSC. `hmscHist_data()` is retained as a
#' backward-compatible alias.
#'
#' @param comm Site by species matrix/data.frame with 0/1 values.
#' @param site_region Data frame with columns `site` and `region`.
#' @param env Optional site-level environmental data with site rownames.
#' @param phy Optional phylogenetic tree object.
#' @param coords Optional site coordinates with site rownames.
#' @param traits Optional species trait data with species rownames or a `species` column.
#' @param x A `HmscEcoEvo_data` object.
#' @param object A `HmscEcoEvo_data` object.
#' @param ... Additional arguments, currently ignored.
#' @return An object of class `c("HmscEcoEvo_data", "hmscHist_data")`.
#' @export
HmscEcoEvo_data <- function(comm, site_region, env = NULL, phy = NULL,
                            coords = NULL, traits = NULL) {
  comm <- .as_matrix01(comm, "comm")

  if (!is.data.frame(site_region)) {
    stop("site_region must be a data.frame with columns 'site' and 'region'.", call. = FALSE)
  }
  .require_cols(site_region, c("site", "region"), "site_region")
  site_region$site <- as.character(site_region$site)
  site_region$region <- as.character(site_region$region)

  missing_sites <- setdiff(rownames(comm), site_region$site)
  if (length(missing_sites) > 0) {
    stop("site_region is missing site(s): ", paste(missing_sites, collapse = ", "), call. = FALSE)
  }
  site_region <- site_region[match(rownames(comm), site_region$site), , drop = FALSE]
  rownames(site_region) <- site_region$site

  env <- .align_rows(env, rownames(comm), "env")
  coords <- .align_rows(coords, rownames(comm), "coords")
  traits <- .align_traits(traits, colnames(comm), "traits")

  out <- list(
    comm = comm,
    site_region = site_region,
    env = env,
    phy = phy,
    coords = coords,
    traits = traits,
    history_tables = list(),
    history_diagnostics = list(messages = character(), warnings = character()),
    available_indices = character()
  )
  class(out) <- c("HmscEcoEvo_data", "hmscHist_data")
  out
}

#' @rdname HmscEcoEvo_data
#' @export
print.HmscEcoEvo_data <- function(x, ...) {
  cat("HmscEcoEvo_data object\n")
  cat("  sites:   ", nrow(x$comm), "\n", sep = "")
  cat("  species: ", ncol(x$comm), "\n", sep = "")
  cat("  regions: ", length(unique(x$site_region$region)), "\n", sep = "")
  cat("  env:     ", if (is.null(x$env)) "no" else paste0(ncol(x$env), " variables"), "\n", sep = "")
  cat("  traits:  ", if (is.null(x$traits)) "no" else paste0(ncol(x$traits), " variables"), "\n", sep = "")
  invisible(x)
}

#' @rdname HmscEcoEvo_data
#' @export
summary.HmscEcoEvo_data <- function(object, ...) {
  print(object)
  if (length(object$history_tables) > 0) {
    cat("\nHistory tables:\n")
    for (nm in names(object$history_tables)) {
      tab <- object$history_tables[[nm]]
      cat("  ", nm, ": ", nrow(tab), " rows\n", sep = "")
    }
  }
  if (length(object$available_indices) > 0) {
    cat("\nAvailable index groups/indices:\n")
    cat("  ", paste(object$available_indices, collapse = ", "), "\n", sep = "")
  }
  invisible(object)
}

#' @rdname HmscEcoEvo_data
#' @export
hmscHist_data <- HmscEcoEvo_data

#' @rdname HmscEcoEvo_data
#' @export
print.hmscHist_data <- print.HmscEcoEvo_data

#' @rdname HmscEcoEvo_data
#' @export
summary.hmscHist_data <- summary.HmscEcoEvo_data
