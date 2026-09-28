#' Create a unified HMSC eco-evolutionary object
#'
#' `hmsc_ecoevo()` standardizes HMSC outputs and auxiliary data into a single
#' S3 object used by all HmscEcoEvo diagnostic functions. Internally, Beta is
#' stored as species x niche-axis, and Gamma is stored as trait x niche-axis.
#' Complex comparative-model outputs such as bayou/SURFACE shifts or OU/BM/EB
#' fits should be supplied through `precomputed`; HmscEcoEvo will not invent
#' these model-based results from insufficient data.
#'
#' @param model Optional fitted Hmsc model. If supplied, the constructor tries
#'   to extract Beta, Gamma, Rho, Omega, posterior draws, species names,
#'   covariate names, traits, phylogeny, XData, and Y when explicit arguments
#'   are missing.
#' @param beta Matrix/data.frame of Beta estimates. Either species x axis or
#'   axis x species; row/column names are required.
#' @param species Optional character vector giving the intended species order.
#'   Supplying this removes ambiguity when Beta orientation cannot be inferred.
#' @param axes Optional character vector giving the intended niche-axis order.
#' @param beta_draws Optional 3D posterior array, draw x species x axis or draw
#'   x axis x species.
#' @param gamma Matrix/data.frame of Gamma estimates. Either trait x axis or
#'   axis x trait; row/column names are required.
#' @param gamma_draws Optional 3D posterior array, draw x trait x axis or draw x
#'   axis x trait.
#' @param rho Optional HMSC rho summary, vector, data.frame, or posterior draws.
#' @param rho_draws Optional rho posterior draws.
#' @param omega Optional residual association matrix.
#' @param omega_draws Optional residual association posterior draws.
#' @param traits Optional species x trait table.
#' @param phylo Optional `ape::phylo` tree with tip labels matching species.
#' @param history Optional `hmsc_history_indices` object or species-level
#'   historical predictor table.
#' @param population Optional population-level niche data. Long form should
#'   contain `species`, `population`, `axis`, and `beta`.
#' @param timeseries Optional time-series table for dynamic eco-evolutionary
#'   diagnostics.
#' @param XData Optional site-level environmental/history predictors.
#' @param Y Optional community matrix.
#' @param coords Optional site coordinates.
#' @param distr Optional HMSC distribution label. Used for response-curve link
#'   functions.
#' @param scaling Optional list describing model scaling.
#' @param clade Optional named clade vector or table with `species`, `clade`.
#' @param regime Optional named regime vector or table with `species`,
#'   `regime`.
#' @param precomputed Optional list of precomputed advanced diagnostic tables.
#' @param metadata Optional free-form metadata list.
#' @return An object of class `hmsc_ecoevo`.
#' @export
#' @examples
#' beta <- matrix(c(0.3, -0.2, 0.1, 0.4), nrow = 2,
#'   dimnames = list(c("sp1", "sp2"), c("temp", "precip")))
#' evo <- hmsc_ecoevo(beta = beta)
#' evo
hmsc_ecoevo <- function(model = NULL,
                        beta = NULL,
                        species = NULL,
                        axes = NULL,
                        beta_draws = NULL,
                        gamma = NULL,
                        gamma_draws = NULL,
                        rho = NULL,
                        rho_draws = NULL,
                        omega = NULL,
                        omega_draws = NULL,
                        traits = NULL,
                        phylo = NULL,
                        history = NULL,
                        population = NULL,
                        timeseries = NULL,
                        XData = NULL,
                        Y = NULL,
                        coords = NULL,
                        distr = NULL,
                        scaling = NULL,
                        clade = NULL,
                        regime = NULL,
                        precomputed = list(),
                        metadata = list()) {
  if (is.null(species) && !is.null(model$spNames)) species <- model$spNames
  if (is.null(axes) && !is.null(model$covNames)) axes <- model$covNames
  if (is.null(traits) && !is.null(model$TrData)) traits <- model$TrData
  if (is.null(phylo) && !is.null(model$phyloTree)) phylo <- model$phyloTree
  if (is.null(XData) && !is.null(model$XData)) XData <- model$XData
  if (is.null(Y) && !is.null(model$Y)) Y <- model$Y

  if (is.null(beta)) {
    beta <- .extract_post_mean(.try_hmsc_post(model, "Beta"))
    beta <- .apply_hmsc_dimnames(beta, model, "Beta")
    if (is.null(beta)) beta <- .get_model_slot(model, c("Beta", "beta"))
  }
  if (is.null(gamma)) {
    gamma <- .extract_post_mean(.try_hmsc_post(model, "Gamma"))
    gamma <- .apply_hmsc_dimnames(gamma, model, "Gamma")
    if (is.null(gamma)) gamma <- .get_model_slot(model, c("Gamma", "gamma"))
  }
  if (is.null(rho)) {
    rho <- .extract_post_mean(.try_hmsc_post(model, "Rho"))
    if (is.null(rho)) rho <- .get_model_slot(model, c("Rho", "rho"))
  }
  if (is.null(omega)) {
    omega <- .extract_post_mean(.try_hmsc_post(model, "Omega"))
    omega <- .apply_hmsc_dimnames(omega, model, "Omega")
    if (is.null(omega)) omega <- .get_model_slot(model, c("Omega", "omega"))
  }

  if (is.null(beta_draws)) beta_draws <- .extract_hmsc_draw_array(model, "Beta")
  if (is.null(beta) && !is.null(beta_draws)) {
    beta <- apply(beta_draws, c(2, 3), mean, na.rm = TRUE)
  }

  species_arg <- species
  species <- NULL
  if (!is.null(beta)) {
    raw_beta <- .ensure_named_matrix(beta, "beta")
    species <- .infer_beta_species(raw_beta, species = species_arg, traits = traits,
                                   phylo = phylo, Y = Y, omega = omega)
    if (is.null(species)) {
      if (!is.null(axes) && all(axes %in% rownames(raw_beta))) {
        species <- colnames(raw_beta)
      } else if (!is.null(axes) && all(axes %in% colnames(raw_beta))) {
        species <- rownames(raw_beta)
      } else {
        warning("Beta orientation is ambiguous; assuming rows are species. ",
                "Pass `species` and/or `axes` to avoid this heuristic.",
                call. = FALSE)
        species <- rownames(raw_beta)
      }
    }
    beta <- .coerce_beta_matrix(raw_beta, species = species, axes = axes, name = "beta")
  }

  if (is.null(beta) && !is.null(beta_draws)) {
    dn <- dimnames(beta_draws)
    if (is.null(dn[[2]]) || is.null(dn[[3]])) {
      stop("beta_draws needs species and axis dimnames when beta is NULL.", call. = FALSE)
    }
    beta <- apply(beta_draws, c(2, 3), mean, na.rm = TRUE)
  }
  if (is.null(beta)) {
    stop("A Beta matrix or beta_draws array is required, or it must be extractable from model.",
         call. = FALSE)
  }

  species <- rownames(beta)
  axes <- colnames(beta)
  beta_draws <- .coerce_beta_draws(beta_draws, species = species, axes = axes,
                                   beta = beta, name = "beta_draws")

  traits <- .align_species_frame(traits, species, "traits")
  history_species <- .get_history_species_table(history, species)

  if (is.null(gamma_draws)) gamma_draws <- .extract_hmsc_draw_array(model, "Gamma")
  if (is.null(gamma) && !is.null(gamma_draws)) {
    gamma <- apply(gamma_draws, c(2, 3), mean, na.rm = TRUE)
  }
  gamma <- .coerce_gamma_matrix(gamma, traits = traits, axes = axes, name = "gamma")
  gamma_draws <- .coerce_gamma_draws(gamma_draws, traits = traits, axes = axes,
                                     gamma = gamma, name = "gamma_draws")

  if (is.null(rho_draws)) rho_draws <- .extract_hmsc_vector_draws(model, "Rho")
  if (is.null(rho) && !is.null(rho_draws)) rho <- colMeans(as.matrix(rho_draws), na.rm = TRUE)

  if (is.null(omega_draws)) omega_draws <- .extract_hmsc_draw_array(model, "Omega")
  if (is.null(omega) && !is.null(omega_draws)) {
    omega <- apply(omega_draws, c(2, 3), mean, na.rm = TRUE)
  }

  if (!is.null(omega)) {
    omega <- .ensure_named_matrix(omega, "omega")
    if (!all(species %in% rownames(omega)) || !all(species %in% colnames(omega))) {
      stop("omega must contain all Beta species as row and column names.", call. = FALSE)
    }
    omega <- omega[species, species, drop = FALSE]
  }
  omega_draws <- .coerce_square_draws(omega_draws, species, "omega_draws")

  phylo <- .align_phylo(phylo, species, "phylo")
  clade <- .align_clade(clade, species)
  regime <- .align_regime(regime, species)

  if (!is.null(Y)) {
    Y <- .ensure_named_matrix(Y, "Y")
    if (!all(species %in% colnames(Y))) {
      stop("Y must contain all Beta species as columns.", call. = FALSE)
    }
    Y <- Y[, species, drop = FALSE]
  }
  XData <- .merge_history_XData(XData, history)
  coords <- .as_df(coords, "coords")
  population <- .standardize_population_beta(population)
  timeseries <- .as_df(timeseries, "timeseries")

  out <- list(
    model = model,
    beta = beta,
    beta_draws = beta_draws,
    gamma = gamma,
    gamma_draws = gamma_draws,
    rho = rho,
    rho_draws = rho_draws,
    omega = omega,
    omega_draws = omega_draws,
    traits = traits,
    phylo = phylo,
    history = history,
    history_species = history_species,
    population = population,
    timeseries = timeseries,
    XData = XData,
    Y = Y,
    coords = coords,
    distr = distr %||% "normal",
    scaling = scaling,
    clade = clade,
    regime = regime,
    precomputed = .as_precomputed_list(precomputed),
    metadata = metadata,
    species = species,
    axes = axes
  )
  attr(out, "scaling") <- scaling
  class(out) <- "hmsc_ecoevo"
  validate_hmsc_ecoevo(out)
  out
}

#' Coerce an object to `hmsc_ecoevo`
#'
#' @param x An object. Existing `hmsc_ecoevo` objects are returned unchanged;
#'   other objects are treated as `model` and passed to [hmsc_ecoevo()].
#' @param ... Arguments passed to [hmsc_ecoevo()].
#' @return A `hmsc_ecoevo` object.
#' @export
as_hmsc_ecoevo <- function(x, ...) {
  if (inherits(x, "hmsc_ecoevo")) return(x)
  hmsc_ecoevo(model = x, ...)
}

#' Extract HMSC outputs into `hmsc_ecoevo`
#'
#' @param model A fitted Hmsc model.
#' @param ... Additional arguments passed to [hmsc_ecoevo()].
#' @return A `hmsc_ecoevo` object.
#' @export
extract_hmsc_ecoevo <- function(model, ...) {
  hmsc_ecoevo(model = model, ...)
}

#' Validate a `hmsc_ecoevo` object
#'
#' @param x Object to validate.
#' @return Invisibly returns `TRUE`.
#' @export
validate_hmsc_ecoevo <- function(x) {
  if (!inherits(x, "hmsc_ecoevo")) stop("x must be a hmsc_ecoevo object.", call. = FALSE)
  B <- .ensure_named_matrix(x$beta, "x$beta")
  if (!identical(rownames(B), x$species)) stop("x$species must match rownames(beta).", call. = FALSE)
  if (!identical(colnames(B), x$axes)) stop("x$axes must match colnames(beta).", call. = FALSE)
  if (!is.null(x$beta_draws)) {
    if (!identical(dimnames(x$beta_draws)[[2]], x$species) ||
        !identical(dimnames(x$beta_draws)[[3]], x$axes)) {
      stop("beta_draws dimnames must match beta species and axes.", call. = FALSE)
    }
  }
  if (!is.null(x$traits) && !identical(rownames(x$traits), x$species)) {
    stop("traits row names must match beta species.", call. = FALSE)
  }
  if (!is.null(x$history_species) && !identical(rownames(x$history_species), x$species)) {
    stop("species-level history row names must match beta species.", call. = FALSE)
  }
  TRUE
}

#' @export
print.hmsc_ecoevo <- function(x, ...) {
  cat("hmsc_ecoevo object\n")
  cat("  species: ", length(x$species), "\n", sep = "")
  cat("  niche axes: ", length(x$axes), "\n", sep = "")
  cat("  beta draws: ", if (is.null(x$beta_draws)) "no" else dim(x$beta_draws)[1], "\n", sep = "")
  cat("  gamma: ", if (is.null(x$gamma)) "no" else paste(dim(x$gamma), collapse = " x "), "\n", sep = "")
  cat("  phylogeny: ", if (is.null(x$phylo)) "no" else "yes", "\n", sep = "")
  cat("  traits: ", if (is.null(x$traits)) "no" else ncol(x$traits), "\n", sep = "")
  cat("  population data: ", if (is.null(x$population)) "no" else nrow(x$population), "\n", sep = "")
  invisible(x)
}

#' @export
summary.hmsc_ecoevo <- function(object, ...) {
  print(object)
  cat("\nNiche axes:\n")
  cat("  ", paste(object$axes, collapse = ", "), "\n", sep = "")
  if (!is.null(object$scaling)) {
    cat("\nScaling metadata is available in attr(x, 'scaling').\n")
  }
  if (length(object$precomputed) > 0) {
    cat("\nPrecomputed diagnostics:\n")
    cat("  ", paste(names(object$precomputed), collapse = ", "), "\n", sep = "")
  }
  invisible(object)
}

.align_regime <- function(regime, species) {
  if (is.null(regime)) return(NULL)
  if (is.data.frame(regime)) {
    .require_cols(regime, c("species", "regime"), "regime")
    z <- setNames(as.character(regime$regime), as.character(regime$species))
  } else {
    z <- as.character(regime)
    if (is.null(names(z))) {
      if (length(z) != length(species)) stop("regime must be named or have one value per species.", call. = FALSE)
      names(z) <- species
    }
  }
  missing <- setdiff(species, names(z))
  if (length(missing) > 0) stop("regime is missing species: ", paste(missing, collapse = ", "), call. = FALSE)
  z[species]
}
