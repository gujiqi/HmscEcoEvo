#' Extract stratified HMSC Beta posterior draws for deep-time projection
#'
#' Extracts genuine retained `Beta` samples from an `Hmsc` object's `postList`
#' without averaging them first.  The output is deliberately long-but-wide:
#' each row is one modern tip lineage in one retained posterior draw, and Beta
#' axes remain columns.  It is suitable as the input to
#' [hee_evolution_ancestral_response_direct()], where each draw is
#' reconstructed separately along a dated tree before projection through a
#' palaeoenvironmental time cube.
#'
#' The default excludes the HMSC intercept.  A modern intercept commonly
#' contains prevalence, sampling, accessibility, and detection information;
#' it should not silently become an ancestral niche intercept.  When
#' `include_intercept = FALSE`, `intercept` is set to zero and the original
#' sampled intercept is retained as `hmsc_intercept_original` for audit.
#'
#' @param model A fitted `Hmsc` object with retained posterior `postList`.
#' @param n_draws Number of retained samples to select without replacement.
#' @param axes Optional environmental Beta axes. Defaults to all non-intercept
#'   HMSC covariates.
#' @param species Optional subset of HMSC species names.
#' @param seed Integer seed used only to choose retained samples.
#' @param include_intercept Logical; include the sampled HMSC intercept in the
#'   projection intercept column. Defaults to `FALSE`.
#' @return A data frame with `lineage`, `response_draw`, `chain`, `iteration`,
#'   `response_weight`, `intercept`, `hmsc_intercept_original`, and one column
#'   for every requested Beta axis. Attributes record the available and
#'   selected posterior sample counts.
#' @export
#' @examples
#' model <- list(
#'   covNames = c("(Intercept)", "temp"),
#'   spNames = c("sp1", "sp2"),
#'   postList = list(list(
#'     list(Beta = matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2)),
#'     list(Beta = matrix(c(0.2, 0.4, 0.1, 0.3), nrow = 2))
#'   ))
#' )
#' hee_hmsc_beta_posterior_draws(model, n_draws = 1, seed = 1)
hee_hmsc_beta_posterior_draws <- function(model,
                                          n_draws = 25L,
                                          axes = NULL,
                                          species = NULL,
                                          seed = 1L,
                                          include_intercept = FALSE) {
  post_list <- model[["postList"]]
  cov_names <- as.character(model[["covNames"]])
  species_all <- as.character(model[["spNames"]])
  if (is.null(post_list) || !length(post_list)) {
    stop("model must contain retained posterior samples in postList.", call. = FALSE)
  }
  if (!length(cov_names) || !length(species_all)) {
    stop("model must contain covNames and spNames for Beta orientation.", call. = FALSE)
  }
  n_draws <- suppressWarnings(as.integer(n_draws)[1L])
  if (!is.finite(n_draws) || n_draws < 1L) {
    stop("n_draws must be a positive integer.", call. = FALSE)
  }
  intercept_axis <- cov_names[grepl("^\\(Intercept\\)$|^Intercept$", cov_names)]
  if (length(intercept_axis) > 1L) {
    stop("model has more than one intercept axis; select axes explicitly.", call. = FALSE)
  }
  if (is.null(axes)) axes <- setdiff(cov_names, intercept_axis)
  axes <- as.character(axes)
  if (!all(axes %in% cov_names)) {
    stop("axes are missing from model$covNames: ",
         paste(setdiff(axes, cov_names), collapse = ", "), call. = FALSE)
  }
  if (!is.null(species)) {
    species <- as.character(species)
    if (!all(species %in% species_all)) {
      stop("species are missing from model$spNames: ",
           paste(setdiff(species, species_all), collapse = ", "), call. = FALSE)
    }
  } else {
    species <- species_all
  }

  sample_index <- do.call(rbind, lapply(seq_along(post_list), function(chain_id) {
    chain <- post_list[[chain_id]]
    if (!is.list(chain)) return(NULL)
    data.frame(chain = chain_id, iteration = seq_along(chain),
               stringsAsFactors = FALSE)
  }))
  if (is.null(sample_index) || !nrow(sample_index)) {
    stop("model$postList does not contain retained chain samples.", call. = FALSE)
  }
  if (n_draws > nrow(sample_index)) {
    stop("n_draws exceeds the number of retained posterior samples (",
         nrow(sample_index), ").", call. = FALSE)
  }

  # Allocate samples across chains as evenly as possible, then draw within each
  # chain. This avoids using an accidental contiguous segment from one chain.
  chain_sizes <- table(sample_index$chain)
  chain_ids <- as.integer(names(chain_sizes))
  base_n <- rep(n_draws %/% length(chain_ids), length(chain_ids))
  extras <- n_draws %% length(chain_ids)
  if (extras > 0L) {
    set.seed(as.integer(seed))
    extra_chain <- sample(seq_along(chain_ids), extras)
    base_n[extra_chain] <- base_n[extra_chain] + 1L
  }
  if (any(base_n > as.integer(chain_sizes))) {
    # Fall back to a globally seeded sample if a chain has too few retained
    # states for strict equal allocation.
    set.seed(as.integer(seed))
    selected_rows <- sort(sample(seq_len(nrow(sample_index)), n_draws))
  } else {
    selected_rows <- integer()
    for (i in seq_along(chain_ids)) {
      idx <- which(sample_index$chain == chain_ids[[i]])
      if (base_n[[i]] > 0L) {
        set.seed(as.integer(seed + chain_ids[[i]]))
        selected_rows <- c(selected_rows, sample(idx, base_n[[i]]))
      }
    }
    selected_rows <- sort(selected_rows)
  }
  selected <- sample_index[selected_rows, , drop = FALSE]
  selected$response_draw <- sprintf("chain%02d_iter%04d",
                                    selected$chain, selected$iteration)
  selected$response_weight <- 1 / nrow(selected)

  species_idx <- match(species, species_all)
  axis_idx <- match(axes, cov_names)
  intercept_idx <- if (length(intercept_axis)) match(intercept_axis, cov_names) else integer()
  out <- vector("list", nrow(selected))
  for (i in seq_len(nrow(selected))) {
    beta <- post_list[[selected$chain[[i]]]][[selected$iteration[[i]]]][["Beta"]]
    if (is.null(beta) || !is.matrix(beta) ||
        nrow(beta) != length(cov_names) || ncol(beta) != length(species_all)) {
      stop("Retained Beta matrix has incompatible dimensions at chain ",
           selected$chain[[i]], ", iteration ", selected$iteration[[i]], ".",
           call. = FALSE)
    }
    values <- as.data.frame(t(beta[axis_idx, species_idx, drop = FALSE]),
                            check.names = FALSE, stringsAsFactors = FALSE)
    names(values) <- axes
    values$lineage <- species
    values$response_draw <- selected$response_draw[[i]]
    values$chain <- selected$chain[[i]]
    values$iteration <- selected$iteration[[i]]
    values$response_weight <- selected$response_weight[[i]]
    hmsc_intercept <- if (length(intercept_idx)) {
      as.numeric(beta[intercept_idx, species_idx])
    } else rep(0, length(species))
    values$hmsc_intercept_original <- hmsc_intercept
    values$intercept <- if (isTRUE(include_intercept)) hmsc_intercept else 0
    out[[i]] <- values[, c("lineage", "response_draw", "chain", "iteration",
                           "response_weight", "intercept",
                           "hmsc_intercept_original", axes), drop = FALSE]
  }
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  attr(out, "n_available_posterior_samples") <- nrow(sample_index)
  attr(out, "n_selected_posterior_samples") <- nrow(selected)
  attr(out, "sampling_design") <- "stratified_by_retained_MCMC_chain"
  out
}
