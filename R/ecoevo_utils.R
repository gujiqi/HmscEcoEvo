# Eco-evolutionary helper functions ----------------------------------------

.require_pkg <- function(pkg, reason = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- paste0("Package '", pkg, "' is required")
    if (!is.null(reason)) msg <- paste0(msg, " for ", reason)
    stop(msg, ". Install it or pass a precomputed result.", call. = FALSE)
  }
  invisible(TRUE)
}

.soft_pkg <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

.as_df <- function(x, name = "x") {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) return(x)
  if (is.matrix(x)) return(as.data.frame(x))
  stop(name, " must be a data.frame or matrix.", call. = FALSE)
}

.ensure_named_matrix <- function(x, name = "x") {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x)) stop(name, " must be a matrix or data.frame.", call. = FALSE)
  storage.mode(x) <- "numeric"
  if (is.null(rownames(x)) || any(rownames(x) == "")) {
    stop(name, " must have row names.", call. = FALSE)
  }
  if (is.null(colnames(x)) || any(colnames(x) == "")) {
    stop(name, " must have column names.", call. = FALSE)
  }
  x
}

.coerce_beta_matrix <- function(beta, species = NULL, axes = NULL,
                                name = "beta") {
  B <- .ensure_named_matrix(beta, name)
  if (is.null(B)) return(NULL)

  if (!is.null(species)) species <- as.character(species)
  if (!is.null(axes)) axes <- as.character(axes)

  row_species <- !is.null(species) && all(species %in% rownames(B))
  col_species <- !is.null(species) && all(species %in% colnames(B))
  row_axes <- !is.null(axes) && all(axes %in% rownames(B))
  col_axes <- !is.null(axes) && all(axes %in% colnames(B))

  if (row_species || col_axes) {
    out <- B
  } else if (col_species || row_axes) {
    out <- t(B)
  } else {
    # HMSC commonly stores Beta as covariates x species. If row names look like
    # covariates and column names like species, transpose; otherwise keep rows.
    if (nrow(B) < ncol(B)) out <- t(B) else out <- B
  }

  if (!is.null(species)) {
    missing <- setdiff(species, rownames(out))
    if (length(missing) > 0) {
      stop(name, " is missing species: ", paste(missing, collapse = ", "),
           call. = FALSE)
    }
    out <- out[species, , drop = FALSE]
  }
  if (!is.null(axes)) {
    missing <- setdiff(axes, colnames(out))
    if (length(missing) > 0) {
      stop(name, " is missing axis/axes: ", paste(missing, collapse = ", "),
           call. = FALSE)
    }
    out <- out[, axes, drop = FALSE]
  }
  out
}

.infer_beta_species <- function(beta, species = NULL, traits = NULL,
                                phylo = NULL, Y = NULL, omega = NULL) {
  if (!is.null(species)) return(as.character(species))
  candidates <- list()
  if (!is.null(traits)) {
    tr <- as.data.frame(traits)
    if ("species" %in% names(tr)) candidates[[length(candidates) + 1L]] <- as.character(tr$species)
    if (!is.null(rownames(tr))) candidates[[length(candidates) + 1L]] <- rownames(tr)
  }
  if (!is.null(phylo) && inherits(phylo, "phylo")) {
    candidates[[length(candidates) + 1L]] <- phylo$tip.label
  }
  if (!is.null(Y) && !is.null(colnames(Y))) {
    candidates[[length(candidates) + 1L]] <- colnames(Y)
  }
  if (!is.null(omega)) {
    if (!is.null(rownames(omega))) candidates[[length(candidates) + 1L]] <- rownames(omega)
    if (!is.null(colnames(omega))) candidates[[length(candidates) + 1L]] <- colnames(omega)
  }
  for (cand in candidates) {
    cand <- cand[!is.na(cand) & cand != ""]
    if (length(cand) > 0 && all(cand %in% rownames(beta))) return(cand)
    if (length(cand) > 0 && all(cand %in% colnames(beta))) return(cand)
  }
  row_has_intercept <- any(.is_intercept_axis(rownames(beta)))
  col_has_intercept <- any(.is_intercept_axis(colnames(beta)))
  if (row_has_intercept && !col_has_intercept) return(colnames(beta))
  if (col_has_intercept && !row_has_intercept) return(rownames(beta))
  NULL
}

.is_intercept_axis <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y %in% c("intercept", "(intercept)", "constant", "const", "beta0", "alpha")
}

.non_intercept_axes <- function(x) {
  x[!.is_intercept_axis(x)]
}

.coerce_gamma_matrix <- function(gamma, traits = NULL, axes = NULL,
                                 name = "gamma") {
  G <- .ensure_named_matrix(gamma, name)
  if (is.null(G)) return(NULL)
  trait_names <- if (is.null(traits)) NULL else {
    if (is.data.frame(traits) || is.matrix(traits)) colnames(traits) else as.character(traits)
  }
  if (!is.null(trait_names) && all(trait_names %in% rownames(G))) {
    out <- G
  } else if (!is.null(trait_names) && all(trait_names %in% colnames(G))) {
    out <- t(G)
  } else if (!is.null(axes) && all(axes %in% colnames(G))) {
    out <- G
  } else if (!is.null(axes) && all(axes %in% rownames(G))) {
    out <- t(G)
  } else {
    out <- G
  }
  if (!is.null(trait_names)) {
    missing <- setdiff(trait_names, rownames(out))
    if (length(missing) > 0) {
      stop(name, " is missing traits: ", paste(missing, collapse = ", "),
           call. = FALSE)
    }
    out <- out[trait_names, , drop = FALSE]
  }
  if (!is.null(axes)) {
    missing <- setdiff(axes, colnames(out))
    if (length(missing) > 0) {
      stop(name, " is missing axis/axes: ", paste(missing, collapse = ", "),
           call. = FALSE)
    }
    out <- out[, axes, drop = FALSE]
  }
  out
}

.coerce_beta_draws <- function(draws, species = NULL, axes = NULL,
                               beta = NULL, name = "beta_draws") {
  if (is.null(draws)) return(NULL)
  if (!is.array(draws) || length(dim(draws)) != 3) {
    stop(name, " must be a 3D array with dimensions draw x species x axis, ",
         "or draw x axis x species.", call. = FALSE)
  }
  dn <- dimnames(draws)
  if (is.null(species) && !is.null(beta)) species <- rownames(beta)
  if (is.null(axes) && !is.null(beta)) axes <- colnames(beta)

  d2 <- if (!is.null(dn[[2]])) dn[[2]] else character()
  d3 <- if (!is.null(dn[[3]])) dn[[3]] else character()

  if (!is.null(species) && all(species %in% d2)) {
    out <- draws
  } else if (!is.null(species) && all(species %in% d3)) {
    out <- aperm(draws, c(1, 3, 2))
  } else if (!is.null(axes) && all(axes %in% d2)) {
    out <- aperm(draws, c(1, 3, 2))
  } else {
    out <- draws
  }

  if (!is.null(species)) {
    if (is.null(dimnames(out)[[2]])) stop(name, " needs species dimnames.", call. = FALSE)
    missing <- setdiff(species, dimnames(out)[[2]])
    if (length(missing) > 0) stop(name, " is missing species: ", paste(missing, collapse = ", "), call. = FALSE)
    out <- out[, species, , drop = FALSE]
  }
  if (!is.null(axes)) {
    if (is.null(dimnames(out)[[3]])) stop(name, " needs axis dimnames.", call. = FALSE)
    missing <- setdiff(axes, dimnames(out)[[3]])
    if (length(missing) > 0) stop(name, " is missing axis/axes: ", paste(missing, collapse = ", "), call. = FALSE)
    out <- out[, , axes, drop = FALSE]
  }
  out
}

.coerce_gamma_draws <- function(draws, traits = NULL, axes = NULL,
                                gamma = NULL, name = "gamma_draws") {
  if (is.null(draws)) return(NULL)
  if (!is.array(draws) || length(dim(draws)) != 3) {
    stop(name, " must be a 3D array with dimensions draw x trait x axis, ",
         "or draw x axis x trait.", call. = FALSE)
  }
  trait_names <- if (!is.null(traits)) {
    if (is.data.frame(traits) || is.matrix(traits)) colnames(traits) else as.character(traits)
  } else if (!is.null(gamma)) rownames(gamma) else NULL
  if (is.null(axes) && !is.null(gamma)) axes <- colnames(gamma)

  dn <- dimnames(draws)
  d2 <- if (!is.null(dn[[2]])) dn[[2]] else character()
  d3 <- if (!is.null(dn[[3]])) dn[[3]] else character()
  if (!is.null(trait_names) && all(trait_names %in% d2)) {
    out <- draws
  } else if (!is.null(trait_names) && all(trait_names %in% d3)) {
    out <- aperm(draws, c(1, 3, 2))
  } else if (!is.null(axes) && all(axes %in% d2)) {
    out <- aperm(draws, c(1, 3, 2))
  } else {
    out <- draws
  }
  if (!is.null(trait_names)) {
    if (is.null(dimnames(out)[[2]])) stop(name, " needs trait dimnames.", call. = FALSE)
    missing <- setdiff(trait_names, dimnames(out)[[2]])
    if (length(missing) > 0) stop(name, " is missing traits: ", paste(missing, collapse = ", "), call. = FALSE)
    out <- out[, trait_names, , drop = FALSE]
  }
  if (!is.null(axes)) {
    if (is.null(dimnames(out)[[3]])) stop(name, " needs axis dimnames.", call. = FALSE)
    missing <- setdiff(axes, dimnames(out)[[3]])
    if (length(missing) > 0) stop(name, " is missing axis/axes: ", paste(missing, collapse = ", "), call. = FALSE)
    out <- out[, , axes, drop = FALSE]
  }
  out
}

.coerce_square_draws <- function(draws, labels, name = "draws") {
  if (is.null(draws)) return(NULL)
  if (!is.array(draws) || length(dim(draws)) != 3) {
    stop(name, " must be a 3D array with dimensions draw x row x column.",
         call. = FALSE)
  }
  dn <- dimnames(draws)
  if (is.null(dn[[2]]) || is.null(dn[[3]])) {
    stop(name, " needs row and column dimnames.", call. = FALSE)
  }
  missing <- setdiff(labels, intersect(dn[[2]], dn[[3]]))
  if (length(missing) > 0) {
    stop(name, " is missing labels: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  draws[, labels, labels, drop = FALSE]
}

.posterior_summary_array <- function(draws, value_name = "value") {
  if (is.null(draws)) return(NULL)
  dims <- dim(draws)
  if (length(dims) != 3) stop("draws must be a 3D array.", call. = FALSE)
  dn <- dimnames(draws)
  rows <- dn[[2]]
  cols <- dn[[3]]
  if (is.null(rows)) rows <- paste0("row", seq_len(dims[2]))
  if (is.null(cols)) cols <- paste0("col", seq_len(dims[3]))
  out <- vector("list", dims[2] * dims[3])
  idx <- 1L
  for (i in seq_len(dims[2])) {
    for (j in seq_len(dims[3])) {
      z <- draws[, i, j]
      p_pos <- mean(z > 0, na.rm = TRUE)
      p_neg <- mean(z < 0, na.rm = TRUE)
      out[[idx]] <- data.frame(
        row = rows[i],
        col = cols[j],
        mean = mean(z, na.rm = TRUE),
        median = stats::median(z, na.rm = TRUE),
        lwr = as.numeric(stats::quantile(z, 0.025, na.rm = TRUE)),
        upr = as.numeric(stats::quantile(z, 0.975, na.rm = TRUE)),
        sd = stats::sd(z, na.rm = TRUE),
        ci_width = diff(as.numeric(stats::quantile(z, c(0.025, 0.975), na.rm = TRUE))),
        p_pos = p_pos,
        p_neg = p_neg,
        support = max(p_pos, p_neg),
        sign = sign(mean(z, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  out <- do.call(rbind, out)
  names(out)[names(out) == "row"] <- "species"
  names(out)[names(out) == "col"] <- "axis"
  names(out)[names(out) == "mean"] <- value_name
  rownames(out) <- NULL
  out
}

.matrix_to_long <- function(x, value_name = "value",
                            row_name = "row", col_name = "col") {
  if (is.null(x)) return(data.frame())
  x <- as.matrix(x)
  data.frame(
    setNames(list(rep(rownames(x), times = ncol(x))), row_name),
    setNames(list(rep(colnames(x), each = nrow(x))), col_name),
    setNames(list(as.vector(x)), value_name),
    stringsAsFactors = FALSE
  )
}

.safe_lm_r2 <- function(y, data) {
  y <- suppressWarnings(as.numeric(y))
  data <- .as_df(data)
  if (is.null(data) || ncol(data) == 0 || length(y) != nrow(data)) return(NA_real_)
  ok <- !is.na(y)
  if (sum(ok) <= 2) return(NA_real_)
  d <- data[ok, , drop = FALSE]
  y <- y[ok]
  keep <- vapply(d, function(z) length(unique(z[!is.na(z)])) > 1, logical(1))
  d <- d[, keep, drop = FALSE]
  if (ncol(d) == 0) return(NA_real_)
  fit <- tryCatch(stats::lm(y ~ ., data = d), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  s <- suppressWarnings(summary(fit))
  as.numeric(s$r.squared)
}

.safe_lm_coefs <- function(y, data) {
  y <- suppressWarnings(as.numeric(y))
  data <- .as_df(data)
  if (is.null(data) || ncol(data) == 0 || length(y) != nrow(data)) return(NULL)
  ok <- !is.na(y)
  if (sum(ok) <= 2) return(NULL)
  d <- data[ok, , drop = FALSE]
  y <- y[ok]
  keep <- vapply(d, function(z) length(unique(z[!is.na(z)])) > 1, logical(1))
  d <- d[, keep, drop = FALSE]
  if (ncol(d) == 0) return(NULL)
  fit <- tryCatch(stats::lm(y ~ ., data = d), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  co <- suppressWarnings(stats::coef(fit))
  co <- co[names(co) != "(Intercept)"]
  data.frame(term = names(co), estimate = as.numeric(co), stringsAsFactors = FALSE)
}

.align_species_frame <- function(x, species, name = "data") {
  if (is.null(x)) return(NULL)
  x <- .as_df(x, name)
  if ("species" %in% names(x)) {
    rownames(x) <- as.character(x$species)
    x$species <- NULL
  }
  if (is.null(rownames(x))) stop(name, " must have species row names or a species column.", call. = FALSE)
  missing <- setdiff(species, rownames(x))
  if (length(missing) > 0) {
    stop(name, " is missing species: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  x[species, , drop = FALSE]
}

.align_clade <- function(clade, species) {
  if (is.null(clade)) return(NULL)
  if (is.data.frame(clade)) {
    .require_cols(clade, c("species", "clade"), "clade")
    z <- setNames(as.character(clade$clade), as.character(clade$species))
  } else {
    z <- as.character(clade)
    if (is.null(names(z))) {
      if (length(z) != length(species)) stop("clade must be named or have one value per species.", call. = FALSE)
      names(z) <- species
    }
  }
  missing <- setdiff(species, names(z))
  if (length(missing) > 0) stop("clade is missing species: ", paste(missing, collapse = ", "), call. = FALSE)
  z[species]
}

.align_phylo <- function(phylo, species, name = "phylo") {
  if (is.null(phylo)) return(NULL)
  .require_pkg("ape", "phylogeny alignment")
  if (!inherits(phylo, "phylo")) stop(name, " must be an ape phylo object.", call. = FALSE)
  missing <- setdiff(species, phylo$tip.label)
  if (length(missing) > 0) {
    stop(name, " is missing species: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  extra <- setdiff(phylo$tip.label, species)
  if (length(extra) > 0) phylo <- ape::drop.tip(phylo, extra)
  phylo$tip.label <- phylo$tip.label
  phylo
}

.phylo_dist <- function(phylo) {
  .require_pkg("ape", "phylogenetic distances")
  ape::cophenetic.phylo(phylo)
}

.phylo_moran_i <- function(z, phylo) {
  d <- .phylo_dist(phylo)
  z <- z[rownames(d)]
  ok <- !is.na(z)
  if (sum(ok) <= 2) return(NA_real_)
  d <- d[ok, ok, drop = FALSE]
  z <- z[ok]
  w <- 1 / d
  diag(w) <- 0
  w[!is.finite(w)] <- 0
  zc <- z - mean(z)
  denom <- sum(zc^2)
  if (denom == 0 || sum(w) == 0) return(NA_real_)
  length(z) / sum(w) * sum(w * tcrossprod(zc)) / denom
}

.local_moran <- function(z, phylo) {
  d <- .phylo_dist(phylo)
  z <- z[rownames(d)]
  w <- 1 / d
  diag(w) <- 0
  w[!is.finite(w)] <- 0
  zc <- z - mean(z, na.rm = TRUE)
  m2 <- mean(zc^2, na.rm = TRUE)
  if (is.na(m2) || m2 == 0) return(rep(NA_real_, length(z)))
  out <- zc * (w %*% zc) / m2
  out <- as.numeric(out)
  names(out) <- names(z)
  out
}

.correlogram_table <- function(B, phylo, bins = 5, permutations = 0) {
  d <- .phylo_dist(phylo)
  d <- d[rownames(B), rownames(B), drop = FALSE]
  pd <- d[upper.tri(d)]
  cuts <- unique(stats::quantile(pd, probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE))
  if (length(cuts) <= 2) cuts <- pretty(pd, n = bins)
  out <- list()
  k <- 1L
  for (axis in colnames(B)) {
    z <- B[, axis]
    zc <- z - mean(z, na.rm = TRUE)
    denom <- sum(zc^2, na.rm = TRUE)
    zz <- as.matrix(stats::dist(z, upper = TRUE))
    bd <- zz[upper.tri(zz)]
    bin <- cut(pd, breaks = cuts, include.lowest = TRUE)
    lev <- levels(bin)
    tab <- lapply(lev, function(lb) {
      keep <- !is.na(bin) & bin == lb
      w <- matrix(0, nrow(d), ncol(d), dimnames = dimnames(d))
      w[upper.tri(w)] <- as.numeric(keep)
      w <- w + t(w)
      moran <- if (denom == 0 || sum(w) == 0) NA_real_ else {
        length(z) / sum(w) * sum(w * tcrossprod(zc), na.rm = TRUE) / denom
      }
      p_value <- NA_real_
      if (is.finite(moran) && permutations > 0) {
        perm_vals <- replicate(permutations, {
          zp <- sample(z)
          zcp <- zp - mean(zp, na.rm = TRUE)
          denom_p <- sum(zcp^2, na.rm = TRUE)
          if (denom_p == 0 || sum(w) == 0) NA_real_ else {
            length(zp) / sum(w) * sum(w * tcrossprod(zcp), na.rm = TRUE) / denom_p
          }
        })
        perm_vals <- perm_vals[is.finite(perm_vals)]
        if (length(perm_vals) > 0) {
          p_value <- (sum(abs(perm_vals) >= abs(moran), na.rm = TRUE) + 1) /
            (length(perm_vals) + 1)
        }
      }
      data.frame(
        distance_bin = lb,
        beta_distance = mean(bd[keep], na.rm = TRUE),
        phylo_distance = mean(pd[keep], na.rm = TRUE),
        moran_i = moran,
        lwr = NA_real_,
        upr = NA_real_,
        p_value = p_value,
        signif = ifelse(is.na(p_value), "",
                  ifelse(p_value < 0.001, "***",
                  ifelse(p_value < 0.01, "**",
                  ifelse(p_value < 0.05, "*", "")))),
        stringsAsFactors = FALSE
      )
    })
    tab <- do.call(rbind, tab)
    tab$axis <- axis
    tab$similarity <- tab$moran_i
    tab$source <- "binned_phylogenetic_moran_i"
    out[[k]] <- tab
    k <- k + 1L
  }
  do.call(rbind, out)
}

.node_signal_table <- function(B, phylo) {
  .require_pkg("ape", "node-level signal")
  ntip <- length(phylo$tip.label)
  nodes <- (ntip + 1):(ntip + phylo$Nnode)
  out <- list()
  k <- 1L
  for (node in nodes) {
    tips <- tryCatch(ape::extract.clade(phylo, node)$tip.label, error = function(e) character())
    tips <- intersect(tips, rownames(B))
    if (length(tips) < 2) next
    for (axis in colnames(B)) {
      z <- B[tips, axis]
      v <- stats::var(z, na.rm = TRUE)
      out[[k]] <- data.frame(
        node = node,
        axis = axis,
        n_species = length(tips),
        mean_beta = mean(z, na.rm = TRUE),
        var_beta = v,
        node_signal = ifelse(is.na(v), NA_real_, 1 / (1 + v)),
        species = paste(tips, collapse = ";"),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  if (length(out) == 0) data.frame() else do.call(rbind, out)
}

.clade_cci <- function(B, clade) {
  if (is.null(clade)) return(data.frame())
  clade <- .align_clade(clade, rownames(B))
  global_dist <- stats::dist(B)
  gmean <- mean(as.numeric(global_dist), na.rm = TRUE)
  out <- lapply(split(names(clade), clade), function(sp) {
    if (length(sp) < 2 || is.na(gmean) || gmean == 0) {
      val <- NA_real_
      within <- NA_real_
    } else {
      within <- mean(as.numeric(stats::dist(B[sp, , drop = FALSE])), na.rm = TRUE)
      val <- 1 - within / gmean
    }
    data.frame(clade = clade[sp[1]], n_species = length(sp),
               mean_within_beta_distance = within,
               CCI = val, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

.optional_table <- function(x, required = character(), name = "table") {
  if (is.null(x)) return(NULL)
  x <- as.data.frame(x)
  if (length(required) > 0) .require_cols(x, required, name)
  x
}

.new_plot_list <- function(plots, title = NULL) {
  plots <- plots[!vapply(plots, is.null, logical(1))]
  structure(list(title = title, plots = plots),
            class = "hmsc_ecoevo_plot_list")
}

#' Print a HmscEcoEvo plot list
#'
#' @param x A plot-list object returned by one of the diagnostic plotting
#'   functions when `patchwork` is not available.
#' @param ... Additional arguments, currently ignored.
#' @return Invisibly returns `x`.
#' @export
print.hmsc_ecoevo_plot_list <- function(x, ...) {
  if (!is.null(x$title)) cat(x$title, "\n", sep = "")
  for (p in x$plots) print(p)
  invisible(x)
}

.combine_plots <- function(plots, title = NULL, ncol = 2) {
  plots <- plots[!vapply(plots, is.null, logical(1))]
  if (.soft_pkg("patchwork")) {
    out <- patchwork::wrap_plots(plots, ncol = ncol)
    if (!is.null(title)) out <- out + patchwork::plot_annotation(title = title)
    out
  } else {
    .new_plot_list(plots, title = title)
  }
}

.gg_ok <- function() .require_pkg("ggplot2", "diagnostic plotting")

.empty_metric <- function(message) {
  list(diagnostics = list(messages = message, warnings = character()))
}

.as_precomputed_list <- function(precomputed) {
  if (is.null(precomputed)) list() else precomputed
}

.support_from_matrix <- function(M, value_name = "value",
                                 row_name = "species", col_name = "axis") {
  out <- .matrix_to_long(M, value_name = value_name,
                         row_name = row_name, col_name = col_name)
  out$median <- out[[value_name]]
  out$lwr <- NA_real_
  out$upr <- NA_real_
  out$sd <- NA_real_
  out$ci_width <- NA_real_
  out$p_pos <- NA_real_
  out$p_neg <- NA_real_
  out$support <- NA_real_
  out$sign <- sign(out[[value_name]])
  out
}

.posterior_ci_vector <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  c(
    mean = mean(x, na.rm = TRUE),
    median = stats::median(x, na.rm = TRUE),
    lwr = as.numeric(stats::quantile(x, 0.025, na.rm = TRUE)),
    upr = as.numeric(stats::quantile(x, 0.975, na.rm = TRUE)),
    sd = stats::sd(x, na.rm = TRUE),
    support = max(mean(x > 0, na.rm = TRUE), mean(x < 0, na.rm = TRUE))
  )
}

.get_history_species_table <- function(history, species) {
  if (is.null(history)) return(NULL)
  if (inherits(history, "hmsc_history_indices")) {
    x <- history$TrData_history
  } else {
    x <- history
  }
  .align_species_frame(x, species, "history")
}

.get_history_site_table <- function(history) {
  if (is.null(history) || !inherits(history, "hmsc_history_indices")) return(NULL)
  .as_df(history$XData_history, "history$XData_history")
}

.merge_history_XData <- function(XData, history) {
  hist_x <- .get_history_site_table(history)
  XData <- .as_df(XData, "XData")
  if (is.null(hist_x)) return(XData)
  if (is.null(XData)) return(hist_x)
  if (!is.null(rownames(XData)) && !is.null(rownames(hist_x))) {
    missing <- setdiff(rownames(XData), rownames(hist_x))
    if (length(missing) > 0) {
      stop("history$XData_history is missing sites found in XData: ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    hist_x <- hist_x[rownames(XData), , drop = FALSE]
  } else if (nrow(XData) != nrow(hist_x)) {
    stop("XData and history$XData_history must have matching row names or the same number of rows.",
         call. = FALSE)
  }
  dup <- intersect(names(hist_x), names(XData))
  if (length(dup) > 0) {
    names(hist_x)[match(dup, names(hist_x))] <- paste0("history_", dup)
  }
  cbind(XData, hist_x)
}

.get_model_slot <- function(model, names) {
  if (is.null(model)) return(NULL)
  for (nm in names) {
    if (!is.null(model[[nm]])) return(model[[nm]])
  }
  NULL
}

.try_hmsc_post <- function(model, par_name) {
  if (is.null(model) || !inherits(model, "Hmsc") || !.soft_pkg("Hmsc")) return(NULL)
  aliases <- unique(c(par_name, if (identical(par_name, "Rho")) "rho" else character()))
  for (nm in aliases) {
    est <- tryCatch(Hmsc::getPostEstimate(model, parName = nm), error = function(e) NULL)
    if (!is.null(est)) return(est)
  }
  NULL
}

.extract_post_mean <- function(est) {
  if (is.null(est)) return(NULL)
  if (is.list(est) && !is.null(est$mean)) return(est$mean)
  if (is.matrix(est) || is.array(est)) return(est)
  NULL
}

.hmsc_post_samples <- function(model) {
  if (is.null(model) || is.null(model$postList)) return(list())
  out <- list()
  k <- 1L
  for (chain in model$postList) {
    if (is.list(chain) && length(chain) > 0) {
      for (sample in chain) {
        if (is.list(sample)) {
          out[[k]] <- sample
          k <- k + 1L
        }
      }
    }
  }
  out
}

.hmsc_matrix_dimnames <- function(model, par_name, dim_value) {
  rn <- cn <- NULL
  if (identical(par_name, "Beta")) {
    rn <- model$covNames
    cn <- model$spNames
  } else if (identical(par_name, "Gamma")) {
    rn <- model$covNames
    cn <- model$trNames
  } else if (par_name %in% c("Omega", "OmegaCor")) {
    rn <- cn <- model$spNames
  }
  if (!is.null(rn) && length(rn) != dim_value[1]) rn <- NULL
  if (!is.null(cn) && length(cn) != dim_value[2]) cn <- NULL
  list(rn = rn, cn = cn)
}

.apply_hmsc_dimnames <- function(x, model, par_name) {
  if (is.null(x) || !(is.matrix(x) || (is.array(x) && length(dim(x)) == 2))) {
    return(x)
  }
  dn <- .hmsc_matrix_dimnames(model, par_name, dim(x))
  if (is.null(rownames(x)) && !is.null(dn$rn)) rownames(x) <- dn$rn
  if (is.null(colnames(x)) && !is.null(dn$cn)) colnames(x) <- dn$cn
  x
}

.extract_hmsc_draw_array <- function(model, par_name, max_draws = Inf) {
  if (is.null(model)) return(NULL)
  if (par_name %in% c("Omega", "OmegaCor")) {
    om <- .extract_hmsc_omega_draw_array(model, cor = identical(par_name, "OmegaCor"),
                                         max_draws = max_draws)
    if (!is.null(om)) return(om)
  }
  samples <- .hmsc_post_samples(model)
  direct_name <- if (identical(par_name, "Rho")) "rho" else par_name
  if (length(samples) > 0 && direct_name %in% names(samples[[1]])) {
    collected <- list()
    for (sample in samples) {
      z <- sample[[direct_name]]
      if (is.matrix(z) || (is.array(z) && length(dim(z)) == 2)) {
        z <- as.matrix(z)
        z <- .apply_hmsc_dimnames(z, model, par_name)
        collected[[length(collected) + 1L]] <- z
      }
      if (length(collected) >= max_draws) break
    }
    if (length(collected) > 0) {
      dims <- unique(vapply(collected, function(x) paste(dim(x), collapse = "x"), character(1)))
      if (length(dims) == 1L) {
        rn <- rownames(collected[[1]])
        cn <- colnames(collected[[1]])
        out <- array(NA_real_,
          dim = c(length(collected), nrow(collected[[1]]), ncol(collected[[1]])),
          dimnames = list(NULL, rn, cn)
        )
        for (i in seq_along(collected)) out[i, , ] <- collected[[i]]
        return(out)
      }
    }
  }
  collected <- list()
  walk <- function(obj, depth = 0L, nm = NULL) {
    if (depth > 8L || length(collected) >= max_draws) return()
    if (!is.null(nm) && identical(nm, par_name) && (is.matrix(obj) || (is.array(obj) && length(dim(obj)) == 2))) {
      collected[[length(collected) + 1L]] <<- as.matrix(obj)
      return()
    }
    if (is.list(obj)) {
      nms <- names(obj)
      if (is.null(nms)) nms <- rep(NA_character_, length(obj))
      for (i in seq_along(obj)) walk(obj[[i]], depth + 1L, nms[i])
    }
  }
  if (!is.null(model$postList)) walk(model$postList, 0L)
  if (length(collected) == 0L) walk(model, 0L)
  if (length(collected) == 0L) return(NULL)
  dims <- unique(vapply(collected, function(x) paste(dim(x), collapse = "x"), character(1)))
  if (length(dims) != 1L) return(NULL)
  rn <- rownames(collected[[1]])
  cn <- colnames(collected[[1]])
  out <- array(NA_real_, dim = c(length(collected), nrow(collected[[1]]), ncol(collected[[1]])),
               dimnames = list(NULL, rn, cn))
  for (i in seq_along(collected)) out[i, , ] <- collected[[i]]
  out
}

.extract_hmsc_vector_draws <- function(model, par_name, max_draws = Inf) {
  if (is.null(model)) return(NULL)
  samples <- .hmsc_post_samples(model)
  direct_name <- if (identical(par_name, "Rho")) "rho" else par_name
  if (length(samples) == 0 || !direct_name %in% names(samples[[1]])) return(NULL)
  collected <- lapply(samples, function(sample) {
    z <- sample[[direct_name]]
    if (is.null(z) || is.matrix(z) || is.array(z) || is.list(z)) return(NULL)
    as.numeric(z)
  })
  collected <- collected[!vapply(collected, is.null, logical(1))]
  if (length(collected) == 0) return(NULL)
  len <- unique(lengths(collected))
  if (length(len) != 1L) return(NULL)
  if (is.finite(max_draws)) collected <- collected[seq_len(min(length(collected), max_draws))]
  out <- do.call(rbind, collected)
  if (ncol(out) == 1) colnames(out) <- direct_name else colnames(out) <- paste0(direct_name, seq_len(ncol(out)))
  out
}

.extract_hmsc_omega_draw_array <- function(model, r = 1, x = NULL,
                                           cor = FALSE, max_draws = Inf) {
  samples <- .hmsc_post_samples(model)
  if (length(samples) == 0 || is.null(model$ranLevels) ||
      length(model$ranLevels) < r) {
    return(NULL)
  }
  collected <- list()
  for (sample in samples) {
    if (is.null(sample$Lambda) || length(sample$Lambda) < r) next
    lambda <- sample$Lambda[[r]]
    if (is.null(lambda)) next
    om <- tryCatch({
      if (is.null(model$ranLevels[[r]]$xDim) || model$ranLevels[[r]]$xDim == 0) {
        crossprod(lambda)
      } else {
        if (is.null(x)) x <- c(1, rep(0, model$ranLevels[[r]]$xDim - 1))
        dim_l <- dim(lambda)
        crossprod(rowSums(lambda * array(rep(x, each = prod(dim_l[1:2])), dim_l),
                          dims = 2))
      }
    }, error = function(e) NULL)
    if (is.null(om)) next
    if (cor) om <- stats::cov2cor(om)
    om <- as.matrix(om)
    if (!is.null(model$spNames) && length(model$spNames) == nrow(om)) {
      dimnames(om) <- list(model$spNames, model$spNames)
    }
    collected[[length(collected) + 1L]] <- om
    if (length(collected) >= max_draws) break
  }
  if (length(collected) == 0) return(NULL)
  dims <- unique(vapply(collected, function(z) paste(dim(z), collapse = "x"), character(1)))
  if (length(dims) != 1L) return(NULL)
  out <- array(NA_real_, dim = c(length(collected), nrow(collected[[1]]), ncol(collected[[1]])),
               dimnames = list(NULL, rownames(collected[[1]]), colnames(collected[[1]])))
  for (i in seq_along(collected)) out[i, , ] <- collected[[i]]
  out
}

.axis_pairs <- function(axes, quadratic = NULL) {
  if (!is.null(quadratic)) {
    q <- as.data.frame(quadratic, stringsAsFactors = FALSE)
    .require_cols(q, c("linear", "quadratic"), "quadratic")
    if (!"axis" %in% names(q)) q$axis <- q$linear
    return(q[, c("axis", "linear", "quadratic"), drop = FALSE])
  }
  out <- list()
  k <- 1L
  for (ax in axes) {
    pats <- c(
      paste0("I\\(", gsub("([\\W])", "\\\\\\1", ax), "\\^2\\)"),
      paste0(ax, "2"),
      paste0(ax, "_sq"),
      paste0(ax, "_squared"),
      paste0(ax, "\\^2")
    )
    hit <- axes[Reduce(`|`, lapply(pats, function(p) grepl(p, axes)))]
    hit <- setdiff(hit, ax)
    if (length(hit) > 0) {
      out[[k]] <- data.frame(axis = ax, linear = ax, quadratic = hit[1],
                             stringsAsFactors = FALSE)
      k <- k + 1L
    }
  }
  if (length(out) == 0) {
    data.frame(axis = character(), linear = character(), quadratic = character())
  } else {
    do.call(rbind, out)
  }
}

.link_response <- function(eta, distr = "normal") {
  d <- tolower(as.character(distr)[1])
  if (grepl("probit|bernoulli|binomial", d)) return(stats::pnorm(eta))
  if (grepl("poisson|count", d)) return(exp(eta))
  eta
}

.standardize_population_beta <- function(population) {
  if (is.null(population)) return(NULL)
  x <- as.data.frame(population)
  .require_cols(x, c("species", "population"), "population")
  if (all(c("axis", "beta") %in% names(x))) return(x)
  value_cols <- setdiff(names(x), c("species", "population", "clade", "region"))
  value_cols <- value_cols[vapply(x[value_cols], is.numeric, logical(1))]
  out <- lapply(value_cols, function(v) {
    data.frame(species = x$species, population = x$population,
               axis = v, beta = x[[v]], stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

.pairwise_beta_phylo <- function(B, phylo = NULL) {
  beta_d <- as.matrix(stats::dist(B))
  sp <- rownames(beta_d)
  phy_d <- if (is.null(phylo)) matrix(NA_real_, length(sp), length(sp),
                                      dimnames = list(sp, sp)) else .phylo_dist(phylo)[sp, sp]
  out <- data.frame(
    species1 = sp[row(beta_d)[upper.tri(beta_d)]],
    species2 = sp[col(beta_d)[upper.tri(beta_d)]],
    beta_distance = beta_d[upper.tri(beta_d)],
    phylo_distance = phy_d[upper.tri(phy_d)],
    stringsAsFactors = FALSE
  )
  out
}
