#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
} else {
  normalizePath("HmscEcoEvo/inst/extdata/ecoevo_cases/generate_ecoevo_cases.R", mustWork = TRUE)
}
case_root <- dirname(script_path)
pkg_dir <- normalizePath(file.path(case_root, "..", "..", ".."), mustWork = TRUE)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(pkg_dir, quiet = TRUE)
} else {
  library(HmscEcoEvo)
}
if (!requireNamespace("ape", quietly = TRUE)) {
  stop("Package 'ape' is required to generate the time-scaled toy phylogeny.", call. = FALSE)
}

set.seed(20260605)

species <- sprintf("sp%02d", 1:15)
linear_axes <- c("temp", "precip", "forest", "soil", "pH")
quadratic_axes <- paste0(linear_axes, "_sq")
axes <- c("(Intercept)", linear_axes, quadratic_axes)
analysis_axes <- c(linear_axes, quadratic_axes)
quadratic_map <- data.frame(linear = linear_axes, quadratic = quadratic_axes,
                            stringsAsFactors = FALSE)
traits_names <- c("body_size", "SLA", "seed_mass", "wood_density", "leaf_N",
                  "root_depth", "height", "flowering_time", "dispersal_mode")
regime_levels <- c("cold_wet", "warm_dry", "mesic_forest", "open_dry",
                   "alpine", "wetland")
case_names <- sprintf("case%02d_%s", 1:7, c(
  "niche_summary",
  "phylo_signal",
  "evo_transition",
  "trait_history_mediation",
  "gamma_evolution",
  "population_evolution",
  "validation_dashboard"
))

root_age_by_case <- c(46, 58, 67, 41, 73, 52, 63)
make_time_tree <- function(case_id) {
  root_age <- root_age_by_case[case_id]
  crown_frac <- c(A = 0.34, B = 0.27, C = 0.42, D = 0.22, E = 0.31) +
    c(A = 0.015, B = -0.010, C = 0.020, D = 0.012, E = -0.018) * case_id
  crown_age <- pmin(root_age - 4, pmax(8, root_age * crown_frac))
  names(crown_age) <- names(crown_frac)
  split_frac <- c(A = 0.42, B = 0.55, C = 0.37, D = 0.62, E = 0.48) +
    c(A = -0.010, B = 0.006, C = 0.008, D = -0.005, E = 0.004) * case_id
  split_age <- pmin(crown_age - 1.2, pmax(1.8, crown_age * split_frac))
  names(split_age) <- names(split_frac)
  clade_tips <- split(species, rep(LETTERS[1:5], each = 3))
  subtrees <- vapply(names(clade_tips), function(cl) {
    sp <- clade_tips[[cl]]
    sprintf("((%s:%.4f,%s:%.4f):%.4f,%s:%.4f):%.4f",
            sp[1], split_age[cl], sp[2], split_age[cl],
            crown_age[cl] - split_age[cl], sp[3], crown_age[cl],
            root_age - crown_age[cl])
  }, character(1))
  tr <- ape::read.tree(text = paste0("(", paste(subtrees, collapse = ","), ");"))
  tr$root.time <- root_age
  tr
}

tree_root_age <- function(phylo) {
  max(ape::node.depth.edgelength(phylo), na.rm = TRUE)
}

branch_time_table <- function(phylo) {
  depth <- ape::node.depth.edgelength(phylo)
  root_age <- tree_root_age(phylo)
  data.frame(
    branch = seq_len(nrow(phylo$edge)),
    parent = phylo$edge[, 1],
    child = phylo$edge[, 2],
    branch_length = phylo$edge.length,
    time_start_Ma = root_age - depth[phylo$edge[, 1]],
    time_end_Ma = root_age - depth[phylo$edge[, 2]],
    stringsAsFactors = FALSE
  )
}

clade <- setNames(rep(LETTERS[1:5], each = 3), species)
regime <- setNames(rep(regime_levels, length.out = length(species)), species)
clade_score <- setNames(rep(seq(-1.2, 1.2, length.out = 5), each = 3), species)

scale_num <- function(x) as.numeric(scale(x))

make_traits <- function(case_id) {
  idx <- seq_along(species)
  cl <- clade_score[species]
  raw <- data.frame(
    body_size = 0.70 * cl + sin(idx / 2) + rnorm(length(idx), 0, 0.12),
    SLA = -0.55 * cl + cos(idx / 3) + rnorm(length(idx), 0, 0.10),
    seed_mass = 0.45 * cl + rep(c(-0.4, 0.2, 0.6), 5) + rnorm(length(idx), 0, 0.10),
    wood_density = 0.35 * cl + rep(c(0.5, -0.1, -0.4), 5) + rnorm(length(idx), 0, 0.09),
    leaf_N = -0.40 * cl + rep(c(0.3, -0.5, 0.4), 5) + rnorm(length(idx), 0, 0.10),
    root_depth = 0.50 * cl + rep(c(-0.5, 0.0, 0.5), 5) + rnorm(length(idx), 0, 0.10),
    height = 0.60 * cl + sin(idx / 4) + rnorm(length(idx), 0, 0.11),
    flowering_time = -0.35 * cl + cos(idx / 5) + rnorm(length(idx), 0, 0.10),
    dispersal_mode = 0.25 * cl + rep(c(-0.7, 0.1, 0.7), 5) + rnorm(length(idx), 0, 0.09)
  )
  raw <- as.data.frame(lapply(raw, scale_num))
  rownames(raw) <- species
  raw + case_id * 0.01
}

make_history <- function(case_id) {
  idx <- seq_along(species)
  cl <- clade_score[species]
  raw <- data.frame(
    paleo_temp = -0.60 * cl + sin(idx / 3) + rnorm(length(idx), 0, 0.08),
    paleo_precip = 0.50 * cl + cos(idx / 4) + rnorm(length(idx), 0, 0.08),
    forest_stability = rep(c(0.8, 0.4, -0.2, -0.6, 0.1), each = 3) +
      rnorm(length(idx), 0, 0.08),
    soil_pc1 = 0.35 * cl + rep(c(-0.5, 0.2, 0.7), 5) + rnorm(length(idx), 0, 0.08),
    pH_legacy = -0.30 * cl + rep(c(0.6, -0.1, -0.4), 5) + rnorm(length(idx), 0, 0.08),
    colonization_age = seq(50, 8, length.out = length(idx)) + rnorm(length(idx), 0, 1.2),
    isolation_index = rep(c(0.2, 0.5, 0.9, 0.3, 0.7), each = 3) +
      rnorm(length(idx), 0, 0.06)
  )
  for (nm in names(raw)) raw[[nm]] <- scale_num(raw[[nm]])
  rownames(raw) <- species
  raw + case_id * 0.005
}

make_beta <- function(traits, history, case_id) {
  trait_coef <- matrix(c(
     0.42, -0.18,  0.25,  0.20, -0.12,
    -0.30,  0.36, -0.10,  0.05,  0.18,
     0.12,  0.10,  0.26,  0.18,  0.08,
    -0.18, -0.22,  0.34,  0.44, -0.24,
     0.10, -0.08, -0.18,  0.06,  0.42,
    -0.28,  0.20,  0.02, -0.32,  0.12,
     0.22, -0.06,  0.16,  0.10, -0.18,
    -0.12,  0.14, -0.20,  0.06,  0.24,
     0.10,  0.18,  0.22, -0.12, -0.06
  ), nrow = length(traits_names), byrow = TRUE,
  dimnames = list(traits_names, linear_axes))
  history_coef <- matrix(c(
    -0.22,  0.10,  0.08, -0.02,  0.04,
     0.12,  0.28, -0.10,  0.05, -0.02,
     0.04, -0.12,  0.30,  0.14, -0.06,
     0.02,  0.06,  0.10,  0.28, -0.12,
    -0.08,  0.04, -0.04, -0.10,  0.30,
     0.05, -0.04,  0.06, -0.02,  0.05,
     0.10, -0.08,  0.04,  0.06, -0.06
  ), nrow = ncol(history), byrow = TRUE,
  dimnames = list(colnames(history), linear_axes))
  clade_shift <- matrix(c(
     0.45, -0.20,  0.30, -0.10,  0.10,
     0.15,  0.35,  0.20,  0.18, -0.25,
    -0.25,  0.12,  0.42,  0.08,  0.16,
    -0.10, -0.35, -0.15,  0.40,  0.30,
     0.34, -0.18, -0.30, -0.22,  0.48
  ), nrow = 5, byrow = TRUE, dimnames = list(LETTERS[1:5], linear_axes))
  B_lin <- as.matrix(traits) %*% trait_coef +
    as.matrix(history) %*% history_coef +
    clade_shift[clade[species], , drop = FALSE]
  B_lin <- B_lin / max(abs(B_lin), na.rm = TRUE) * 1.05
  hidden_species <- rep(c(-1.15, 0.05, 1.15), 5)
  hidden_axis <- c(temp = 0.42, precip = -0.34, forest = 0.38, soil = -0.36, pH = 0.32)
  B_lin <- B_lin + outer(hidden_species, hidden_axis)
  B_lin <- B_lin + matrix(rnorm(length(B_lin), 0, 0.07 + case_id * 0.002),
                          nrow = nrow(B_lin), dimnames = dimnames(B_lin))
  B_quad <- -0.20 - 0.04 * abs(B_lin) +
    matrix(rnorm(length(B_lin), 0, 0.025), nrow = nrow(B_lin),
           dimnames = list(species, quadratic_axes))
  colnames(B_quad) <- quadratic_axes
  intercept <- -0.15 + 0.12 * traits$body_size - 0.10 * history$isolation_index +
    rnorm(length(species), 0, 0.04)
  B <- cbind("(Intercept)" = intercept, B_lin, B_quad)
  B[, axes, drop = FALSE]
}

make_gamma <- function(case_id) {
  G_lin <- matrix(c(
     0.63, -0.41,  0.24,  0.39, -0.18,
    -0.48,  0.56, -0.33, -0.05,  0.21,
     0.21,  0.14,  0.37,  0.28,  0.12,
    -0.22, -0.31,  0.41,  0.62, -0.27,
     0.12, -0.08, -0.19,  0.07,  0.58,
    -0.36,  0.29,  0.03, -0.42,  0.11,
     0.28, -0.11,  0.22,  0.18, -0.16,
    -0.18,  0.21, -0.24,  0.10,  0.35,
     0.16,  0.27,  0.31, -0.18, -0.09
  ), nrow = length(traits_names), byrow = TRUE,
  dimnames = list(traits_names, linear_axes))
  G_quad <- -0.18 * abs(G_lin) + matrix(rnorm(length(G_lin), 0, 0.025),
                                        nrow = nrow(G_lin),
                                        dimnames = list(traits_names, quadratic_axes))
  colnames(G_quad) <- quadratic_axes
  intercept <- matrix(rnorm(length(traits_names), 0, 0.04), ncol = 1,
                      dimnames = list(traits_names, "(Intercept)"))
  G <- cbind(intercept, G_lin, G_quad)
  G + matrix(rnorm(length(G), 0, 0.012 + case_id * 0.001),
             nrow = nrow(G), dimnames = dimnames(G))
}

make_draws <- function(M, n_draws = 120, sd = 0.05) {
  out <- array(NA_real_, dim = c(n_draws, nrow(M), ncol(M)),
               dimnames = list(paste0("draw", seq_len(n_draws)),
                               rownames(M), colnames(M)))
  for (i in seq_len(n_draws)) {
    out[i, , ] <- M + matrix(rnorm(length(M), 0, sd), nrow = nrow(M),
                             dimnames = dimnames(M))
  }
  out
}

make_rho <- function(case_id, n_draws = 120) {
  mu <- c(temp = 0.56, precip = 0.36, forest = 0.58, soil = 0.27, pH = 0.20) +
    rnorm(length(linear_axes), 0, 0.02 + case_id * 0.001)
  draws <- matrix(NA_real_, n_draws, length(linear_axes),
                  dimnames = list(paste0("draw", seq_len(n_draws)), linear_axes))
  for (i in seq_along(linear_axes)) draws[, i] <- rnorm(n_draws, mu[i], 0.07)
  rho <- data.frame(axis = linear_axes,
                    mean = colMeans(draws),
                    lwr = apply(draws, 2, quantile, 0.025),
                    upr = apply(draws, 2, quantile, 0.975),
                    stringsAsFactors = FALSE)
  list(rho = rho, rho_draws = draws)
}

make_omega <- function(B) {
  D <- as.matrix(stats::dist(B[, linear_axes, drop = FALSE]))
  O <- exp(-D / max(D, na.rm = TRUE) * 2.2)
  diag(O) <- 1
  dimnames(O) <- list(species, species)
  O
}

make_omega_draws <- function(O, n_draws = 60) {
  out <- array(NA_real_, dim = c(n_draws, nrow(O), ncol(O)),
               dimnames = list(paste0("draw", seq_len(n_draws)), species, species))
  for (i in seq_len(n_draws)) {
    z <- O + matrix(rnorm(length(O), 0, 0.025), nrow = nrow(O))
    z <- (z + t(z)) / 2
    diag(z) <- 1
    out[i, , ] <- pmax(pmin(z, 1), -1)
  }
  out
}

make_population <- function(B, traits, case_id) {
  focal <- species[1:5]
  rows <- list()
  k <- 1L
  env_grid <- seq(-1.8, 1.8, length.out = 5)
  for (sp in focal) {
    for (pop_id in seq_len(6)) {
      pop <- paste0("pop", LETTERS[pop_id])
      genotype <- paste0("G", pop_id)
      pop_shift <- (pop_id - 3.5) * 0.045 + rnorm(1, 0, 0.012)
      for (axis in linear_axes) {
        env <- env_grid[pop_id %% length(env_grid) + 1]
        genetic_value <- traits[sp, "body_size"] * 0.18 + pop_shift + rnorm(1, 0, 0.03)
        species_beta <- B[sp, axis]
        beta <- species_beta + pop_shift + 0.10 * env + rnorm(1, 0, 0.025)
        rows[[k]] <- data.frame(
          species = sp,
          population = pop,
          axis = axis,
          species_beta = species_beta,
          beta = beta,
          local_adaptation = abs(pop_shift) + 0.07 + 0.01 * case_id,
          plasticity = abs(0.10 * env) + 0.06 + 0.006 * case_id,
          genetic_value = genetic_value,
          environment = env,
          response = beta + 0.26 * env + (pop_id - 3) * 0.06 * env,
          genotype = genotype,
          stringsAsFactors = FALSE
        )
        k <- k + 1L
      }
    }
  }
  do.call(rbind, rows)
}

make_timeseries <- function(B, traits, case_id) {
  rows <- list()
  k <- 1L
  for (sp in species[1:5]) {
    for (axis in linear_axes) {
      beta0 <- B[sp, axis]
      for (time in seq(0, 20, by = 2)) {
        climate <- 0.20 * sin(time / 4) + 0.01 * time
        drift <- 0.015 * time * ifelse(clade[sp] %in% c("A", "E"), 1, -0.6)
        rows[[k]] <- data.frame(
          species = sp,
          axis = axis,
          time = time,
          beta = beta0 + drift + climate * 0.15 + rnorm(1, 0, 0.018),
          trait_value = traits[sp, "body_size"] + 0.02 * time + rnorm(1, 0, 0.025),
          ecological_turnover = 0.28 + 0.025 * time + 0.04 * sin(time / 3),
          evolutionary_turnover = 0.12 + 0.018 * time + 0.03 * cos(time / 4),
          community_metric = 0.35 + 0.03 * time + beta0 * 0.12 + rnorm(1, 0, 0.02),
          stringsAsFactors = FALSE
        )
        k <- k + 1L
      }
    }
  }
  do.call(rbind, rows)
}

make_site_data <- function(B, case_id) {
  n_site <- 28
  site <- paste0("site", sprintf("%02d", seq_len(n_site)))
  XData <- data.frame(
    temp = scale_num(seq(-2, 2, length.out = n_site) + rnorm(n_site, 0, 0.10)),
    precip = scale_num(sin(seq(0, 2 * pi, length.out = n_site)) + rnorm(n_site, 0, 0.09)),
    forest = scale_num(rep(c(-1.2, -0.4, 0.6, 1.1), length.out = n_site) + rnorm(n_site, 0, 0.10)),
    soil = scale_num(cos(seq(0, 2 * pi, length.out = n_site)) + rnorm(n_site, 0, 0.08)),
    pH = scale_num(rep(seq(-1.4, 1.4, length.out = 7), 4) + rnorm(n_site, 0, 0.09)),
    history_index = scale_num(seq(1.2, -1.2, length.out = n_site) + rnorm(n_site, 0, 0.12))
  )
  rownames(XData) <- site
  coords <- data.frame(x = seq_len(n_site),
                       y = sin(seq_len(n_site) / 3) * 4 + case_id * 0.15,
                       row.names = site)
  X <- as.matrix(XData[, linear_axes])
  X2 <- X^2
  colnames(X2) <- quadratic_axes
  eta <- matrix(B[, "(Intercept)"], nrow = n_site, ncol = length(species),
                byrow = TRUE) +
    X %*% t(B[, linear_axes, drop = FALSE]) +
    X2 %*% t(B[, quadratic_axes, drop = FALSE])
  prob <- pmin(pmax(pnorm(eta), 0.02), 0.98)
  Y <- matrix(rbinom(n_site * length(species), 1, prob), nrow = n_site,
              dimnames = list(site, species))
  list(XData = XData, coords = coords, Y = Y)
}

pairwise_beta_phylo_local <- function(B, phylo) {
  sp <- rownames(B)
  pairs <- utils::combn(sp, 2)
  pd <- ape::cophenetic.phylo(phylo)
  out <- data.frame(
    species1 = pairs[1, ],
    species2 = pairs[2, ],
    beta_distance = as.numeric(sqrt(colSums((B[pairs[1, ], , drop = FALSE] -
      B[pairs[2, ], , drop = FALSE])^2))),
    phylo_distance = as.numeric(pd[pairs[1, ], pairs[2, ]]),
    stringsAsFactors = FALSE
  )
  out$beta_cosine_similarity <- vapply(seq_len(nrow(out)), function(i) {
    a <- B[out$species1[i], ]
    b <- B[out$species2[i], ]
    sum(a * b) / sqrt(sum(a^2) * sum(b^2))
  }, numeric(1))
  out
}

make_phylo_precomputed <- function(case_id, phylo) {
  axes_rep <- linear_axes
  global <- function(metric, base) {
    data.frame(axis = axes_rep, estimate = base + rnorm(length(axes_rep), 0, 0.025),
               lwr = pmax(base - 0.18, -1), upr = pmin(base + 0.18, 1),
               source = paste0("precomputed_", metric), stringsAsFactors = FALSE)
  }
  max_distance <- max(ape::cophenetic.phylo(phylo), na.rm = TRUE)
  correlogram <- do.call(rbind, lapply(axes_rep, function(axis) {
    d <- seq(max_distance * 0.06, max_distance * 0.94, length.out = 7)
    sim <- 0.75 * exp(-d / (max_distance * 0.38)) - 0.18 * (d / max_distance) +
      rnorm(length(d), 0, 0.025) + match(axis, axes_rep) * 0.01
    p_value <- pmin(0.99, pmax(0.001, exp(-abs(sim) * 6.2)))
    data.frame(axis = axis, phylo_distance = d, similarity = sim, moran_i = sim,
               lwr = sim - 0.15, upr = sim + 0.15,
               p_value = p_value,
               signif = ifelse(p_value < 0.001, "***",
                         ifelse(p_value < 0.01, "**",
                         ifelse(p_value < 0.05, "*", ""))),
               source = "precomputed_binned_moran_i", stringsAsFactors = FALSE)
  }))
  list(
    pagel_lambda = global("pagel_lambda", c(0.58, 0.34, 0.52, 0.24, 0.18)),
    blomberg_k = global("blomberg_k", c(0.62, 0.38, 0.60, 0.28, 0.20)),
    abouheif_cmean = global("abouheif_cmean", c(0.46, 0.30, 0.42, 0.22, 0.14)),
    phylogenetic_correlogram = correlogram
  )
}

make_precomputed <- function(B, G, traits, history, pop, ts, case_id, phylo) {
  edge_meta <- branch_time_table(phylo)
  root_age <- tree_root_age(phylo)
  nedge <- nrow(phylo$edge)
  branch_shifts <- expand.grid(branch = seq_len(nedge), axis = linear_axes,
                               KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  branch_shifts <- merge(branch_shifts, edge_meta, by = "branch", all.x = TRUE)
  branch_shifts$P_shift <- plogis(-1.25 + 0.10 * branch_shifts$branch +
    0.35 * branch_shifts$branch_length / max(branch_shifts$branch_length) +
    0.42 * (branch_shifts$axis == "temp") +
    0.34 * (branch_shifts$axis == "forest") +
    0.22 * (branch_shifts$axis == "pH") +
    rnorm(nrow(branch_shifts), 0, 0.10))
  branch_shifts$delta_theta <- rnorm(nrow(branch_shifts), 0, 0.20) +
    ifelse(branch_shifts$P_shift > 0.65, branch_shifts$P_shift * 0.85, 0)
  branch_shifts$axis_specific_shift <- branch_shifts$delta_theta
  branch_shifts$branch_contrast_outlier_score <- abs(branch_shifts$delta_theta) * 2.2
  branch_prob <- aggregate(P_shift ~ branch, branch_shifts, max, na.rm = TRUE)
  shift_count_posterior <- data.frame(
    draw = seq_len(800),
    number_of_niche_shifts = vapply(seq_len(800), function(i) {
      sum(stats::rbinom(nrow(branch_prob), size = 1, prob = branch_prob$P_shift))
    }, numeric(1)),
    source = "poisson_binomial_approximation_from_branch_P_shift",
    stringsAsFactors = FALSE
  )

  branch_rate_shifts <- edge_meta
  branch_rate_shifts$rate_ratio <- pmax(0.15,
    0.65 + 1.80 * branch_rate_shifts$branch_length / max(branch_rate_shifts$branch_length) +
      0.20 * sin(branch_rate_shifts$branch + case_id))
  branch_rate_shifts$P_rate_shift <- pmin(0.98, pmax(0.02,
    0.12 + abs(branch_rate_shifts$rate_ratio - 1) * 0.42 +
      0.05 * cos(branch_rate_shifts$time_start_Ma / 8)))
  branch_rate_shifts$source <- "precomputed_branch_rate_shift_probability"

  rate_shifts <- expand.grid(axis = linear_axes, clade = c("global", LETTERS[1:5]),
                             KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  clade_rate <- c(global = 1.00, A = 2.10, B = 1.35, C = 0.86, D = 0.62, E = 1.75)
  axis_rate <- c(temp = 1.18, precip = 1.05, forest = 0.95, soil = 0.88, pH = 1.10)
  rate_shifts$rate_ratio <- clade_rate[rate_shifts$clade] * axis_rate[rate_shifts$axis] +
    rnorm(nrow(rate_shifts), 0, 0.05)
  rate_shifts$axis_specific_rate <- 0.12 * rate_shifts$rate_ratio
  rate_shifts$P_rate_shift <- pmin(0.97, pmax(0.05, 0.18 + abs(rate_shifts$rate_ratio - 1) * 0.38))

  dtt <- do.call(rbind, lapply(linear_axes, function(axis) {
    time <- seq(0, root_age, length.out = 30)
    data.frame(axis = axis, time = time,
               disparity = 0.15 + 0.018 * time +
                 0.15 * sin(time / 12 + match(axis, linear_axes)) +
                 0.05 * match(axis, linear_axes),
               bm_mean = 0.12 + 0.014 * time,
               bm_lwr = 0.02 + 0.010 * time,
               bm_upr = 0.25 + 0.020 * time,
               stringsAsFactors = FALSE)
  }))
  mdi <- data.frame(axis = linear_axes,
                    MDI = c(1.10, 0.55, 0.70, -0.15, 0.25) + rnorm(length(linear_axes), 0, 0.03),
                    stringsAsFactors = FALSE)
  early_burst <- data.frame(axis = linear_axes,
                            early_burst_parameter = c(0.45, 0.15, -0.22, -0.55, 0.30) +
                              rnorm(length(linear_axes), 0, 0.03),
                            stringsAsFactors = FALSE)
  peaks <- data.frame(species = species,
                      peak = c("cold_wet", "cold_wet", "warm_dry",
                               "warm_dry", "warm_dry", "mesic",
                               "mesic", "open_dry", "open_dry",
                               "alpine", "alpine", "alkaline",
                               "cold_wet", "alkaline", "alkaline"),
                      stringsAsFactors = FALSE)
  pairwise <- pairwise_beta_phylo_local(B[, linear_axes, drop = FALSE], phylo)
  pairwise$convergence_score <- as.numeric(pairwise$phylo_distance >
    stats::quantile(pairwise$phylo_distance, 0.72) &
    pairwise$beta_distance < stats::quantile(pairwise$beta_distance, 0.30))
  pairwise$source <- "precomputed_convergence_screen"

  clade_gamma <- do.call(rbind, lapply(LETTERS[1:5], function(cl) {
    d <- as.data.frame(as.table(G[, linear_axes, drop = FALSE]), stringsAsFactors = FALSE)
    names(d) <- c("trait", "axis", "gamma")
    s <- c(A = 0.16, B = 0.08, C = -0.04, D = -0.10, E = 0.13)[cl]
    d$gamma <- d$gamma + s * c(temp = 1, precip = -0.4, forest = 0.6, soil = 0.35, pH = -0.5)[d$axis]
    d$clade <- cl
    d[, c("clade", "trait", "axis", "gamma")]
  }))
  regime_gamma <- do.call(rbind, lapply(regime_levels, function(rg) {
    d <- as.data.frame(as.table(G[, linear_axes, drop = FALSE]), stringsAsFactors = FALSE)
    names(d) <- c("trait", "axis", "gamma")
    d$gamma <- d$gamma + rnorm(nrow(d), 0, 0.08)
    d$regime <- rg
    d[, c("regime", "trait", "axis", "gamma")]
  }))
  gamma_shift_probability <- as.data.frame(as.table(abs(G[, linear_axes, drop = FALSE])),
                                           stringsAsFactors = FALSE)
  names(gamma_shift_probability) <- c("trait", "axis", "gamma_shift_probability")
  gamma_shift_probability$gamma_shift_probability <- pmin(0.98,
    0.10 + gamma_shift_probability$gamma_shift_probability * 0.85)
  gamma_branch_shift_probability <- edge_meta
  gamma_branch_shift_probability$gamma_shift_probability <- pmin(0.98, pmax(0.03,
    0.10 + 0.72 * branch_prob$P_shift[match(edge_meta$branch, branch_prob$branch)] +
      0.06 * sin(edge_meta$time_start_Ma / 6 + case_id)))
  gamma_branch_shift_probability$source <- "precomputed_branch_gamma_shift_probability"

  local_plasticity_partition <- aggregate(cbind(local_adaptation, plasticity) ~ species + axis,
                                          pop, mean)
  total <- local_plasticity_partition$local_adaptation +
    local_plasticity_partition$plasticity + 0.18
  local_plasticity_partition$species_mean_contribution <- 0.18 / total
  local_plasticity_partition$local_adaptation_contribution <-
    local_plasticity_partition$local_adaptation / total
  local_plasticity_partition$plasticity_contribution <-
    local_plasticity_partition$plasticity / total
  local_plasticity_partition$residual_contribution <- 0.18 / total
  local_plasticity_partition$LA_PL_ratio <-
    local_plasticity_partition$local_adaptation_contribution /
    local_plasticity_partition$plasticity_contribution

  feedback_strength <- data.frame(
    path = c("community_to_selection", "selection_to_traits",
             "traits_to_beta_divergence", "beta_divergence_to_community",
             "feedback_loop"),
    feedback_strength = c(0.46, 0.38, 0.41, 0.37, 0.29) + case_id * 0.002,
    source = "precomputed_path_coefficients_not_causal",
    stringsAsFactors = FALSE
  )
  contribution_time <- do.call(rbind, lapply(seq(0, 20, by = 2), function(t) {
    v <- c(trait_evolution = 0.25 + 0.03 * sin(t / 4),
           environment_climate = 0.30 + 0.04 * cos(t / 5),
           biotic_interactions = 0.20 + 0.02 * sin(t / 3),
           space_barriers = 0.13 + 0.02 * cos(t / 4),
           unexplained = 0.12)
    v <- v / sum(v)
    data.frame(time = t, component = names(v), contribution = as.numeric(v),
               stringsAsFactors = FALSE)
  }))

  trace <- do.call(rbind, lapply(c("beta_shift", "beta_rate", "niche_breadth",
                                   "gamma", "lambda_res"), function(par) {
    data.frame(draw = seq_len(400),
               parameter = par,
               value = cumsum(rnorm(400, 0, 0.015)) + rnorm(1, 0, 0.08),
               stringsAsFactors = FALSE)
  }))
  tree_sensitivity <- do.call(rbind, lapply(c("niche_shift_magnitude", "rate_ratio",
                                              "niche_breadth", "gamma",
                                              "residual_phylo_signal"), function(metric) {
    data.frame(metric = metric, value = rnorm(80, 0, 0.35 + runif(1, 0, 0.12)),
               tree = paste0("tree", seq_len(80)), stringsAsFactors = FALSE)
  }))
  simulation_recovery <- do.call(rbind, lapply(c("shift_magnitude", "rate_ratio",
                                                 "niche_breadth", "gamma",
                                                 "residual_phylo_signal"), function(metric) {
    true <- rnorm(180, 0, 0.85)
    data.frame(metric = metric, true_value = true,
               estimated_value = true * 0.92 + rnorm(length(true), 0, 0.25),
               stringsAsFactors = FALSE)
  }))
  block_cv <- expand.grid(block = paste0("Fold ", 1:5),
                          metric = c("AUC", "RMSE", "log_score", "R2"),
                          KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  block_cv$value <- c(0.89, 0.91, 0.87, 0.90, 0.88,
                      0.63, 0.58, 0.67, 0.60, 0.66,
                     -0.78, -0.72, -0.83, -0.74, -0.81,
                      0.61, 0.65, 0.59, 0.63, 0.60)
  posterior_predictive_check <- data.frame(
    metric = c("PPC_p_value", "PIT_KS_p_value", "mean_richness_error", "beta_dispersion_error"),
    value = c(0.34, 0.41, 0.08, 0.10),
    stringsAsFactors = FALSE
  )
  mcmc_diagnostics <- data.frame(
    metric = rep(c("ESS", "Rhat"), each = 5),
    parameter = rep(c("beta_shift", "beta_rate", "niche_breadth", "gamma", "lambda_res"), 2),
    value = c(980, 840, 760, 690, 720, 1.003, 1.006, 1.008, 1.010, 1.004),
    stringsAsFactors = FALSE
  )

  trait_omission_sensitivity <- calc_trait_mediation_metrics(
    hmsc_ecoevo(beta = B, gamma = G, traits = traits, history = history,
                phylo = phylo, clade = clade, regime = regime),
    axes = linear_axes
  )$trait_omission_sensitivity

  c(list(
    branch_shifts = branch_shifts,
    shift_count_posterior = shift_count_posterior,
    branch_rate_shifts = branch_rate_shifts,
    rate_shifts = rate_shifts,
    dtt = dtt,
    mdi = mdi,
    early_burst = early_burst,
    convergence = pairwise,
    peaks = peaks,
    clade_gamma = clade_gamma,
    regime_gamma = regime_gamma,
    gamma_shift_probability = gamma_shift_probability,
    gamma_branch_shift_probability = gamma_branch_shift_probability,
    local_plasticity_partition = local_plasticity_partition,
    reaction_norms = pop[, c("species", "population", "genotype", "environment", "response")],
    trait_omission_sensitivity = trait_omission_sensitivity,
    history_explained_r2 = data.frame(axis = linear_axes, predictor_set = "history",
                                      r2 = c(0.38, 0.32, 0.45, 0.28, 0.22),
                                      stringsAsFactors = FALSE),
    eco_evolutionary_feedback_path = data.frame(
      from = c("community_state_t", "selection_t1", "traits_t1",
               "beta_divergence_t1", "community_state_t1"),
      to = c("selection_t1", "traits_t1", "beta_divergence_t1",
             "community_state_t1", "selection_t2"),
      path = feedback_strength$path,
      weight = feedback_strength$feedback_strength,
      source = "precomputed_path_coefficients_not_causal",
      stringsAsFactors = FALSE
    ),
    beta_t_shift_rate = aggregate(beta ~ species + axis, ts, function(z) mean(abs(diff(z)))),
    trait_evolution_contribution = data.frame(
      trait = traits_names,
      contribution = c(0.16, 0.12, 0.10, 0.15, 0.08, 0.18, 0.11, 0.05, 0.06),
      stringsAsFactors = FALSE
    ),
    community_driven_selection_index = data.frame(species = species,
      index = scale_num(seq_along(species) + rnorm(length(species), 0, 0.6)),
      stringsAsFactors = FALSE),
    feedback_strength = feedback_strength,
    eco_evo_turnover_ratio = data.frame(metric = "eco_evo_turnover_ratio",
      value = mean(ts$evolutionary_turnover) / mean(ts$ecological_turnover),
      stringsAsFactors = FALSE),
    lagged_niche_shift = data.frame(lag = 1:3, correlation = c(0.36, 0.24, 0.11),
                                    stringsAsFactors = FALSE),
    interaction_mediated_evolution = data.frame(interaction = c("competition", "facilitation",
                                                                "enemy_release"),
      contribution = c(0.18, 0.11, 0.08), stringsAsFactors = FALSE),
    contribution_time = contribution_time,
    tree_uncertainty_sensitivity = tree_sensitivity,
    environment_omission_check = data.frame(
      metric = c("without_temperature", "without_precipitation", "without_forest",
                 "without_soil", "without_pH"),
      value = c(38, 42, 61, 29, 34), stringsAsFactors = FALSE),
    prior_sensitivity = data.frame(metric = c("beta_prior", "gamma_prior",
                                              "rho_prior", "omega_prior"),
                                   value = c(-4, -3, 5, 2),
                                   stringsAsFactors = FALSE),
    simulation_recovery = simulation_recovery,
    block_cross_validation = block_cv,
    posterior_predictive_check = posterior_predictive_check,
    mcmc_diagnostics = mcmc_diagnostics,
    trace = trace
  ), make_phylo_precomputed(case_id, phylo))
}

save_case <- function(case_id, case_name) {
  tree <- make_time_tree(case_id)
  root_age <- tree_root_age(tree)
  traits <- make_traits(case_id)
  history <- make_history(case_id)
  B <- make_beta(traits, history, case_id)
  G <- make_gamma(case_id)
  beta_draws <- make_draws(B, n_draws = 120, sd = 0.045 + case_id * 0.002)
  gamma_draws <- make_draws(G, n_draws = 120, sd = 0.035 + case_id * 0.001)
  rho_obj <- make_rho(case_id)
  omega <- make_omega(B)
  omega_draws <- make_omega_draws(omega)
  pop <- make_population(B, traits, case_id)
  ts <- make_timeseries(B, traits, case_id)
  site <- make_site_data(B, case_id)
  pre <- make_precomputed(B, G, traits, history, pop, ts, case_id, tree)

  evo <- hmsc_ecoevo(
    beta = B,
    beta_draws = beta_draws,
    gamma = G,
    gamma_draws = gamma_draws,
    rho = rho_obj$rho,
    rho_draws = rho_obj$rho_draws,
    omega = omega,
    omega_draws = omega_draws,
    traits = traits,
    phylo = tree,
    history = history,
    population = pop,
    timeseries = ts,
    XData = site$XData,
    Y = site$Y,
    coords = site$coords,
    distr = "probit",
    clade = clade,
    regime = regime,
    precomputed = pre,
    metadata = list(
      case_id = case_id,
      case_name = case_name,
      analysis_axes = analysis_axes,
      linear_axes = linear_axes,
      quadratic_map = quadratic_map,
      tree_time_unit = "Ma",
      root_age = root_age,
      figure_template = paste0("Figure ", case_id),
      note = "Complete deterministic toy data for exercising all HmscEcoEvo diagnostic modules. Advanced evolutionary tables are precomputed demonstration inputs, not empirical causal evidence."
    )
  )

  dir <- file.path(case_root, case_name)
  if (dir.exists(dir)) unlink(dir, recursive = TRUE)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(evo, file.path(dir, "evo.rds"))
  saveRDS(pre, file.path(dir, "precomputed.rds"))
  saveRDS(beta_draws, file.path(dir, "beta_draws.rds"))
  saveRDS(gamma_draws, file.path(dir, "gamma_draws.rds"))
  saveRDS(omega_draws, file.path(dir, "omega_draws.rds"))
  saveRDS(rho_obj$rho_draws, file.path(dir, "rho_draws.rds"))
  ape::write.tree(tree, file.path(dir, "phylo.tre"))
  utils::write.csv(B, file.path(dir, "beta.csv"))
  utils::write.csv(G, file.path(dir, "gamma.csv"))
  utils::write.csv(rho_obj$rho, file.path(dir, "rho.csv"), row.names = FALSE)
  utils::write.csv(traits, file.path(dir, "traits.csv"))
  utils::write.csv(history, file.path(dir, "history_species.csv"))
  utils::write.csv(omega, file.path(dir, "omega.csv"))
  utils::write.csv(pop, file.path(dir, "population.csv"), row.names = FALSE)
  utils::write.csv(ts, file.path(dir, "timeseries.csv"), row.names = FALSE)
  utils::write.csv(site$XData, file.path(dir, "XData.csv"))
  utils::write.csv(site$coords, file.path(dir, "coords.csv"))
  utils::write.csv(site$Y, file.path(dir, "Y.csv"))
  utils::write.csv(data.frame(species = names(clade), clade = unname(clade),
                              regime = unname(regime)),
                   file.path(dir, "species_groups.csv"), row.names = FALSE)
  utils::write.csv(data.frame(
    case_id = case_id,
    case_name = case_name,
    linear_axes = paste(linear_axes, collapse = ";"),
    analysis_axes = paste(analysis_axes, collapse = ";"),
    tree_time_unit = "Ma",
    root_age = root_age
  ), file.path(dir, "metadata.csv"), row.names = FALSE)
  invisible(dir)
}

dirs <- Map(save_case, seq_along(case_names), case_names)
cat("Generated ", length(dirs), " Nature-template HmscEcoEvo cases in:\n",
    case_root, "\n", sep = "")
