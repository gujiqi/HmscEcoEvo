make_toy_evo <- function() {
  beta <- matrix(c(
    0.8, -0.3, -0.25,
    0.4,  0.2, -0.15,
   -0.2,  0.6, -0.05,
   -0.6, -0.4, -0.20
  ), nrow = 4, byrow = TRUE)
  rownames(beta) <- paste0("sp", 1:4)
  colnames(beta) <- c("temp", "precip", "temp_sq")
  set.seed(1)
  beta_draws <- array(rnorm(60 * length(beta), rep(as.vector(beta), each = 60), 0.05),
                      dim = c(60, nrow(beta), ncol(beta)),
                      dimnames = list(NULL, rownames(beta), colnames(beta)))
  traits <- data.frame(size = c(1, 2, 2.5, 4), dispersal = c(0.2, 0.5, 0.7, 1.1),
                       row.names = rownames(beta))
  gamma <- matrix(c(0.3, -0.2, 0.1, 0.4, -0.1, 0.2),
                  nrow = 2, byrow = TRUE,
                  dimnames = list(colnames(traits), colnames(beta)))
  phy <- ape::read.tree(text = "((sp1:1,sp2:1):1,(sp3:1,sp4:1):1);")
  history <- data.frame(hist_age = c(1, 2, 3, 4), hist_loss = c(0.1, 0.2, 0.3, 0.4),
                        row.names = rownames(beta))
  pop <- data.frame(
    species = rep(rownames(beta), each = 4),
    population = rep(c("p1", "p2"), times = 8),
    axis = rep(c("temp", "precip"), each = 2, times = 4),
    beta = rep(c(0.7, 0.9, -0.2, -0.4), times = 4),
    local_adaptation = rep(c(0.6, 0.7, 0.2, 0.3), times = 4),
    plasticity = rep(c(0.3, 0.2, 0.4, 0.5), times = 4),
    genetic_value = rep(c(1, 2, 1, 2), times = 4),
    trait = rep("body_size", 16),
    trait_value = rep(c(1.0, 1.3, 0.8, 1.1), times = 4),
    genotype = rep(c("g1", "g2"), times = 8),
    environment = rep(c(-1, 1), times = 8),
    response = rep(c(0.2, 0.8), times = 8)
  )
  ts <- data.frame(species = "sp1", axis = "temp", time = 1:4,
                   beta = c(0.1, 0.2, 0.25, 0.4),
                   ecological_turnover = c(1, 1.2, 1.1, 1.4),
                   evolutionary_turnover = c(0.2, 0.3, 0.25, 0.4),
                   community_metric = c(0.1, 0.2, 0.3, 0.5))
  hmsc_ecoevo(
    beta = beta,
    beta_draws = beta_draws,
    gamma = gamma,
    rho = data.frame(axis = "temp", mean = 0.4, lwr = 0.1, upr = 0.7),
    traits = traits,
    phylo = phy,
    history = history,
    population = pop,
    timeseries = ts,
    clade = c(sp1 = "A", sp2 = "A", sp3 = "B", sp4 = "B"),
    regime = c(sp1 = "r1", sp2 = "r1", sp3 = "r2", sp4 = "r2"),
    XData = data.frame(temp = 1:4, precip = 4:1),
    coords = data.frame(x = 1:4, y = c(1, 1, 2, 2))
  )
}

test_that("hmsc_ecoevo constructor aligns core data", {
  evo <- make_toy_evo()
  expect_s3_class(evo, "hmsc_ecoevo")
  expect_equal(dim(evo$beta), c(4, 3))
  expect_equal(dim(evo$beta_draws)[2:3], c(4, 3))
  expect_equal(rownames(evo$traits), evo$species)
})

test_that("hmsc_ecoevo reorders Gamma by trait and axis names", {
  beta <- matrix(c(1, 2, 3, 4), nrow = 2,
                 dimnames = list(c("sp1", "sp2"), c("temp", "precip")))
  traits <- data.frame(size = c(1, 2), dispersal = c(0.2, 0.8),
                       row.names = rownames(beta))
  gamma <- matrix(c(10, 20, 30, 40), nrow = 2,
                  dimnames = list(c("dispersal", "size"),
                                  c("precip", "temp")))
  evo <- hmsc_ecoevo(beta = beta, traits = traits, gamma = gamma)
  expect_equal(rownames(evo$gamma), c("size", "dispersal"))
  expect_equal(colnames(evo$gamma), c("temp", "precip"))
  expect_equal(evo$gamma["size", "temp"], 40)
  expect_equal(evo$gamma["dispersal", "precip"], 10)
})

test_that("niche metrics include all basic Beta outputs", {
  evo <- make_toy_evo()
  m <- calc_niche_metrics(evo)
  expect_s3_class(m, "hmsc_niche_metrics")
  expect_true(all(c("beta", "support", "lwr", "upr") %in% names(m$beta_summary)))
  expect_true("beta_magnitude" %in% names(m$magnitude))
  expect_true(nrow(m$optimum_breadth) > 0)
  expect_equal(dim(m$cosine_similarity), c(4, 4))
  expect_equal(dim(m$angle), c(4, 4))
  expect_true(all(c("response_lwr", "response_upr") %in% names(m$response_curves)))
})

test_that("quadratic niche optimum and breadth use Gaussian-equivalent formulas", {
  beta <- matrix(c(1, -0.5), nrow = 1,
                 dimnames = list("sp1", c("temp", "temp_sq")))
  m <- calc_niche_metrics(hmsc_ecoevo(beta = beta, species = "sp1",
                                      axes = c("temp", "temp_sq")))
  ob <- m$optimum_breadth
  expect_equal(ob$optimum, 1)
  expect_equal(ob$breadth_index, 1)
})

test_that("P_shift is used for model-based niche shift counts", {
  evo <- make_toy_evo()
  pre <- list(branch_shifts = data.frame(
    branch = rep(seq_len(nrow(evo$phylo$edge)), each = 1),
    axis = "temp",
    P_shift = c(0.95, 0.10, 0.80, 0.20, 0.72, 0.65)
  ))
  m <- calc_evo_transition_metrics(evo, axes = "temp",
                                   p_shift_threshold = 0.7,
                                   precomputed = pre)
  expect_equal(m$shift_counts$number_of_niche_shifts[m$shift_counts$axis == "temp"], 3)
  expect_true(all(c("mean_P_shift", "max_P_shift", "sum_P_shift") %in% names(m$shift_support)))
})

test_that("shift-count posterior from branch probabilities is deterministic", {
  evo <- make_toy_evo()
  pre <- list(branch_shifts = data.frame(
    branch = seq_len(nrow(evo$phylo$edge)),
    axis = "temp",
    P_shift = c(0.95, 0.10, 0.80, 0.20, 0.72, 0.65)
  ))
  m1 <- calc_evo_transition_metrics(evo, axes = "temp", precomputed = pre)
  m2 <- calc_evo_transition_metrics(evo, axes = "temp", precomputed = pre)
  expect_identical(m1$shift_count_posterior, m2$shift_count_posterior)
  expect_equal(
    mean(m1$shift_count_posterior$number_of_niche_shifts),
    sum(pre$branch_shifts$P_shift),
    tolerance = 0.01
  )
  expect_true(all(m1$shift_count_posterior$number_of_niche_shifts >= 0))
})

test_that("evolutionary metrics include distance-corrected Beta similarity", {
  evo <- make_toy_evo()
  m <- calc_evo_transition_metrics(evo)
  expect_true("distance_corrected_similarity" %in% names(m))
  expect_true(all(c("expected_beta_cosine_similarity",
                    "distance_corrected_similarity",
                    "source") %in% names(m$distance_corrected_similarity)))
  expect_true(nrow(m$distance_corrected_similarity) > 0)
})

test_that("local phylogenetic signal remains matched to species names", {
  beta <- matrix(c(1.2, 0.8, -0.4, -0.9), ncol = 1,
                 dimnames = list(c("sp1", "sp2", "sp3", "sp4"), "temp"))
  phy <- ape::read.tree(text = "((sp3:1,sp4:1):1,(sp1:1,sp2:1):1);")
  evo1 <- hmsc_ecoevo(beta = beta, phylo = phy)
  evo2 <- hmsc_ecoevo(beta = beta[c("sp4", "sp2", "sp1", "sp3"), , drop = FALSE],
                      phylo = phy)
  l1 <- calc_phylo_signal_metrics(evo1)$local_signal
  l2 <- calc_phylo_signal_metrics(evo2)$local_signal
  l1 <- l1[order(l1$species), c("species", "local_moran_i")]
  l2 <- l2[order(l2$species), c("species", "local_moran_i")]
  rownames(l1) <- rownames(l2) <- NULL
  expect_equal(l1, l2, tolerance = 1e-10)
})

test_that("Gamma specialization varies by axis and honors axis filtering", {
  evo <- make_toy_evo()
  m <- calc_gamma_evolution_metrics(evo, axes = c("temp", "precip"))
  s <- m$trait_axis_specialization
  expect_false(any(s$axis == "temp_sq"))
  expect_false(any(s$axis == "(Intercept)"))
  expect_true(length(unique(round(s$trait_axis_specialization, 4))) > 2)
})

test_that("HMSC-like postList posterior draws are extracted and oriented", {
  species <- c("sp1", "sp2")
  axes <- c("(Intercept)", "temp", "precip", "forest", "temp_sq")
  traits <- data.frame(size = c(1, 2), mass = c(0.3, 0.7),
                       row.names = species)
  samples <- lapply(1:4, function(i) {
    list(
      Beta = matrix(seq_len(length(axes) * length(species)) / 10 + i / 100,
                    nrow = length(axes), ncol = length(species)),
      Gamma = matrix(seq_len(length(axes) * ncol(traits)) / 20 + i / 100,
                     nrow = length(axes), ncol = ncol(traits)),
      rho = 0.2 + i / 100,
      Lambda = list(matrix(c(0.4, 0.1, 0.2, 0.5) + i / 100, nrow = 2))
    )
  })
  model <- list(
    spNames = species,
    covNames = axes,
    trNames = colnames(traits),
    TrData = traits,
    XData = data.frame(site_env = c(1, 2, 3), row.names = paste0("site", 1:3)),
    Y = matrix(0, nrow = 3, ncol = 2,
               dimnames = list(paste0("site", 1:3), species)),
    postList = list(samples),
    ranLevels = list(list(xDim = 0))
  )
  hist <- list(
    XData_history = data.frame(hist_env = c(0.1, 0.2, 0.3),
                               row.names = paste0("site", 1:3)),
    TrData_history = data.frame(hist_trait = c(0.4, 0.8),
                                row.names = species)
  )
  class(hist) <- c("HmscEcoEvo_history_indices", "hmsc_history_indices")
  evo <- hmsc_ecoevo(model = model, history = hist)
  expect_equal(dim(evo$beta), c(2, 5))
  expect_equal(rownames(evo$beta), species)
  expect_equal(colnames(evo$beta), axes)
  expect_equal(dim(evo$beta_draws), c(4, 2, 5))
  expect_equal(dim(evo$gamma), c(2, 5))
  expect_equal(dim(evo$gamma_draws), c(4, 2, 5))
  expect_equal(nrow(evo$rho_draws), 4)
  expect_equal(dim(evo$omega_draws), c(4, 2, 2))
  expect_true("hist_env" %in% names(evo$XData))
  expect_true("hist_trait" %in% names(evo$history_species))
})

test_that("Beta orientation uses supplied axes when axes outnumber species", {
  B <- matrix(seq_len(10), nrow = 5,
              dimnames = list(c("(Intercept)", "temp", "precip", "forest", "temp_sq"),
                              c("sp1", "sp2")))
  evo <- hmsc_ecoevo(beta = B, axes = rownames(B))
  expect_equal(dim(evo$beta), c(2, 5))
  expect_equal(rownames(evo$beta), c("sp1", "sp2"))
  expect_equal(colnames(evo$beta), rownames(B))
})

test_that("all calc functions return metric objects", {
  evo <- make_toy_evo()
  expect_s3_class(calc_phylo_signal_metrics(evo), "hmsc_phylo_signal_metrics")
  expect_s3_class(calc_evo_transition_metrics(evo), "hmsc_evo_transition_metrics")
  expect_s3_class(calc_trait_mediation_metrics(evo), "hmsc_trait_mediation_metrics")
  expect_s3_class(calc_gamma_evolution_metrics(evo), "hmsc_gamma_evolution_metrics")
  expect_s3_class(calc_population_evolution_metrics(evo), "hmsc_population_evolution_metrics")
  expect_s3_class(calc_validation_metrics(evo), "hmsc_validation_metrics")
})

test_that("MCMC ESS diagnostics retain response parameter labels", {
  evo <- make_toy_evo()
  diagnostics <- calc_validation_metrics(evo)$mcmc_diagnostics
  ess <- diagnostics[diagnostics$metric == "ESS", , drop = FALSE]

  expect_equal(nrow(ess), prod(dim(evo$beta_draws)[-1]))
  expect_true(all(grepl("^beta_", ess$parameter)))
})

test_that("population metrics include population trait shift", {
  evo <- make_toy_evo()
  m <- calc_population_evolution_metrics(evo)
  expect_true("population_trait_shift" %in% names(m))
  expect_true(nrow(m$population_trait_shift) > 0)
  expect_true(all(c("population_trait_shift", "pairwise_trait_divergence") %in%
                    names(m$population_trait_shift)))
})

test_that("phylogenetic correlogram reports Moran's I rather than negative distance", {
  evo <- make_toy_evo()
  m <- calc_phylo_signal_metrics(evo)
  expect_true(all(c("moran_i", "similarity", "source") %in% names(m$correlogram)))
  expect_true(all(m$correlogram$source %in%
                    c("binned_phylogenetic_moran_i", "phylosignal::phyloCorrelogram")))
  expect_equal(m$correlogram$similarity, m$correlogram$moran_i)
})

test_that("phytools phylosig list results are converted to numeric estimates", {
  testthat::skip_if_not_installed("phytools")
  evo <- make_toy_evo()
  m <- calc_phylo_signal_metrics(evo)
  opt <- m$global_signal[m$global_signal$source == "phytools::phylosig", ,
                         drop = FALSE]
  expect_true(nrow(opt) > 0)
  expect_type(opt$estimate, "double")
  expect_true(all(is.finite(opt$estimate) | is.na(opt$estimate)))
})

test_that("plot functions return printable objects", {
  evo <- make_toy_evo()
  plots <- list(
    plot_niche_summary(evo),
    plot_phylo_signal(evo),
    plot_evo_transition(evo),
    plot_trait_mediation(evo),
    plot_gamma_evolution(evo),
    plot_population_evolution(evo),
    plot_validation_dashboard(evo)
  )
  for (p in plots) {
    expect_true(inherits(p, "ggplot") || inherits(p, "patchwork") ||
                  inherits(p, "hmsc_ecoevo_plot_list") ||
                  inherits(p, "hmsc_ecoevo_nature_page"))
  }
})

test_that("old hmscHist aliases remain backward compatible", {
  comm <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
  rownames(comm) <- c("site1", "site2")
  colnames(comm) <- c("sp1", "sp2")
  site_region <- data.frame(site = c("site1", "site2"), region = c("A", "A"))
  proj <- hmscHist_data(comm, site_region)
  expect_s3_class(proj, "HmscEcoEvo_data")
  expect_s3_class(proj, "hmscHist_data")
  srh <- data.frame(species = c("sp1", "sp2"), region = c("A", "A"),
                    entry_time = c(10, 20), retained = c(1, 1))
  proj <- build_history_tables(proj, species_region_history = srh, method = "manual")
  hist <- calc_history_indices(proj)
  expect_s3_class(hist, "HmscEcoEvo_history_indices")
  expect_true("colonization_age" %in% names(hist$XData_history))
})
