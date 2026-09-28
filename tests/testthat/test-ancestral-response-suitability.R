test_that("direct ancestral response reconstructs beta rather than averaging maps", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((A:5,B:5):7,C:12);")
  beta_draws <- array(
    c(
      2, 0.8, 0.1,
      -1, 0.2, 0.0,
      0.5, -0.4, 0.2,
      2.2, 0.7, 0.0,
      -0.8, 0.3, 0.1,
      0.6, -0.5, 0.1
    ),
    dim = c(2, 3, 3),
    dimnames = list(
      paste0("s", 1:2),
      c("A", "B", "C"),
      c("(Intercept)", "bio1", "bio12")
    )
  )

  response <- hee_evolution_ancestral_response_direct(
    tip_beta_draws = beta_draws,
    tree = tree,
    times = c(8, 0),
    basis_cols = c("bio1", "bio12"),
    branch_uncertainty = "none"
  )

  expect_true(all(c("lineage", "time_ma", "response_draw", "bio1", "bio12") %in%
                    names(response)))
  expect_false(any(response$lineage %in% c("A", "B") & response$time_ma == 8))
  expect_true(all(c("A", "B", "C") %in% response$lineage[response$time_ma == 0]))
  expect_equal(unique(response$intercept), 0)
  expect_match(unique(response$response_source),
               "direct_BM_ancestral_environmental_response",
               fixed = TRUE)

  modern_map <- data.frame(
    species = c("A", "B"),
    cell_id = c("c1", "c1"),
    suitability = c(0.8, 0.2)
  )
  expect_error(
    hee_evolution_ancestral_response_direct(
      tip_beta_draws = modern_map,
      tree = tree,
      times = 8,
      basis_cols = "bio1"
    ),
    "tip_beta_draws"
  )
})

test_that("ancestral environmental suitability uses ancestral beta and palaeo X", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((A:5,B:5):7,C:12);")
  beta_draws <- array(
    c(
      0, 1,
      0, 0,
      0, -1,
      0, 0.8,
      0, 0.2,
      0, -0.8
    ),
    dim = c(2, 3, 2),
    dimnames = list(
      paste0("s", 1:2),
      c("A", "B", "C"),
      c("(Intercept)", "bio1")
    )
  )
  earth <- expand.grid(
    cell_id = c("c1", "c2"),
    time_ma = c(8, 0),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  earth$region <- "R1"
  earth$land <- 1
  earth$bio1 <- ifelse(earth$cell_id == "c1", 1, -1)
  E <- hee_paleo_earth_state(earth, env_cols = "bio1")

  out <- hee_environmental_filtering_ancestral_suitability(
    earth_state = E,
    tip_beta_draws = beta_draws,
    tree = tree,
    times = c(8, 0),
    basis_cols = "bio1",
    link = "probit",
    branch_uncertainty = "none"
  )

  expect_s3_class(out, "hee_ancestral_suitability")
  expect_true(all(out$draws$suitability >= 0 & out$draws$suitability <= 1))
  expect_true(all(c("suitability_mean", "suitability_q025",
                    "suitability_q975", "n_response_draws") %in%
                    names(out$summary)))
  expect_false(any(out$draws$lineage %in% c("A", "B") &
                     out$draws$time_ma == 8))
  expect_true(all(out$summary$suitability_mean >= 0 &
                    out$summary$suitability_mean <= 1))
  for (model in c("OU", "EB")) {
    sensitivity <- hee_environmental_filtering_ancestral_suitability(
      earth_state = E, tip_beta_draws = beta_draws, tree = tree,
      times = c(8, 0), basis_cols = "bio1", link = "probit",
      evolution_model = model, branch_uncertainty = "none")
    expect_true(all(is.finite(sensitivity$draws$suitability)))
    expect_true(all(sensitivity$draws$suitability >= 0 &
                    sensitivity$draws$suitability <= 1))
    expect_equal(unique(attr(sensitivity$responses,
                             "model_comparison")$model), model)
  }
})

test_that("direct ancestral response supports quick subsets after dropping extra tree tips", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "(((A:3,B:3):2,(C:2,D:2):3):5,(E:4,F:4):6);")
  beta <- data.frame(
    species = c("A", "B", "C", "D"),
    response_draw = "s1",
    bio1 = c(0.2, 0.4, -0.1, -0.3),
    bio12 = c(1.0, 0.8, 0.2, 0.1),
    stringsAsFactors = FALSE
  )

  out <- hee_evolution_ancestral_response_direct(
    tip_beta_draws = beta,
    tree = tree,
    times = c(4, 0),
    basis_cols = c("bio1", "bio12"),
    branch_uncertainty = "none",
    drop_tree_extra = TRUE
  )

  expect_true(nrow(out) > 0)
  expect_true(all(is.finite(out$bio1)))
  expect_true(all(is.finite(out$bio12)))
  expect_false(any(out$lineage %in% c("E", "F")))
})

test_that("direct ancestral response keeps all near-ultrametric tips at 0 Ma", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((A:10.0000000,B:9.9999990):0.0000000,C:10.0000000);")
  beta <- data.frame(
    species = c("A", "B", "C"),
    response_draw = "s1",
    bio1 = c(0.1, 0.2, 0.3),
    stringsAsFactors = FALSE
  )

  out <- hee_evolution_ancestral_response_direct(
    tip_beta_draws = beta,
    tree = tree,
    times = 0,
    basis_cols = "bio1",
    branch_uncertainty = "none"
  )

  expect_setequal(out$lineage[out$time_ma == 0], c("A", "B", "C"))
})

test_that("polytomies use Brownian tree covariance instead of the tip mean", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((A:2,B:2,C:2):3,(D:3,E:3):2);")
  expect_false(ape::is.binary(tree))
  beta <- data.frame(species = tree$tip.label, response_draw = "s1",
    bio1 = c(0, 0, 6, 10, 10))
  out <- hee_evolution_ancestral_response_direct(
    beta, tree, times = c(5, 2, 0), basis_cols = "bio1",
    branch_uncertainty = "none")
  tip_cov <- ape::vcv.phylo(tree)
  weights <- solve(tip_cov, rep(1, length(tree$tip.label)))
  expected_root <- sum(weights * beta$bio1) / sum(weights)
  root <- out$bio1[out$lineage_type == "root"]
  expect_equal(root, expected_root, tolerance = 1e-8)
  expect_gt(abs(root - mean(beta$bio1)), .01)
  expect_true(all(is.finite(out$bio1)))
})

test_that("binary Brownian conditional nodes agree with ape ancestral means", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((A:5,B:5):7,C:12);")
  values <- c(A = 0, B = 2, C = 4)
  beta <- data.frame(species = names(values), response_draw = "s1",
    bio1 = as.numeric(values))
  out <- hee_evolution_ancestral_response_direct(
    beta, tree, times = 12, basis_cols = "bio1",
    branch_uncertainty = "none")
  fit <- ape::ace(values, tree, type = "continuous", method = "REML")
  expect_equal(out$bio1[out$lineage_type == "root"],
               unname(fit$ace[[1L]]), tolerance = 1e-4)
})

test_that("OU and EB covariance have their intended temporal meaning", {
  tree <- ape::read.tree(text = "((A:5,B:5):7,(C:5,D:5):7,E:12);")
  geom <- HmscEcoEvo:::.hee_response_tree_geometry(tree)
  bm <- HmscEcoEvo:::.hee_response_covariance(geom, "BM", 0)
  ou <- HmscEcoEvo:::.hee_response_covariance(geom, "OU", 1e-9)
  eb <- HmscEcoEvo:::.hee_response_covariance(geom, "EB", -1e-9)
  expect_equal(ou, bm, tolerance = 1e-6)
  expect_equal(eb, bm, tolerance = 1e-6)
  strong_ou <- HmscEcoEvo:::.hee_response_covariance(geom, "OU", 0.3)
  expect_lt(strong_ou[1, 2], bm[1, 2])
  expect_equal(strong_ou[1, 2],
    exp(-0.3 * 5)^2 * (-expm1(-2 * 0.3 * 7) / (2 * 0.3)),
    tolerance = 1e-8)
  early <- HmscEcoEvo:::.hee_response_interval_variance(0, 2, "EB", -0.2)
  late <- HmscEcoEvo:::.hee_response_interval_variance(10, 12, "EB", -0.2)
  expect_gt(early, late)
  expect_equal(early / late, exp(2), tolerance = 1e-8)
})

test_that("BM OU and EB reconstruct responses with auditable model fits", {
  tree <- ape::read.tree(text = "((A:5,B:5):7,(C:5,D:5):7,E:12);")
  beta <- data.frame(species = tree$tip.label, response_draw = "s1",
                     bio1 = c(-0.9, 0.2, 1.1, 0.7, -0.4),
                     bio12 = c(0.7, -0.2, 0.3, 1.2, -0.8))
  outputs <- lapply(c("BM", "OU", "EB", "AICc"), function(model) {
    hee_evolution_ancestral_response_direct(
      beta, tree, times = c(12, 8, 4, 0),
      basis_cols = c("bio1", "bio12"),
      evolution_model = model, branch_uncertainty = "none")
  })
  names(outputs) <- c("BM", "OU", "EB", "AICc")
  expect_true(all(vapply(outputs, function(x) all(is.finite(x$bio1)), logical(1))))
  expect_equal(nrow(attr(outputs$AICc, "model_comparison")), 6L)
  comparison <- attr(outputs$AICc, "model_comparison")
  expect_setequal(unique(comparison$model), c("BM", "OU", "EB"))
  expect_true(all(comparison$boundary_status %in%
    c("none", "lower_parameter_bound", "upper_parameter_bound",
      "constant_tip_values")))
  expect_equal(sum(comparison$selected), 2L)
  for (ax in c("bio1", "bio12")) {
    z <- comparison[comparison$axis == ax, ]
    expect_equal(z$model[z$selected], z$model[which.min(z$aicc)])
  }
  expect_true(all(comparison$parameter_per_ma[comparison$model == "OU"] > 0))
  expect_true(all(comparison$parameter_per_ma[comparison$model == "EB"] < 0))
  expect_equal(outputs$BM$bio1[outputs$BM$time_ma == 0],
               beta$bio1[match(outputs$BM$lineage[outputs$BM$time_ma == 0],
                               beta$species)])
  expect_gt(max(abs(outputs$OU$bio1 - outputs$BM$bio1)), 1e-6)
  eb_parameter <- attr(outputs$EB, "model_comparison")$parameter_per_ma[1L]
  if (abs(eb_parameter * 12) < 1e-4) {
    expect_equal(outputs$EB$bio1, outputs$BM$bio1, tolerance = 1e-5)
  } else {
    expect_gt(max(abs(outputs$EB$bio1 - outputs$BM$bio1)), 1e-6)
  }
})

test_that("flat tip responses retain a defined boundary diagnostic", {
  tree <- ape::read.tree(text = "((A:5,B:5):7,(C:5,D:5):7,E:12);")
  beta <- data.frame(species = tree$tip.label, response_draw = "s1",
                     bio1 = rep(0.2, 5))
  out <- hee_evolution_ancestral_response_direct(
    beta, tree, times = c(12, 6, 0), basis_cols = "bio1",
    evolution_model = "AICc", branch_uncertainty = "brownian_bridge",
    seed = 4)
  expect_equal(out$bio1, rep(0.2, nrow(out)), tolerance = 1e-8)
  fit <- attr(out, "model_comparison")
  expect_equal(fit$boundary_status[fit$model == "OU"],
               "constant_tip_values")
  expect_equal(fit$boundary_status[fit$model == "EB"],
               "constant_tip_values")
  expect_equal(fit$model[fit$selected], "BM")
})

test_that("model-specific stochastic bridges are reproducible and time coherent", {
  tree <- ape::read.tree(text = "((A:5,B:5):7,(C:5,D:5):7,E:12);")
  beta <- data.frame(species = tree$tip.label, response_draw = "s1",
                     bio1 = c(-1, 0.2, 1, 0.8, -0.4))
  for (model in c("BM", "OU", "EB")) {
    run <- function() hee_evolution_ancestral_response_direct(
      beta, tree, times = c(12, 10, 8, 6, 4, 2, 0), basis_cols = "bio1",
      evolution_model = model, branch_uncertainty = "brownian_bridge",
      seed = 42)
    first <- run()
    expect_equal(first, run())
    modern_name <- hee_evolution_ancestral_response_direct(
      beta, tree, times = c(12, 10, 8, 6, 4, 2, 0), basis_cols = "bio1",
      evolution_model = model, branch_uncertainty = "process_bridge",
      seed = 42)
    expect_equal(first, modern_name)
    expect_equal(first$bio1[first$time_ma == 0],
                 beta$bio1[match(first$lineage[first$time_ma == 0],
                                 beta$species)])
    expect_true(all(is.finite(first$bio1)))
  }
})
