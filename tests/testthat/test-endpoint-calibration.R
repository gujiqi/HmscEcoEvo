test_that("occupancy-weighted diversity leaves an empty cell empty", {
  x <- rbind(
    occupied_even = c(0.5, 0.5),
    occupied_uneven = c(0.9, 0.1),
    empty = c(0, 0)
  )
  out <- hee_occupancy_weighted_diversity(x)

  expect_equal(out$occupancy_weighted_shannon_entropy[3], 0)
  expect_equal(out$occupancy_weighted_gini_simpson[3], 0)
  expect_equal(out$occupancy_weighted_effective_lineages_shannon[3], 0)
  expect_equal(out$occupancy_weighted_effective_lineages_simpson[3], 0)
  expect_false(out$diversity_defined[3])
  expect_gt(out$occupancy_weighted_shannon_entropy[1],
            out$occupancy_weighted_shannon_entropy[2])
  expect_equal(out$occupancy_weighted_effective_lineages_simpson[1], 2)
})

test_that("modern endpoint scoring reports calibration and richness agreement", {
  p <- rbind(c(0.9, 0.1), c(0.2, 0.8), c(0.1, 0.2))
  y <- rbind(c(1, 0), c(0, 1), c(0, 0))
  rownames(p) <- rownames(y) <- paste0("c", 1:3)
  colnames(p) <- colnames(y) <- c("sp1", "sp2")
  out <- hee_modern_endpoint_score(p, y, data_role = "spatial_holdout")

  expect_equal(out$summary$data_role, "spatial_holdout")
  expect_equal(out$summary$n_cells, 3)
  expect_equal(out$summary$n_species, 2)
  expect_lt(out$summary$brier_score, 0.1)
  expect_equal(nrow(out$species), 2)
  expect_equal(nrow(out$cells), 3)
})

test_that("endpoint particle weighting prefers the better aligned particle", {
  y <- matrix(c(1, 0, 1, 0), nrow = 2,
              dimnames = list(c("c1", "c2"), c("sp1", "sp2")))
  particles <- list(
    aligned = matrix(c(0.9, 0.1, 0.9, 0.1), nrow = 2,
                     dimnames = dimnames(y)),
    misaligned = matrix(c(0.1, 0.9, 0.1, 0.9), nrow = 2,
                        dimnames = dimnames(y))
  )
  out <- hee_endpoint_particle_weights(
    particles, y, data_role = "spatial_holdout"
  )

  expect_equal(sum(out$weight), 1)
  expect_gt(out$weight[out$particle_id == "aligned"],
            out$weight[out$particle_id == "misaligned"])
  expect_equal(out$scientific_role[[1L]], "endpoint_conditioned_particle_weight")
})
