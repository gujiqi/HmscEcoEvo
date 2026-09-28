test_that("BioGeoBEARS comparison and LRT are computed from supplied tables", {
  models <- data.frame(model = c("DEC", "DEC+J", "DIVALIKE"),
                       logLik = c(-20, -18, -22),
                       n_parameters = c(2, 3, 2))
  cmp <- hee_bgb_compare(models, n = 30)
  expect_true(all(c("AIC", "AICc", "delta_AICc", "AICc_weight") %in% names(cmp)))
  expect_equal(sum(cmp$AICc_weight), 1, tolerance = 1e-8)
  lrt <- hee_bgb_lrt(cmp[cmp$model == "DEC", ], cmp[cmp$model == "DEC+J", ])
  expect_true(lrt$LR > 0)
})

test_that("BioGeoBEARS AICc refuses undefined small-sample weights", {
  models <- data.frame(model = c("DEC", "DEC+J"),
                       logLik = c(-10, -9),
                       n_parameters = c(2, 3))
  expect_error(hee_bgb_compare(models, n = 4), "AICc is undefined")
  expect_error(hee_bgb_compare(models, n = NA_real_), "finite positive")
})

test_that("BioGeoBEARS comparison rejects duplicated candidate model names", {
  models <- data.frame(model = c("DEC", "DEC"),
                       logLik = c(-10, -9),
                       n_parameters = c(2, 2))
  expect_error(hee_bgb_compare(models, n = 20),
               "Duplicate BioGeoBEARS model name")
})

test_that("BioGeoBEARS LRT rejects non-nested standard model pairs", {
  dec <- data.frame(model = "DEC", logLik = -20, n_parameters = 2)
  diva <- data.frame(model = "DIVALIKE", logLik = -18, n_parameters = 2)
  expect_error(hee_bgb_lrt(dec, diva), "recognised nested pairs")
  bad_df <- data.frame(model = "DEC+J", logLik = -18, n_parameters = 2)
  expect_error(hee_bgb_lrt(dec, bad_df), "positive")
})

test_that("BSM event tables preserve the standard audit schema", {
  events <- data.frame(
    simulation_id = 1,
    model = "DEC",
    tree_id = "tree1",
    time_ma = 10,
    branch_or_node = "branch_1",
    event_type = "dispersal",
    from_area = "A",
    to_area = "B",
    ancestor_range = "A",
    descendant_range = "AB",
    stringsAsFactors = FALSE
  )
  out <- hee_bgb_bsm(events)
  expect_true(all(c("simulation_id", "model", "tree_id", "time_ma",
                    "branch_or_node", "event_type", "from_area", "to_area",
                    "ancestor_range", "descendant_range") %in% names(out)))
})

test_that("BioGeoBEARS accessibility model averaging normalises supplied weights", {
  acc <- list(
    DEC = data.frame(species = "sp1", region = "A", time_ma = 0,
                     accessibility = 0.8),
    DIVALIKE = data.frame(species = "sp1", region = "A", time_ma = 0,
                          accessibility = 0.6)
  )
  out <- hee_bgb_accessibility(acc, weights = c(DEC = 10, DIVALIKE = 1))
  expect_equal(out$accessibility, (10 * 0.8 + 1 * 0.6) / 11,
               tolerance = 1e-12)
  expect_true(out$accessibility >= 0 && out$accessibility <= 1)
})

test_that("BioGeoBEARS accessibility data.frame input is validated and clipped", {
  acc <- data.frame(species = "sp1", region = "A", time_ma = 0,
                    accessibility = 2)
  out <- hee_bgb_accessibility(acc)
  expect_equal(out$accessibility, 1)

  alt <- hee_bgb_accessibility(data.frame(species = "sp1", region = "A",
                                          time_ma = 0, A_BGB = -1))
  expect_equal(alt$accessibility, 0)

  dup <- data.frame(species = c("sp1", "sp1"), region = "A", time_ma = 0,
                    accessibility = c(0.2, 0.8))
  expect_error(hee_bgb_accessibility(dup), "Duplicate key")

  bad_time <- data.frame(species = "sp1", region = "A", time_ma = -1,
                         accessibility = 0.5)
  expect_error(hee_bgb_accessibility(bad_time), "negative ages")
})

test_that("BioGeoBEARS accessibility rejects missing or duplicated model weights", {
  acc <- list(
    DEC = data.frame(species = "sp1", region = "A", time_ma = 0,
                     accessibility = 0.8),
    DIVALIKE = data.frame(species = "sp1", region = "A", time_ma = 0,
                          accessibility = 0.6)
  )
  expect_error(hee_bgb_accessibility(acc, weights = c(DEC = 1)),
               "missing model")
  bad <- acc
  bad$DEC <- rbind(bad$DEC, bad$DEC)
  expect_error(hee_bgb_accessibility(bad), "Duplicate key")
})
