test_that("basic workflow computes indices without event tables", {
  comm <- matrix(c(1,1,0, 0,1,1), nrow = 2, byrow = TRUE)
  rownames(comm) <- c("site1", "site2")
  colnames(comm) <- c("sp1", "sp2", "sp3")
  site_region <- data.frame(site = c("site1", "site2"), region = c("A", "B"))
  srh <- data.frame(
    species = c("sp1", "sp2", "sp2", "sp3"),
    region = c("A", "A", "B", "B"),
    entry_time = c(80, 60, 50, 30),
    source_region = c("C", "B", "A", "A"),
    entry_prob = c(1, 1, 1, 1),
    retained = c(1, 1, 1, 0),
    loss_time = c(NA, NA, NA, 10)
  )

  proj <- hmscHist_data(comm, site_region)
  proj <- build_history_tables(proj, species_region_history = srh, method = "manual")
  hist <- calc_history_indices(proj)

  expect_true("colonization_age" %in% names(hist$XData_history))
  expect_true("prob_weighted_source_diversity" %in% names(hist$XData_history))
  expect_true("dominant_source_region" %in% names(hist$random_history))
  expect_true("lineage_retention_index" %in% names(hist$XData_history))
})

test_that("as_hmsc_xdata combines env and selected variables", {
  comm <- matrix(c(1,1,0, 0,1,1), nrow = 2, byrow = TRUE)
  rownames(comm) <- c("site1", "site2")
  colnames(comm) <- c("sp1", "sp2", "sp3")
  site_region <- data.frame(site = c("site1", "site2"), region = c("A", "B"))
  env <- data.frame(temp = c(1,2), row.names = c("site1", "site2"))
  srh <- data.frame(
    species = c("sp1", "sp2", "sp2", "sp3"),
    region = c("A", "A", "B", "B"),
    entry_time = c(80, 60, 50, 30),
    source_region = c("C", "B", "A", "A"),
    entry_prob = c(1, 1, 1, 1),
    retained = c(1, 1, 1, 0),
    loss_time = c(NA, NA, NA, 10)
  )
  proj <- hmscHist_data(comm, site_region, env = env)
  proj <- build_history_tables(proj, species_region_history = srh, method = "manual")
  hist <- calc_history_indices(proj)
  X <- as_hmsc_xdata(hist, env = env, strategy = "minimal")
  expect_true("temp" %in% names(X))
  expect_equal(rownames(X), rownames(comm))
})

test_that("gain loss turnover predictors are non-complementary summaries", {
  comm <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
  rownames(comm) <- c("site1", "site2")
  colnames(comm) <- c("sp1", "sp2")
  site_region <- data.frame(site = c("site1", "site2"), region = c("A", "B"))
  history_events <- data.frame(
    event_type = c("dispersal", "dispersal", "loss", "loss"),
    from_region = c("B", "B", NA, NA),
    to_region = c("A", "A", NA, NA),
    region = c(NA, NA, "A", "B"),
    event_prob = c(2, 3, 1, 10)
  )

  proj <- hmscHist_data(comm, site_region)
  proj <- build_history_tables(proj, history_events = history_events, method = "manual")
  hist <- calc_history_indices(proj, compute = "occupancy_dynamics")
  X <- hist$XData_history

  expect_equal(X["site1", "historical_gain_fraction"], 1)
  expect_equal(X["site2", "historical_gain_fraction"], 0)
  expect_equal(X["site1", "historical_loss_fraction"], 1 / 11)
  expect_equal(X["site2", "historical_loss_fraction"], 10 / 11)
  expect_equal(X["site1", "range_loss_fraction"], 1 / 6)
  expect_equal(X["site2", "range_loss_fraction"], 1)
  expect_equal(X["site1", "historical_turnover_predictor"], 6 / 16)
  expect_equal(X["site2", "historical_turnover_predictor"], 10 / 16)
  expect_false(all(X$historical_turnover_predictor == 1))
})

test_that("as_hmsc_xdata skips zero-variance history variables and falls back", {
  comm <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
  rownames(comm) <- c("site1", "site2")
  colnames(comm) <- c("sp1", "sp2")
  site_region <- data.frame(site = c("site1", "site2"), region = c("A", "A"))
  env <- data.frame(temp = c(5, 7), row.names = rownames(comm))
  srh <- data.frame(
    species = c("sp1", "sp2"),
    region = c("A", "A"),
    entry_time = c(10, 20),
    entry_prob = c(1, 1),
    retained = c(1, 1)
  )

  proj <- hmscHist_data(comm, site_region, env = env)
  proj <- build_history_tables(proj, species_region_history = srh, method = "manual")
  hist <- calc_history_indices(proj)
  X <- as_hmsc_xdata(hist, env = env, strategy = "minimal")

  expect_true("colonization_age" %in% names(X))
  expect_false("lineage_retention_index" %in% names(X))
  expect_false("historical_gain_fraction" %in% names(X))
})
