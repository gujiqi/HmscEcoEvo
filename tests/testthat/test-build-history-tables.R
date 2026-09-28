test_that("manual history tables are standardized", {
  comm <- matrix(c(1,0,1, 0,1,1), nrow = 2, byrow = TRUE)
  rownames(comm) <- c("site1", "site2")
  colnames(comm) <- c("sp1", "sp2", "sp3")
  site_region <- data.frame(site = c("site1", "site2"), region = c("A", "B"))
  proj <- hmscHist_data(comm, site_region)
  srh <- data.frame(species = c("sp1", "sp2"), region = c("A", "B"), entry_time = c(10, 20))
  proj <- build_history_tables(proj, species_region_history = srh, method = "manual")
  expect_true("species_region_history" %in% names(proj$history_tables))
  expect_true("entry_prob" %in% names(proj$history_tables$species_region_history))
  expect_true("colonization_age" %in% proj$available_indices)
})

test_that("bsm-like events can be converted", {
  comm <- matrix(c(1,0,1, 0,1,1), nrow = 2, byrow = TRUE)
  rownames(comm) <- c("site1", "site2")
  colnames(comm) <- c("sp1", "sp2", "sp3")
  site_region <- data.frame(site = c("site1", "site2"), region = c("A", "B"))
  proj <- hmscHist_data(comm, site_region)
  bsm <- data.frame(
    event_type = c("dispersal", "loss", "speciation"),
    species = c("sp1", "sp1", "sp2"),
    from_region = c("A", NA, NA),
    to_region = c("B", NA, NA),
    region = c(NA, "B", "A"),
    time = c(30, 20, 40),
    event_prob = c(0.8, 0.6, 0.7),
    descendant_species = c(NA, NA, "sp2")
  )
  tip_ranges <- matrix(c(1,0,0, 0,1,0, 0,0,1), nrow = 3, byrow = TRUE)
  rownames(tip_ranges) <- c("sp1", "sp2", "sp3")
  colnames(tip_ranges) <- c("A", "B", "C")
  proj <- build_history_tables(proj, bsm_events = bsm, tip_ranges = tip_ranges)
  expect_true(all(c("species_region_history", "history_events", "speciation_events") %in% names(proj$history_tables)))
})
