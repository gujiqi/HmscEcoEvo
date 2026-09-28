test_that("process attribution decomposes CTMC occupancy sources exactly", {
  state <- data.frame(
    posterior_draw_id = 1,
    lineage_id = "sp1",
    cell_id = "c1",
    time_ma = c(1, 0),
    q = c(0.25, 0.4),
    lambda_C = c(0.3, 0.3),
    lambda_L = c(0.1, 0.1),
    cell_area_km2 = 10
  )

  out <- hee_process_attribution(state, target = "occupancy")

  budget <- out$state_budget
  persistent <- budget$estimate[budget$component == "persistent_occupancy"]
  new <- budget$estimate[budget$component == "new_colonisation"]
  recol <- budget$estimate[budget$component == "recolonisation"]

  r <- 0.3 + 0.1
  p01 <- 0.3 / r * (1 - exp(-r))
  p11 <- 0.3 / r + 0.1 / r * exp(-r)
  expect_equal(persistent, 0.25 * exp(-0.1), tolerance = 1e-8)
  expect_equal(new, (1 - 0.25) * p01, tolerance = 1e-8)
  expect_equal(recol, 0.25 * (p11 - exp(-0.1)), tolerance = 1e-8)
  expect_equal(persistent + new + recol,
               out$state_budget$estimate[
                 out$state_budget$component == "persistent_occupancy"] +
                 out$state_budget$estimate[
                   out$state_budget$component == "new_colonisation"] +
                 out$state_budget$estimate[
                   out$state_budget$component == "recolonisation"],
               tolerance = 1e-12)
  source_rows <- out$summary$contribution_type %in%
    c("retained_occupancy", "new_colonisation", "recolonisation")
  expect_equal(sum(out$summary$share_within_layer[source_rows], na.rm = TRUE),
               1, tolerance = 1e-8)
})

test_that("process attribution uses positive delta_t independent of input order", {
  state <- data.frame(
    posterior_draw_id = 1,
    lineage_id = "sp1",
    cell_id = "c1",
    time_ma = c(0, 2, 1),
    q = c(0.5, 0.1, 0.3),
    lambda_C = 0.2,
    lambda_L = 0.05,
    cell_area_km2 = 1
  )
  out <- hee_process_attribution(state, target = "occupancy")

  expect_true(all(out$time_series$delta_t > 0))
  expect_equal(sort(unique(out$time_series$time_start_ma), decreasing = TRUE),
               c(2, 1))
})

test_that("upstream processes are NA without a counterfactual simulator", {
  state <- data.frame(
    posterior_draw_id = 1,
    lineage_id = "sp1",
    cell_id = "c1",
    time_ma = c(1, 0),
    q = c(0.2, 0.3),
    lambda_C = 0.1,
    lambda_L = 0.05,
    cell_area_km2 = 1
  )
  out <- hee_process_attribution(state)
  up <- out$summary[out$summary$layer == "upstream_mechanism", ]

  expect_true(all(is.na(up$estimate)))
  expect_true(any(grepl("no_counterfactual_simulator", up$status)))
  expect_true(any(up$process_slug == "biotic_filtering" &
                    grepl("not_estimated", up$status)))
})

test_that("exact Shapley attribution splits interaction terms", {
  state <- data.frame(
    posterior_draw_id = 1,
    lineage_id = "sp1",
    cell_id = "c1",
    time_ma = c(1, 0),
    q = c(0.2, 0.3),
    lambda_C = 0.1,
    lambda_L = 0.05,
    cell_area_km2 = 1
  )
  simulator <- function(process_switches, reference, seed) {
    sw <- process_switches
    env <- isTRUE(sw[["environmental_filtering"]])
    disp <- isTRUE(sw[["dispersal"]])
    evo <- isTRUE(sw[["evolution"]])
    10 + (if (env) 1 else 0) +
      (if (disp) 2 else 0) +
      (if (evo) 4 else 0) +
      (if (env && disp) 3 else 0)
  }
  obj <- list(
    state = state,
    simulator = simulator,
    metadata = list(
      process_status = c(
        environmental_filtering = "data_estimated",
        dispersal = "data_estimated",
        evolution = "data_estimated",
        biotic_filtering = "not_estimated"
      )
    )
  )

  out <- hee_process_attribution(
    obj,
    upstream_processes = c("environmental_filtering", "dispersal",
                           "evolution")
  )
  up <- out$summary[out$summary$layer == "upstream_mechanism", ]
  est <- stats::setNames(up$estimate, up$process_slug)

  expect_equal(est[["environmental_filtering"]], 2.5, tolerance = 1e-8)
  expect_equal(est[["dispersal"]], 3.5, tolerance = 1e-8)
  expect_equal(est[["evolution"]], 4.0, tolerance = 1e-8)
})

test_that("lineage event deltas are reported separately from local loss", {
  state <- data.frame(
    posterior_draw_id = 1,
    lineage_id = "sp1",
    cell_id = "c1",
    time_ma = c(1, 0),
    q = c(0.2, 0.3),
    lambda_C = 0.1,
    lambda_L = 0.05,
    cell_area_km2 = 1
  )
  obj <- list(
    state = state,
    events = list(
      speciation_events = data.frame(time_ma = 0.8, delta_metric = 2),
      lineage_extinction_events = data.frame(time_ma = 0.2,
                                             delta_metric = -1)
    )
  )
  out <- hee_process_attribution(obj)
  ev <- out$summary[out$summary$layer == "lineage_event", ]

  expect_equal(ev$estimate[ev$process_slug == "speciation"], 2)
  expect_equal(ev$estimate[ev$process_slug == "extinction"], -1)
  expect_true(any(out$summary$contribution_type == "local_loss" &
                    out$summary$process_slug == "state_transition"))
})

test_that("targets that require extra structure fail explicitly", {
  state <- data.frame(
    posterior_draw_id = 1,
    lineage_id = "sp1",
    cell_id = "c1",
    time_ma = c(1, 0),
    q = c(0.2, 0.3),
    lambda_C = 0.1,
    lambda_L = 0.05,
    cell_area_km2 = 1
  )

  expect_error(
    hee_process_attribution(state, target = "phylogenetic_diversity"),
    "requires joint binary state draws"
  )
  expect_error(
    hee_process_attribution(state, target = "range_area",
                            normalise = "capacity"),
    "requires target_args\\$capacity"
  )
})

test_that("function catalog includes process attribution", {
  catalog <- hee_function_catalog()
  expect_true("hee_process_attribution" %in% catalog$canonical_name)
  expect_equal(catalog$layer[catalog$canonical_name ==
                               "hee_process_attribution"],
               "result_diagnostic")
})

test_that("interval-state attribution supports Case05 streaming state", {
  interval <- data.frame(
    posterior_draw_id = c(1, 1, 2, 2),
    lineage_id = c("sp1", "sp2", "sp1", "sp2"),
    cell_id = c("c1", "c1", "c2", "c2"),
    time_start_ma = c(10, 10, 10, 10),
    time_end_ma = c(5, 5, 5, 5),
    q0 = c(0.2, 0.6, 0.3, 0.5),
    q_observed_next = c(0.4, 0.5, 0.5, 0.4),
    lambda_C = c(0.08, 0.05, 0.1, 0.06),
    lambda_L = c(0.03, 0.04, 0.02, 0.05),
    cell_area_km2 = c(100, 100, 120, 120)
  )
  obj <- list(
    intervals = interval,
    metadata = list(
      process_status = c(
        environmental_filtering = "model_derived",
        dispersal = "scenario_parameterised",
        state_transition = "estimated_exact_ctmc_budget"
      )
    )
  )

  out <- hee_process_attribution(
    obj,
    target = "lineage_richness",
    integrate_time = TRUE,
    normalise = "time",
    return = c("summary", "time_series", "maps", "draws")
  )

  expect_s3_class(out, "hee_process_attribution")
  expect_true(isTRUE(out$metadata$interval_input))
  expect_true(any(out$summary$process_slug == "state_transition"))
  expect_true(any(out$summary$contribution_type == "local_loss"))
  expect_true(nrow(out$time_series) > 0)
  expect_true(nrow(out$maps) > 0)
})
