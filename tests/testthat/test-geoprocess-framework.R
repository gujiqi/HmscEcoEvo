expect_prob01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  expect_false(anyNA(x))
  expect_true(all(is.finite(x)))
  expect_true(all(x >= 0 & x <= 1))
}

test_that("network metrics detect fragmentation and corridor persistence", {
  conn <- expand.grid(
    from_region = c("A", "B", "C"),
    to_region = c("A", "B", "C"),
    time_ma = c(10, 0),
    stringsAsFactors = FALSE
  )
  conn$connectivity <- ifelse(conn$from_region == conn$to_region, 1, 0)
  conn$connectivity[conn$time_ma == 10 &
                      conn$from_region %in% c("A", "B") &
                      conn$to_region %in% c("A", "B")] <- 0.8
  conn$connectivity[conn$time_ma == 0] <- 0.8

  out <- hee_network_metrics(conn, threshold = 0.2)
  expect_true(all(c("node_metrics", "network_summary",
                    "corridor_persistence") %in% names(out)))
  expect_prob01(out$network_summary$fragmentation_index)
  expect_prob01(out$network_summary$largest_component_fraction)
  expect_prob01(out$corridor_persistence$corridor_persistence)

  old <- out$network_summary[out$network_summary$time_ma == 10, ]
  now <- out$network_summary[out$network_summary$time_ma == 0, ]
  expect_gt(old$fragmentation_index, now$fragmentation_index)
  expect_lt(old$largest_component_fraction, now$largest_component_fraction)
})

test_that("extinction layers separate local, regional, and lineage proxies", {
  dyn <- expand.grid(
    species = c("sp1", "sp2"),
    region = c("R1", "R2"),
    time_ma = c(5, 0),
    stringsAsFactors = FALSE
  )
  dyn$occupancy_probability <- c(0.9, 0.1, 0.8, 0.2, 0, 0, 0, 0)
  ext <- dyn[, c("species", "region", "time_ma")]
  ext$extinction_probability <- c(0.1, 0.7, 0.2, 0.8, 1, 1, 1, 1)

  out <- hee_extinction_layers(dyn, ext, threshold = 0.5)
  expect_prob01(out$local_extinction$local_extinction_risk)
  expect_prob01(out$regional_extinction$regional_extinction_proxy)
  expect_prob01(out$lineage_extinction$lineage_extinction_proxy)
  expect_true(any(out$lineage_extinction$lineage_extinction_proxy == 1))
})

test_that("cradle museum grave classification uses explicit region-time keys", {
  pool <- data.frame(
    region = c("A", "B", "C", "D"),
    time_ma = c(10, 10, 10, 10),
    species_pool_size = c(2, 10, 1, 0),
    origination = c(5, 0, 0, 0),
    regional_extinction = c(0, 0, 4, 0),
    source_sink_balance = c(0, 0.5, -0.5, 0),
    stringsAsFactors = FALSE
  )
  out <- hee_cradle_museum_grave(pool)
  expect_prob01(out$cradle_score)
  expect_prob01(out$museum_score)
  expect_prob01(out$grave_score)
  expect_equal(out$dominant_region_status[out$region == "A"], "cradle")
  expect_equal(out$dominant_region_status[out$region == "C"], "grave")
  expect_equal(out$dominant_region_status[out$region == "D"], "undetermined")

  dup <- rbind(pool, pool[1, ])
  expect_error(hee_cradle_museum_grave(dup), "Duplicate key")
})

test_that("model family and mechanism hypotheses expose M0-M8 and H1-H8", {
  fam <- hee_dynamic_model_family()
  hyp <- hee_geoprocess_hypotheses()
  expect_equal(fam$model_id, paste0("M", 0:8))
  expect_equal(hyp$hypothesis_id, paste0("H", 1:8))
  expect_true(all(c("required_components", "implementation_status") %in%
                    names(fam)))
  expect_true(all(c("framework_role", "core_process_mapping",
                    "interpretation_boundary") %in% names(fam)))
  expect_true(all(fam$framework_role == "scenario_diagnostic_not_core_process"))
  expect_true(any(grepl("not a speciation-rate model",
                        fam$core_process_mapping, fixed = TRUE)))
  expect_true(all(c("prediction", "proxy_outputs",
                    "interpretation_limit") %in% names(hyp)))
  expect_true("core_process_mapping" %in% names(hyp))
  expect_true(all(grepl("^Drivers|^evolution|^Post-event",
                        hyp$core_process_mapping)))
})

test_that("geoprocess model comparison flags impossible geography", {
  m1 <- data.frame(
    species = "sp1", cell_id = c("c1", "c2"), time_ma = 0,
    probability = c(0.3, 0.4), geographic_existence = c(1, 0)
  )
  m2 <- data.frame(
    species = "sp1", cell_id = c("c1", "c2"), time_ma = 0,
    probability = c(0.2, 0), geographic_existence = c(1, 0)
  )
  cmp <- hee_geoprocess_model_comparison(list(M1 = m1, M2 = m2))
  expect_equal(cmp$model_id, c("M1", "M2"))
  expect_gt(cmp$impossible_fraction[cmp$model_id == "M1"], 0)
  expect_equal(cmp$impossible_fraction[cmp$model_id == "M2"], 0)
  expect_error(
    hee_geoprocess_model_comparison(list(data.frame(probability = 0.5))),
    "named list"
  )
})
