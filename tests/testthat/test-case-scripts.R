read_case_script <- function(name) {
  path <- c(
    file.path("scripts", name),
    file.path("..", "..", "scripts", name)
  )
  path <- path[file.exists(path)]
  skip_if(length(path) == 0L, "case scripts are not installed in this test context")
  readLines(path[1], warn = FALSE)
}

extract_case_function <- function(lines, name) {
  start <- grep(paste0("^", name, " <- function"), lines)
  expect_length(start, 1L)
  balance <- 0L
  out <- character()
  for (i in seq.int(start, length(lines))) {
    out <- c(out, lines[i])
    balance <- balance +
      lengths(regmatches(lines[i], gregexpr("\\{", lines[i], fixed = FALSE))) -
      lengths(regmatches(lines[i], gregexpr("\\}", lines[i], fixed = FALSE)))
    if (length(out) > 1L && balance <= 0L) {
      break
    }
  }
  paste(out, collapse = "\n")
}

test_that("case scripts accept dashed and plain key-value CLI arguments", {
  case01 <- read_case_script("case01_simulated_full_540Ma.R")
  case02 <- read_case_script("case02_realdata_full_540Ma_template.R")

  expect_true(any(grepl('pat <- paste0("^(--)?", name, "=")', case01,
                        fixed = TRUE)))
  expect_true(any(grepl('pat <- paste0("^(--)?", key, "=")', case02,
                        fixed = TRUE)))
  expect_true(any(grepl('args %in% c(key, paste0("--", key))', case02,
                        fixed = TRUE)))
})

test_that("case scripts use safe joins for high-risk projection keys", {
  case01 <- read_case_script("case01_simulated_full_540Ma.R")
  case02 <- read_case_script("case02_realdata_full_540Ma_template.R")

  expect_true(any(grepl("safe_left_join <- function", case01, fixed = TRUE)))
  expect_true(any(grepl("safe_left_join <- function", case02, fixed = TRUE)))
  expect_true(any(grepl("label = \"region_cube\"", case02, fixed = TRUE)))
  expect_false(any(grepl("merge(suit, unique(region_grid", case02,
                         fixed = TRUE)))
})

test_that("Case 1 M1-M5 projection loop calls the package combination formula", {
  case01 <- read_case_script("case01_simulated_full_540Ma.R")
  expect_true(sum(grepl("do.call(hee_combine, combine_args)", case01,
                       fixed = TRUE)) >= 2)
  expect_false(any(grepl("P\\$probability <- P\\$probability \\* P\\$accessibility",
                         case01)))
})

test_that("Case 1 occupied-area summaries keep zero-area time slices", {
  case01 <- read_case_script("case01_simulated_full_540Ma.R")
  env <- new.env(parent = globalenv())
  eval(parse(text = extract_case_function(case01, "area_summary")), env)

  P_zero <- data.frame(
    cell_id = c("c1", "c2", "c1", "c2"),
    time_ma = c(540, 540, 0, 0),
    probability = c(0, 0.2, 0.4, 0.49),
    land_area_km2 = c(10, 20, 10, 20)
  )
  out_zero <- env$area_summary(P_zero, "M5_env_phylo_bgb_dynamic_dispersal")
  expect_equal(nrow(out_zero), 2L)
  expect_setequal(out_zero$time_ma, c(0, 540))
  expect_true(all(out_zero$land_area_km2 == 0))
  expect_true(all(out_zero$model_id == "M5_env_phylo_bgb_dynamic_dispersal"))

  P_partial <- P_zero
  P_partial$probability[P_partial$cell_id == "c1" & P_partial$time_ma == 0] <- 0.6
  out_partial <- env$area_summary(P_partial, "M4_env_phylo_bgb_static_dispersal")
  out_partial <- out_partial[order(out_partial$time_ma), ]
  expect_equal(out_partial$time_ma, c(0, 540))
  expect_equal(out_partial$land_area_km2, c(10, 0))
})

test_that("Case 5 no-BioGeoBEARS dispatcher documents its model contract", {
  case05 <- read_case_script("case05_plant200_global_dynamic_no_bgb.R")

  expect_true(any(grepl("forward_scenario", case05, fixed = TRUE)))
  expect_true(any(grepl("endpoint_conditioned_reconstruction", case05,
                        fixed = TRUE)))
  expect_true(any(grepl("all_processes_demo", case05, fixed = TRUE)))
  expect_true(any(grepl("case05_removed_bgb_inputs.csv", case05,
                        fixed = TRUE)))
  expect_true(any(grepl("case05_function_coverage.csv", case05,
                        fixed = TRUE)))
  expect_true(any(grepl("endpoint_conditioned_reconstruction requires an independent modern endpoint",
                        case05, fixed = TRUE)))
  expect_true(any(grepl("BioGeoBEARS/BSM is intentionally cancelled",
                        case05, fixed = TRUE)))
  expect_true(any(grepl("not_used_in_no_bgb_case05", case05, fixed = TRUE)))
  expect_true(any(grepl("empirical_core_off_zero_effect", case05,
                        fixed = TRUE)))
  expect_true(any(grepl("run_biotic_filtering", case05, fixed = TRUE)))
  expect_true(any(grepl("run_speciation_demo", case05, fixed = TRUE)))
  expect_true(any(grepl("run_lineage_extinction_demo", case05,
                        fixed = TRUE)))
  expect_true(any(grepl("case05_eight_process_inventory.csv", case05,
                        fixed = TRUE)))
  expect_true(any(grepl("biotic_competition_pressure", case05,
                        fixed = TRUE)))
  expect_true(any(grepl("speciation_inheritance_footprint", case05,
                        fixed = TRUE)))
  expect_true(any(grepl("scenario_lineage_extinction_pressure", case05,
                        fixed = TRUE)))
})

test_that("Case 5 all-processes engine exposes scenario process metrics", {
  engine <- read_case_script("case04_plant200_global_dynamic_no_bgb_final.R")

  expect_true(any(grepl("calc_biotic_filtering <- function", engine,
                        fixed = TRUE)))
  expect_true(any(grepl("run_biotic_filtering", engine, fixed = TRUE)))
  expect_true(any(grepl("biotic_competition_pressure", engine,
                        fixed = TRUE)))
  expect_true(any(grepl("biotic_facilitation_pressure", engine,
                        fixed = TRUE)))
  expect_true(any(grepl("speciation_event_summary.csv", engine,
                        fixed = TRUE)))
  expect_true(any(grepl("speciation_inheritance_footprint", engine,
                        fixed = TRUE)))
  expect_true(any(grepl("scenario_lineage_extinction_pressure", engine,
                        fixed = TRUE)))
  expect_true(any(grepl("SCENARIO_DEMO_NOT_EMPIRICAL", engine,
                        fixed = TRUE)))
})

test_that("Case 5 three-dispersal comparison preserves its controlled contract", {
  case05 <- read_case_script("case05_plant200_eight_process_three_dispersal.R")
  engine <- read_case_script("case04_plant200_global_dynamic_no_bgb_final.R")

  expect_true(any(grepl('"forward_geographic"', case05, fixed = TRUE)))
  expect_true(any(grepl('"topographic_limited"', case05, fixed = TRUE)))
  expect_true(any(grepl('"particle_topographic"', case05, fixed = TRUE)))
  expect_true(any(grepl("hee_dispersal_particle_kernel", engine, fixed = TRUE)))
  expect_true(any(grepl("hee_dispersal_rate_kernel", engine, fixed = TRUE)))
  expect_true(any(grepl("hee_projection_global_grid_step", engine, fixed = TRUE)))
  expect_true(any(grepl("max_hazard_per_step", engine, fixed = TRUE)))
  expect_true(any(grepl("persistence_reference_myr", engine, fixed = TRUE)))
  expect_true(any(grepl("diagnostic_points_per_time", engine, fixed = TRUE)))
  expect_true(any(grepl("case04_diagnostic_points_by_time.csv", engine,
                        fixed = TRUE)))
  expect_true(any(grepl("case04_diagnostic_point_process_trace.csv", engine,
                        fixed = TRUE)))
  expect_true(any(grepl("full_1deg_grid_used_for_model_and_maps", engine,
                        fixed = TRUE)))
  expect_true(any(grepl("case04_checkpoint_after_draw.rds", engine,
                        fixed = TRUE)))
  expect_true(any(grepl("case04_posterior_draw_tasks.csv", engine,
                        fixed = TRUE)))
  expect_true(any(grepl("streaming_complete_1deg_time_slices", engine,
                        fixed = TRUE)))
  expect_true(any(grepl("palaeo_slice_index", engine, fixed = TRUE)))
  expect_true(any(grepl("topographic_slice_index", engine, fixed = TRUE)))
  expect_true(any(grepl("load_earth_time <- function", engine, fixed = TRUE)))
  expect_true(any(grepl("movement_kernel_cache", engine, fixed = TRUE)))
  expect_true(any(grepl("utils::head(inside_times, n_time_slices)", engine,
                        fixed = TRUE)))
  expect_true(any(grepl("compact_shared_root_patch", engine, fixed = TRUE)))
  expect_true(any(grepl("case04_root_initialisation_audit.csv", engine,
                        fixed = TRUE)))
  expect_true(any(grepl("metric_storage", engine, fixed = TRUE)))
  expect_true(any(grepl("Random long-distance-dispersal destinations have been removed", engine,
                        fixed = TRUE)))
  expect_true(any(grepl("time_selection", engine, fixed = TRUE)))
  expect_true(any(grepl("evenly_spaced", engine, fixed = TRUE)))
  expect_true(any(grepl("selected_complete_1deg_time_slices.csv", engine,
                        fixed = TRUE)))
  expect_true(any(grepl("movement_kernel_cache_dir", engine, fixed = TRUE)))
  expect_true(any(grepl("same_inputs_same_tree_same_HMSC_draws", case05,
                        fixed = TRUE)))
  expect_true(any(grepl("case05_shared_metric_limits.csv", case05,
                        fixed = TRUE)))
  expect_true(any(grepl("case05_final_validation.csv", case05,
                        fixed = TRUE)))
  expect_true(any(grepl("case05_colonisation_zero_arrival_audit.csv", case05,
                        fixed = TRUE)))
  expect_true(any(grepl("case05_difference_state_codebook.csv", case05,
                        fixed = TRUE)))
  expect_true(any(grepl("response_history", case05, fixed = TRUE)))
  expect_true(any(grepl("comparison_state_tif", case05, fixed = TRUE)))
  expect_true(any(grepl("no-BioGeoBEARS", case05, fixed = TRUE)))
  expect_true(any(grepl("run_biotic_filtering", case05, fixed = TRUE)))
  expect_true(any(grepl("run_speciation_demo", case05, fixed = TRUE)))
  expect_true(any(grepl("run_lineage_extinction_demo", case05, fixed = TRUE)))
  expect_true(any(grepl("diagnostic_points_per_time", case05, fixed = TRUE)))
  expect_true(any(grepl("palaeo_slice_index", case05, fixed = TRUE)))
  expect_true(any(grepl("topographic_slice_index", case05, fixed = TRUE)))
})

test_that("posterior draw shard runner preserves full-grid posterior semantics", {
  runner <- read_case_script("case05_run_parallel_posterior_shards.R")
  merger <- read_case_script("case04_merge_posterior_draw_shards.R")

  expect_true(any(grepl("response_draw_ids", runner, fixed = TRUE)))
  expect_true(any(grepl("workers", runner, fixed = TRUE)))
  expect_true(any(grepl("shared_movement_kernel_cache", runner, fixed = TRUE)))
  expect_true(any(grepl("time_selection", runner, fixed = TRUE)))
  expect_true(any(grepl("diagnostic_points_per_time=40", runner, fixed = TRUE)))
  expect_true(any(grepl("make_png=false", runner, fixed = TRUE)))
  expect_true(any(grepl("streaming_complete_1deg_time_slices", merger,
                        fixed = TRUE)))
  expect_true(any(grepl("identical cell_id x time_ma keys", merger,
                        fixed = TRUE)))
  expect_true(any(grepl("No metric-table row was spatially sampled", merger,
                        fixed = TRUE)))
})

test_that("all Case01-Case07 R scripts parse and current cases use six-process APIs", {
  dirs <- c("scripts", file.path("..", "..", "scripts"))
  dirs <- dirs[dir.exists(dirs)]
  skip_if(length(dirs) == 0L, "case scripts are not installed in this test context")
  files <- list.files(dirs[[1L]], pattern = "^case0[1-7].*\\.R$",
                      full.names = TRUE)
  expect_gt(length(files), 0L)
  for (file in files) {
    expect_error(parse(file = file), NA, info = basename(file))
  }

  current <- c("case04_plant200_ancestral_suitability_no_bgb.R",
               "case05_v8_terrestrial_three_scheme_4deg.R",
               "case06_plant200_response_model_study.R",
               "case07_patch_evidence_ladder.R")
  for (name in current) {
    lines <- read_case_script(name)
    expect_false(any(grepl("hee_(colonisation|persistence)_", lines)),
                 info = name)
  }
})
