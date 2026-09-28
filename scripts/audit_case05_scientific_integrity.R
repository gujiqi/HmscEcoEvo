#!/usr/bin/env Rscript

# Case05 scientific-integrity audit
#
# This script audits an already-computed Case05 output without modifying its
# posterior shards, maps, or checkpoints. It is deliberately conservative:
# failed gates mean that the output is a process-scenario diagnostic, not a
# palaeodistribution reconstruction suitable for ecological interpretation.

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  key <- paste0("--", name, "=")
  hit <- args[startsWith(args, key)]
  if (!length(hit)) return(default)
  sub(paste0("^", key), "", hit[[length(hit)]])
}

case_root <- arg_value(
  "case_root",
  "C:/Users/Google/Documents/HMSC-HIST/outputs/HmscEcoEvo/case05_plant200_arias_tardis_preflight1_fullgrid_20260922"
)
transport_index <- arg_value(
  "transport_index",
  "C:/Users/Google/Documents/HMSC-HIST/HmscEcoEvo/derived_inputs/case05_paleomap_h3_grid_transport_325_0Ma_5Myr_1deg_v2_targetcentred/paleomap_h3_transport_index.csv"
)
output <- arg_value(
  "output",
  "C:/Users/Google/Documents/HMSC-HIST/outputs/HmscEcoEvo/case05_scientific_audit_20260922"
)

dir.create(output, recursive = TRUE, showWarnings = FALSE)
write_csv <- function(x, file) utils::write.csv(x, file, row.names = FALSE, na = "")
read_csv <- function(file) utils::read.csv(file, check.names = FALSE,
                                           stringsAsFactors = FALSE)
issue <- function(id, severity, status, component, finding, consequence,
                  required_correction, evidence) {
  data.frame(
    issue_id = id,
    severity = severity,
    status = status,
    component = component,
    finding = finding,
    scientific_consequence = consequence,
    required_correction = required_correction,
    evidence = evidence,
    stringsAsFactors = FALSE
  )
}

if (!dir.exists(case_root)) stop("Missing Case05 root: ", case_root)
if (!file.exists(transport_index)) stop("Missing transport index: ", transport_index)

scenario_dirs <- list.dirs(file.path(case_root, "02_scenario_runs"),
                           recursive = FALSE, full.names = TRUE)
if (!length(scenario_dirs)) stop("No scenario directories under: ", case_root)

transport <- read_csv(transport_index)
required_transport <- c("time_from_ma", "time_to_ma", "n_source_land_cells",
                        "n_target_land_cells", "n_mapped_source_cells",
                        "n_mapped_target_cells", "source_coverage_fraction",
                        "target_coverage_fraction")
missing_transport <- setdiff(required_transport, names(transport))
if (length(missing_transport)) {
  stop("Transport index is missing columns: ", paste(missing_transport, collapse = ", "))
}
transport$source_coverage_fraction <- as.numeric(transport$source_coverage_fraction)
transport$target_coverage_fraction <- as.numeric(transport$target_coverage_fraction)
transport$unmapped_source_fraction <- 1 - transport$source_coverage_fraction
transport$unmapped_target_fraction <- 1 - transport$target_coverage_fraction
transport$source_conservative_gate <- ifelse(
  is.finite(transport$source_coverage_fraction) & transport$source_coverage_fraction >= 0.99,
  "PASS", "FAIL"
)
transport$target_interpolation_gate <- ifelse(
  is.finite(transport$target_coverage_fraction) & transport$target_coverage_fraction >= 0.99,
  "PASS", "FAIL"
)
write_csv(transport, file.path(output, "plate_transport_interval_audit.csv"))

# A target-normalised interpolation table may look complete from the target
# perspective while being unusable as forward occupancy transport.  Audit the
# response to a one-cell source pulse: its occupied area should be carried to
# approximately one source-cell area (apart from documented submergence), not
# deleted or copied into hundreds of target cells.
transport_source_rows <- list()
transport_example_rows <- list()
for (ii in seq_len(nrow(transport))) {
  tr_file <- transport$file[[ii]]
  if (!is.character(tr_file) || !file.exists(tr_file)) next
  tr_object <- readRDS(tr_file)
  links <- if (is.list(tr_object) && !is.null(tr_object$transport)) {
    tr_object$transport
  } else {
    tr_object
  }
  links <- as.data.frame(links, stringsAsFactors = FALSE)
  required_links <- c("source_cell_id", "target_cell_id", "target_weight")
  if (!all(required_links %in% names(links))) next
  links$source_cell_id <- as.character(links$source_cell_id)
  links$target_cell_id <- as.character(links$target_cell_id)
  links$target_weight <- as.numeric(links$target_weight)
  source_area <- if ("source_area_km2" %in% names(links)) {
    stats::aggregate(as.numeric(links$source_area_km2),
                     list(source_cell_id = links$source_cell_id),
                     function(x) x[[which(is.finite(x))[1L]]])
  } else {
    data.frame(source_cell_id = unique(links$source_cell_id), x = 1)
  }
  names(source_area)[[2L]] <- "source_area_km2"
  target_area <- if ("target_area_km2" %in% names(links)) {
    as.numeric(links$target_area_km2)
  } else {
    rep(1, nrow(links))
  }
  by_source <- stats::aggregate(
    cbind(
      outgoing_weight = links$target_weight,
      n_target_cells = rep(1, nrow(links)),
      carried_area_km2 = links$target_weight * target_area
    ),
    list(source_cell_id = links$source_cell_id), sum
  )
  by_source <- merge(source_area, by_source, by = "source_cell_id", all = TRUE)
  by_source$n_target_cells[is.na(by_source$n_target_cells)] <- 0
  by_source$outgoing_weight[is.na(by_source$outgoing_weight)] <- 0
  by_source$carried_area_km2[is.na(by_source$carried_area_km2)] <- 0
  by_source$area_transfer_ratio <- by_source$carried_area_km2 /
    pmax(by_source$source_area_km2, .Machine$double.eps)
  # Add the genuinely unmapped source cells, which otherwise do not appear in
  # the link table and are the principal source-loss failure.
  source_ids <- if (is.list(tr_object) && !is.null(tr_object$source_cells)) {
    as.character(tr_object$source_cells)
  } else {
    unique(links$source_cell_id)
  }
  missing_source <- setdiff(source_ids, by_source$source_cell_id)
  if (length(missing_source)) {
    by_source <- rbind(
      by_source,
      data.frame(
        source_cell_id = missing_source,
        source_area_km2 = NA_real_,
        outgoing_weight = 0,
        n_target_cells = 0,
        carried_area_km2 = 0,
        area_transfer_ratio = 0,
        stringsAsFactors = FALSE
      )
    )
  }
  finite_ratio <- by_source$area_transfer_ratio[is.finite(by_source$area_transfer_ratio)]
  transport_source_rows[[length(transport_source_rows) + 1L]] <- data.frame(
    time_from_ma = transport$time_from_ma[[ii]],
    time_to_ma = transport$time_to_ma[[ii]],
    n_source_cells = length(source_ids),
    n_unmapped_source_cells = sum(by_source$n_target_cells == 0),
    fraction_unmapped_source_cells = mean(by_source$n_target_cells == 0),
    n_source_cells_with_multiple_targets = sum(by_source$n_target_cells > 1),
    fraction_source_cells_with_multiple_targets = mean(by_source$n_target_cells > 1),
    max_targets_from_one_source = max(by_source$n_target_cells),
    q50_targets_from_one_source = as.numeric(stats::quantile(by_source$n_target_cells, 0.50)),
    q99_targets_from_one_source = as.numeric(stats::quantile(by_source$n_target_cells, 0.99)),
    q50_area_transfer_ratio = if (length(finite_ratio)) as.numeric(stats::quantile(finite_ratio, 0.50)) else NA_real_,
    q99_area_transfer_ratio = if (length(finite_ratio)) as.numeric(stats::quantile(finite_ratio, 0.99)) else NA_real_,
    max_area_transfer_ratio = if (length(finite_ratio)) max(finite_ratio) else NA_real_,
    stringsAsFactors = FALSE
  )
  if (transport$time_from_ma[[ii]] %in% c(40, 5)) {
    top <- by_source[order(-by_source$n_target_cells, -by_source$area_transfer_ratio), , drop = FALSE]
    top <- utils::head(top, 25L)
    top$time_from_ma <- transport$time_from_ma[[ii]]
    top$time_to_ma <- transport$time_to_ma[[ii]]
    transport_example_rows[[length(transport_example_rows) + 1L]] <- top
  }
}
if (length(transport_source_rows)) {
  transport_source_audit <- do.call(rbind, transport_source_rows)
  write_csv(transport_source_audit,
            file.path(output, "plate_transport_source_pulse_audit.csv"))
}
if (length(transport_example_rows)) {
  write_csv(do.call(rbind, transport_example_rows),
            file.path(output, "plate_transport_source_pulse_examples.csv"))
}

time_rows <- list()
metric_rows <- list()
config_rows <- list()
for (scenario_dir in scenario_dirs) {
  scenario_id <- basename(scenario_dir)
  summary_file <- file.path(scenario_dir, "22_report", "case04_final_time_summary.csv")
  metric_file <- file.path(scenario_dir, "16_state_space_inference", "cell_time_all_metrics.csv")
  config_file <- file.path(scenario_dir, "00_config", "case04_final_config.csv")
  if (file.exists(summary_file)) {
    x <- read_csv(summary_file)
    x$scenario_id <- scenario_id
    time_rows[[scenario_id]] <- x
  }
  if (file.exists(metric_file)) {
    x <- read_csv(metric_file)
    x$scenario_id <- scenario_id
    metric_rows[[scenario_id]] <- x
  }
  if (file.exists(config_file)) {
    x <- read_csv(config_file)
    x$scenario_id <- scenario_id
    config_rows[[scenario_id]] <- x
  }
}
if (!length(time_rows)) stop("No Case05 final time summaries were found.")
time_summary <- do.call(rbind, time_rows)
write_csv(time_summary, file.path(output, "case05_time_summary_audit.csv"))
if (length(config_rows)) write_csv(do.call(rbind, config_rows),
                                   file.path(output, "case05_config_audit.csv"))

anchor_cols <- c("time_ma", "scenario_id",
                 "global_mean_expected_lineage_richness",
                 "global_mean_forward_raw_expected_lineage_richness",
                 "global_mean_terminal_anchor_difference_expected_lineage_richness",
                 "global_mean_terminal_anchor_weight",
                 "global_mean_hmsc_fixed_effect_nowcast_richness")
anchor_cols <- intersect(anchor_cols, names(time_summary))
anchor <- time_summary[time_summary$time_ma %in% c(10, 5, 0), anchor_cols,
                       drop = FALSE]
write_csv(anchor, file.path(output, "terminal_anchor_discontinuity_audit.csv"))

metric_semantics <- data.frame(
  metric = c(
    "relative_environmental_support_sum",
    "expected_lineage_richness",
    "binary_lineage_richness",
    "shannon_diversity",
    "simpson_diversity",
    "mean_arrival_probability_reporting_interval",
    "mean_colonisation_probability_reporting_interval",
    "mean_persistence_probability",
    "mean_local_extinction_probability",
    "scenario_lineage_extinction_pressure",
    "speciation_inheritance_footprint"
  ),
  valid_interpretation = c(
    "Intercept-free ancestral environmental-response support score; diagnostic only.",
    "Sum of occupancy probabilities over active sampled-surviving lineages.",
    "Count above an explicitly chosen occupancy threshold; threshold-sensitive.",
    "Occupancy-weighted entropy across selected active lineages, not observed community Shannon diversity.",
    "Occupancy-weighted Gini-Simpson quantity across selected active lineages, not observed community Simpson diversity.",
    "Hazard converted to probability over the declared reporting interval; source-limited process diagnostic.",
    "Arrival times establishment hazard converted to probability over the declared reporting interval; scenario-dependent.",
    "Conditional probability of local retention over its stated biological reference interval.",
    "One minus conditional local retention; not global lineage extinction.",
    "Optional scenario-only pressure diagnostic; not an estimated historical lineage-extinction rate.",
    "Tree-node inheritance diagnostic; not an inferred geographic speciation location."
  ),
  prohibited_interpretation = c(
    "Expected species richness or final palaeo-occupancy.",
    "Full historical plant species richness or observed palaeocommunity richness.",
    "Observed species richness.",
    "Observed palaeocommunity Shannon diversity.",
    "Observed palaeocommunity Simpson diversity.",
    "Empirically estimated historical dispersal probability.",
    "Empirically estimated historical colonisation probability.",
    "Empirically estimated historical persistence probability without calibration data.",
    "Historical global extinction probability.",
    "Historical extinction rate.",
    "Historical speciation location or rate."
  ),
  display_requirement = c(
    "Separate diagnostic atlas and its own scale; never share a richness legend.",
    "Raw GIS layer plus scientific scale and a companion adaptive inspection scale.",
    "State threshold and sensitivity; never substitute for the expected value.",
    "Use only as an occupancy-weighted diagnostic with a separate legend.",
    "Use only as an occupancy-weighted diagnostic with a separate legend.",
    "Use a process-specific scale, not 0-1 if values occupy a narrow range.",
    "Use a process-specific scale, not 0-1 if values occupy a narrow range.",
    "State the reference interval in title and legend.",
    "State the reference interval in title and legend.",
    "Keep out of empirical-core result figures.",
    "Keep out of geographic diversification claims."
  ),
  stringsAsFactors = FALSE
)
write_csv(metric_semantics, file.path(output, "metric_semantics_and_display_audit.csv"))

lineage_rows <- list()
for (scenario_dir in scenario_dirs) {
  active_file <- file.path(scenario_dir, "03_lineage_time_tree", "active_lineage_counts.csv")
  if (file.exists(active_file)) {
    x <- read_csv(active_file)
    x$scenario_id <- basename(scenario_dir)
    lineage_rows[[basename(scenario_dir)]] <- x
  }
}
if (length(lineage_rows)) {
  lineage <- do.call(rbind, lineage_rows)
  lineage$temporal_domain_interpretation <-
    "Before a named modern tip originated, a species-specific map must be NA; its ancestral lineage may have a valid lineage map."
  write_csv(lineage, file.path(output, "active_lineage_time_audit.csv"))
}

map_scale_rows <- list()
if (length(metric_rows)) {
  metrics <- do.call(rbind, metric_rows)
  numeric_names <- names(metrics)[vapply(metrics, is.numeric, logical(1))]
  excluded <- c("time_ma", "lon", "lat", "H_state", "land_area_km2")
  numeric_names <- setdiff(numeric_names, excluded)
  for (nm in numeric_names) {
    values <- metrics[[nm]]
    if (!is.numeric(values) || !any(is.finite(values))) next
    map_scale_rows[[nm]] <- data.frame(
      metric = nm,
      observed_min = min(values, na.rm = TRUE),
      observed_q01 = as.numeric(stats::quantile(values, 0.01, na.rm = TRUE)),
      observed_q50 = as.numeric(stats::quantile(values, 0.50, na.rm = TRUE)),
      observed_q99 = as.numeric(stats::quantile(values, 0.99, na.rm = TRUE)),
      observed_max = max(values, na.rm = TRUE),
      n_finite = sum(is.finite(values)),
      stringsAsFactors = FALSE
    )
  }
}
if (length(map_scale_rows)) write_csv(do.call(rbind, map_scale_rows),
                                      file.path(output, "map_value_distribution_audit.csv"))

anchor_issue <- if (all(c("global_mean_terminal_anchor_difference_expected_lineage_richness",
                          "global_mean_forward_raw_expected_lineage_richness") %in% names(time_summary))) {
  zero <- time_summary[time_summary$time_ma == 0, , drop = FALSE]
  ratio <- with(zero,
                global_mean_terminal_anchor_difference_expected_lineage_richness /
                  pmax(global_mean_forward_raw_expected_lineage_richness, .Machine$double.eps))
  any(is.finite(ratio) & ratio > 0.25)
} else FALSE

findings <- do.call(rbind, list(
  issue(
    "C05-001", "critical", "FAIL", "terminal endpoint",
    "The final 0 Ma occupancy map includes a convex HMSC fixed-effect terminal anchor.",
    "A last-slice change cannot be attributed to Earth history, dispersal, colonisation, persistence, or evolution.",
    "Set terminal_anchor_mode=none for all process trajectories. Retain the 0 Ma HMSC map only as a separate training diagnostic. Implement likelihood-based endpoint conditioning with held-out or independent data before claiming reconstruction.",
    if (anchor_issue) "Anchor contribution exceeds 25% of the raw forward endpoint in at least one scenario." else "Inspect terminal-anchor audit."
  ),
  issue(
    "C05-002", "critical", "FAIL", "plate carriage",
    "The supplied PALEOMAP target-centred transport tables are target-normalised interpolation tables, not source-conservative occupancy transport.",
    "Unmapped source land cells can disappear solely because no outgoing plate link was created, producing artificial continental cut-offs and local loss.",
    "Do not use these tables for forward occupancy carriage. Build source-complete forward plate-overlap links, require source coverage >=0.99 after verified land-loss accounting, and verify source-area / occupancy conservation on synthetic fields.",
    paste0("Source coverage ranges from ",
           formatC(min(transport$source_coverage_fraction, na.rm = TRUE), digits = 3, format = "f"),
           " to ", formatC(max(transport$source_coverage_fraction, na.rm = TRUE), digits = 3, format = "f"),
           "; current engine only gates target coverage.")
  ),
  issue(
    "C05-013", "critical", "FAIL", "plate carriage directionality",
    "The target-centred nearest-track table maps one source cell to multiple target cells while leaving many other source cells without an outgoing link.",
    "A local occupied source patch can be copied to many cells or deleted before dispersal. This is neither continental carriage nor biological movement, and it can create striped, abruptly truncated, or unrealistically redistributed occupancy maps.",
    "Do not transpose or renormalise this table ad hoc. Rebuild transport from plate polygons or validated source-to-target overlap, with explicit source-area weights, true land loss, and a conservation test on uniform and local-pulse fields.",
    "See plate_transport_source_pulse_audit.csv; any source with multiple targets or zero targets invalidates this table as a forward occupancy operator."
  ),
  issue(
    "C05-003", "critical", "FAIL", "inferential design",
    "The runner is a forward deterministic mean-field scenario with a hand-set root patch and hand-set demographic/dispersal rates, not an endpoint-conditioned ancestral-range inference.",
    "It cannot be presented as an Arias-style reconstruction or a posterior palaeodistribution estimate.",
    "Rename retained runs as forward process scenarios. For inferential output, fit terminal-range likelihood / stochastic mapping or SMC with explicit root and rate priors, using independent or cross-fitted endpoint information and fossil likelihoods.",
    "Configuration declares forward_scenario_matrix_streaming_not_particle_smoothing."
  ),
  issue(
    "C05-004", "high", "FAIL", "P2 calibration",
    "Movement, establishment, persistence, root range, and LDD parameters are scenario constants rather than data-calibrated estimates.",
    "Low arrival and colonisation can create nearly empty maps even on valid land, while a final anchor hides endpoint mismatch.",
    "Use explicit priors and sensitivity grids, calibrate against held-out modern ranges / spatial blocks and fossil clade support, and report posterior or scenario envelopes rather than a single parameter set.",
    "The default source emigration, establishment, persistence, root occupancy, and root patch size are fixed in the runner configuration."
  ),
  issue(
    "C05-005", "high", "FAIL", "Arias/TARDIS interpretation",
    "Arias diffusion and TARDIS least-cost paths are being used as ingredients of a generic forward occupancy CTMC without their respective conditioning/inference layers.",
    "The method is not equivalent to either article and must not use their names as a claim of methodological equivalence.",
    "Use Arias-style spherical diffusion inside a terminal-range/tree likelihood or stochastic-mapping module. Use TARDIS paths as branch-conditioned route diagnostics after spatial node locations are inferred, not as an uncalibrated global occupancy update.",
    "See local source PDFs: 扩散函数.pdf and s41559-025-02739-y.pdf."
  ),
  issue(
    "C05-006", "high", "FAIL", "metric semantics",
    "Relative environmental support is plotted beside diversity and occupancy products despite being an intercept-free response score.",
    "Users can incorrectly expect it to match dynamic expected lineage richness at 0 Ma.",
    "Move it to a P1 diagnostic section with its own scale and label. Do not call it richness or compare it numerically to final occupancy.",
    "Runner configuration calls it diagnostic only; map atlas still groups it with richness metrics."
  ),
  issue(
    "C05-007", "high", "FAIL", "map communication",
    "A common 0-200 richness scale and common 0-1 probability scale visually flatten maps whose dynamic values occupy a tiny fraction of those ranges.",
    "Valid but small signals look like empty land or broken rasters, masking whether an ecological pattern exists.",
    "Keep raw GIS values and an absolute comparison atlas, but add a second, clearly labelled per-metric robust-quantile inspection atlas. Never apply a shared richness scale to non-richness scores.",
    "Map-value distribution audit quantifies the mismatch."
  ),
  issue(
    "C05-008", "medium", "NEEDS_IMPLEMENTATION", "lineage identity",
    "Aggregate maps retain active lineages, but no user-facing species-to-ancestral-lineage time index is emitted.",
    "A map may be mistaken for a modern species before that tip exists; zero may be confused with outside-lineage-time NA.",
    "Write species_ancestral_path.csv and lineage_time_index.csv. Pre-tip species maps must be NA with the contemporaneous ancestral lineage ID, while aggregate active-lineage maps remain valid.",
    "Active lineage counts are available, but the requested identity tables are absent from scenario output."
  ),
  issue(
    "C05-009", "medium", "NEEDS_IMPLEMENTATION", "cladogenesis",
    "Copy-then-diverge is a valid default scenario, but its spatial inheritance is not constrained by independent historical geographic evidence in the no-BioGeoBEARS model.",
    "Speciation-footprint maps cannot be interpreted as inferred geographic speciation locations.",
    "Keep only as a labelled inheritance diagnostic or move to sensitivity analyses; do not include it in empirical-core conclusions.",
    "Runner process status itself labels the speciation layer as optional demonstration diagnostic."
  ),
  issue(
    "C05-010", "medium", "NEEDS_IMPLEMENTATION", "diversity scope",
    "Shannon, Simpson, PD, and functional layers are calculated from selected sampled-surviving lineage occupancies rather than a complete fossil-inclusive flora.",
    "They are not historical global plant-community diversity estimates.",
    "Rename all aggregate products as sampled-surviving-lineage diversity, calculate nonlinear metrics inside posterior state draws, and retain a clear sampled-tree temporal-domain NA mask before the root age.",
    "The dated tree contains selected extant tips and their ancestral branches only."
  ),
  issue(
    "C05-014", "high", "FAIL", "tree-time integration",
    "The dynamic loop uses only available palaeo-environment slice times. Dated-tree node ages and fossil age boundaries are not inserted into the process timeline.",
    "A split occurring between two 5 Ma environmental layers is applied at the next layer, so daughter lineages can inherit, disperse, and persist for too long or too short an interval. This can create abrupt last-slice changes in lineage counts and diversity.",
    "Construct master_time = unique(environment times union tree-node ages union fossil age boundaries). Interpolate only the declared Earth covariates between environmental layers, apply lineage identity changes at exact tree-node ages, and retain a node-timing audit.",
    "case05_plant200_arias_tardis.R constructs inside_times directly from all_env_times and then maps current lineages to the previous environmental slice."
  ),
  issue(
    "C05-011", "high", "FAIL", "nonlinear diversity metrics",
    "Shannon and Simpson are calculated after normalising cell-level lineage occupancies even where total expected occupancy is effectively zero.",
    "A cell with many numerically tiny probabilities can display high entropy despite having no meaningful reconstructed assemblage.",
    "Mask nonlinear composition metrics below a declared expected-lineage-richness support threshold, calculate them within posterior state draws, and show the support mask alongside every map.",
    "The map-value audit can show high Shannon/Simpson quantiles while expected lineage richness has a zero median."
  ),
  issue(
    "C05-012", "high", "FAIL", "disabled process maps",
    "The main scenario emits maps for optional speciation and lineage-extinction diagnostics even when their process modules are disabled, producing all-zero rasters.",
    "An all-zero map can be misread as evidence that no speciation or extinction occurred.",
    "Do not write or render a process map when its module is disabled. Write a one-row NOT_RUN metadata record instead; retain scenario-only maps in a separate appendix when explicitly enabled.",
    "The map-value audit reports zero maxima for speciation_inheritance_footprint and scenario_lineage_extinction_pressure."
  )
))
write_csv(findings, file.path(output, "case05_scientific_integrity_findings.csv"))

summary_lines <- c(
  "# Case05 科学完整性审计",
  "",
  "本审计不删除或改变 Case05 的任何结果、posterior shard、地图或 checkpoint。它只判断当前运行是否可以解释为古分布重建。",
  "",
  "## 结论",
  "",
  "**不要恢复或发布已经暂停的正式运行。** 当前 target-centred 板块表和 0 Ma 终点锚定，使过程层面的解释失效。",
  "",
  "## 必须完成的重构顺序",
  "",
  "1. 用来源完整、面积感知的正向板块搬运取代 target-centred 插值；除非能证实陆地沉没，否则拒绝任何来源格网无解释丢失的时间区间。",
  "2. 从动态轨迹中移除凸组合 0 Ma 终点锚定。HMSC 现代预测只能作为独立的 P1 训练诊断图。",
  "3. 实现终点条件推断：使用末端范围、空间留出记录和化石给历史轨迹加权，而不是混合最后一张地图。",
  "4. 将 Arias 球面扩散作为拟合或条件化的系统地理过程；在节点空间位置推断完成后，才用 TARDIS 最小成本路径做路线诊断。",
  "5. 校准或显式包络根范围、移动、建立和存续参数；单一参数组合只能称为敏感性情景。",
  "6. 使用 environment ∪ tree nodes ∪ fossil boundaries 的主时间轴；分化必须在精确节点时间发生。",
  "7. 分开指标语义和色标；增加物种到祖先谱系的时间索引，并在根年龄之前和物种形成之前明确输出 NA。",
  "",
  "## 输出文件",
  "",
  "- `case05_scientific_integrity_findings.csv`：按严重程度排列的问题和必需修正。",
  "- `plate_transport_interval_audit.csv`：每个时间区间的来源和目标覆盖率。",
  "- `plate_transport_source_pulse_audit.csv`：单一来源脉冲在当前板块表下是否被删除或复制的守恒审计。",
  "- `terminal_anchor_discontinuity_audit.csv`：10、5、0 Ma 的终点组成与跳变。",
  "- `map_value_distribution_audit.csv`：每个地图指标的原始数值分布。",
  "- `metric_semantics_and_display_audit.csv`：每个指标可写与不可写的解释及制图要求。",
  "",
  "## 科学边界",
  "",
  "本次被审计的 Case05 输出仅保留为正向过程情景诊断。它们不是经过验证的历史植物分布、扩散率、定殖率、存续率、灭绝率或古多样性重建。"
)
writeLines(summary_lines, file.path(output, "README_zh.md"), useBytes = TRUE)

cat("Case05 scientific-integrity audit written to:\n", normalizePath(output, winslash = "/", mustWork = FALSE), "\n", sep = "")
