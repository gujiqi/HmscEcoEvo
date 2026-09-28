#!/usr/bin/env Rscript

# Native-time Plant200 maps. Refugia use corrected Case05 patch results;
# movement origins/routes are forward-predictive, not posterior edge use.

`%||%` <- function(x, y) if (is.null(x)) y else x
parse_args <- function(args) {
  answer <- list()
  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) next
    z <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    answer[[z[[1L]]]] <- paste(z[-1L], collapse = "=")
  }
  answer
}
args <- parse_args(commandArgs(trailingOnly = TRUE))
index_only <- identical(args$index_only, "true")
script <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
if (!length(script)) stop("Run this case with Rscript.", call. = FALSE)
pkg_root <- normalizePath(file.path(dirname(script[[1L]]), ".."), winslash = "/")
root <- normalizePath(file.path(pkg_root, ".."), winslash = "/")
out <- normalizePath(args$output %||% file.path(root, "outputs", "HmscEcoEvo",
  paste0("case07_patch_evidence_ladder_", format(Sys.Date(), "%Y%m%d"))),
  winslash = "/", mustWork = TRUE)
comparison <- file.path(root, "outputs", "HmscEcoEvo",
  "case05_v8_four_scheme_comparison_20260926")
corrected <- file.path(comparison, "10_patch_history_results_corrected")
terrain <- file.path(root, "outputs", "HmscEcoEvo",
  "case05_v8_full200_terrain_resistance_20260926")
cache <- file.path(pkg_root, "derived_inputs",
  "case05_v2_paleomap_h3_carrier_325_0Ma_1deg",
  "movement_cache_arias_topographic_rate0p4")
for (package in c("HmscEcoEvo", "ggplot2", "terra")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Missing ", package, call. = FALSE)
}
if (utils::packageVersion("HmscEcoEvo") < package_version("1.0.4") ||
    sum(HmscEcoEvo::hee_core_process_catalog()$layer == "core_process") != 6L) {
  stop("Case07 maps require HmscEcoEvo >= 1.0.4 with the six-process API.",
       call. = FALSE)
}
if (!file.exists(file.path(out, "05_quality", "case07_quality_gates.csv"))) {
  stop("Run case07_patch_evidence_ladder.R first.", call. = FALSE)
}
prior_quality <- utils::read.csv(file.path(out, "05_quality",
  "case07_quality_gates.csv"), stringsAsFactors = FALSE)
if (!all(prior_quality$pass)) stop("Case07 base results have failed quality gates.")
source_quality <- utils::read.csv(file.path(corrected,
  "case05_patch_quality_gates.csv"), stringsAsFactors = FALSE)
if (!all(source_quality$pass)) stop("Corrected Case05 source has failed quality gates.")

map_dir <- file.path(out, "06_real_time_maps")
dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)
field <- readRDS(file.path(terrain, "05_shared_scale_maps",
  "case05_v8_map_values.rds"))
times <- as.numeric(names(field$lineage_location_density))
if (length(times) != 66L || !identical(times, sort(times, decreasing = TRUE))) {
  stop("Expected all 66 native Case05 time slices in descending Ma order.")
}
selected <- as.numeric(args$max_times %||% length(times))
if (!is.finite(selected) || selected < 1 || selected > length(times)) {
  stop("max_times must be between 1 and 66.")
}
times <- head(times, as.integer(selected))
focal <- c("Ilex_aquifolium", "Camellia_sinensis", "Jasminum_officinale")
lookup <- utils::read.csv(file.path(corrected,
  "four_scheme_location_mass_in_supported_area.csv"), stringsAsFactors = FALSE)
lookup <- lookup[lookup$scenario == "terrain_resistance" &
  lookup$focal_tip %in% focal, , drop = FALSE]
if (anyDuplicated(lookup[, c("focal_tip", "time_ma")])) {
  stop("Duplicate focal-tip/time lineage lookup.")
}
members <- setNames(lapply(focal, function(tip) utils::read.csv(file.path(corrected,
  paste0(tip, "_patch_cells.csv")), stringsAsFactors = FALSE)), focal)
candidates <- setNames(lapply(focal, function(tip) utils::read.csv(file.path(corrected,
  paste0(tip, "_candidate_refugia.csv")), stringsAsFactors = FALSE)), focal)
cache_index <- utils::read.csv(file.path(cache,
  "case05_v2_movement_cache_index.csv"), stringsAsFactors = FALSE)

grid <- data.frame(cell_id = sprintf("g4_%04d", seq_len(4050L)),
  lon = rep(seq(-178, 178, by = 4), times = 45),
  lat = rep(seq(-88, 88, by = 4), each = 90), stringsAsFactors = FALSE)
lat_s <- grid$lat - 2; lat_n <- grid$lat + 2
grid$cell_area_km2 <- 6371.0088^2 * (4 * pi / 180) *
  (sin(lat_n * pi / 180) - sin(lat_s * pi / 180))
bin4 <- function(lon, lat) {
  ix <- pmin(89L, pmax(0L, floor((as.numeric(lon) + 180) / 4)))
  iy <- pmin(44L, pmax(0L, floor((as.numeric(lat) + 90) / 4)))
  as.integer(iy * 90L + ix + 1L)
}
base_cells <- data.frame(cell_id = grid$cell_id, time_ma = 0,
  lon = grid$lon, lat = grid$lat, elevation_m = 0,
  topographic_resistance = 0, movement_allowed = 1)
base_edges <- HmscEcoEvo::hee_dispersal_connectivity_graph(
  base_cells, resolution_deg = 4, neighbours = 8, omega_topo = 0
)$edges[, c("from_cell_id", "to_cell_id")]
offsets <- expand.grid(dx = -3:3, dy = -3:3)
offsets <- offsets[pmax(abs(offsets$dx), abs(offsets$dy)) >= 2 &
  pmax(abs(offsets$dx), abs(offsets$dy)) <= 3, , drop = FALSE]
ix <- rep(0:89, times = 45); iy <- rep(0:44, each = 90)
jump <- lapply(seq_len(nrow(offsets)), function(k) {
  yy <- iy + offsets$dy[[k]]
  good <- yy >= 0 & yy < 45
  from <- which(good)
  to <- yy[good] * 90 + ((ix[good] + offsets$dx[[k]]) %% 90) + 1L
  data.frame(from_cell_id = grid$cell_id[from],
             to_cell_id = grid$cell_id[to])
})
base_edges <- rbind(base_edges, do.call(rbind, jump))
edge_from <- match(base_edges$from_cell_id, grid$cell_id)
edge_to <- match(base_edges$to_cell_id, grid$cell_id)
edge_km <- HmscEcoEvo::hee_great_circle_distance_km(
  grid$lon[edge_from], grid$lat[edge_from],
  grid$lon[edge_to], grid$lat[edge_to])
if (anyNA(edge_from) || anyNA(edge_to) || any(!is.finite(edge_km))) {
  stop("Movement edge grid has invalid endpoints or distances.")
}

carrier_at <- function(time_ma) {
  path <- cache_index$file[match(time_ma, cache_index$time_ma)]
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("Missing real ancient Earth carrier at ", time_ma, " Ma.")
  }
  readRDS(path)$carrier
}
permeability_at <- function(carrier) {
  active <- !is.na(carrier$active_land) & carrier$active_land &
    is.finite(carrier$land_area_km2) & carrier$land_area_km2 > 0 &
    is.finite(carrier$paleo_lon) & is.finite(carrier$paleo_lat) &
    is.finite(carrier$topographic_permeability)
  z <- carrier[active, , drop = FALSE]
  b <- bin4(z$paleo_lon, z$paleo_lat)
  numerator <- rowsum(z$land_area_km2 * z$topographic_permeability,
    b, reorder = FALSE)
  denominator <- rowsum(z$land_area_km2, b, reorder = FALSE)
  value <- numeric(nrow(grid))
  index <- as.integer(rownames(numerator))
  value[index] <- numerator[, 1L] /
    denominator[match(rownames(numerator), rownames(denominator)), 1L]
  pmin(1, pmax(0, value))
}

if (!index_only) {
diagnostics <- list(); rows <- list()
fields <- setNames(lapply(focal, function(.) list()), focal)
max_value <- 0
for (k in seq_along(times)) {
  time_ma <- times[[k]]
  key <- as.character(time_ma)
  land_fraction <- field$landscape_weight[[match(key,
    names(field$lineage_location_density))]]
  if (length(land_fraction) != nrow(grid) || anyNA(land_fraction)) {
    stop("Invalid native-time land mask at ", time_ma, " Ma.")
  }
  land <- land_fraction > 0
  permeability <- permeability_at(carrier_at(time_ma))
  edge_cost <- edge_km * (1 + 0.5 *
    (1 - pmin(1, pmax(0,
      (permeability[edge_from] + permeability[edge_to]) / 2))))
  edges <- base_edges
  edges$time_ma <- time_ma
  edges$distance_km <- edge_km
  edges$effective_cost_km <- edge_cost
  landscape_cells <- data.frame(cell_id = grid$cell_id, time_ma = time_ma,
    lon = grid$lon, lat = grid$lat, movement_weight = land_fraction,
    cell_area_km2 = grid$cell_area_km2)
  landscape <- HmscEcoEvo::hee_dispersal_landscape_weights(
    landscape_cells, movement_weight_col = "movement_weight",
    cell_area_col = "cell_area_km2")
  kernel <- HmscEcoEvo::hee_dispersal_spherical_kernel(
    cells = landscape, edges = edges, delta_t_myr = 1,
    diffusion_variance_rad2_per_myr = 0.0002,
    source_emigration_rate_per_myr = 0.4, time_ma = time_ma,
    edge_cost_col = "effective_cost_km", spatial_domain = "global_explicit")
  move <- kernel$edges
  from <- match(move$from_cell_id, grid$cell_id)
  to <- match(move$to_cell_id, grid$cell_id)
  if (anyNA(from) || anyNA(to) || any(!land[from] | !land[to])) {
    stop("Snapshot movement kernel includes missing or marine endpoints at ",
         time_ma, " Ma.")
  }
  move_probability <- 1 - exp(-0.4)
  for (tip in focal) {
    record <- lookup[lookup$focal_tip == tip & lookup$time_ma == time_ma, , drop = FALSE]
    if (nrow(record) != 1L) stop("Missing active ancestor lookup for ", tip,
                                " at ", time_ma, " Ma.")
    lineage <- record$active_ancestral_lineage[[1L]]
    densities <- field$lineage_location_density[[key]]
    if (!lineage %in% colnames(densities)) {
      stop("Case05 map lacks active ancestor ", lineage, " at ", time_ma, " Ma.")
    }
    p <- as.numeric(densities[, lineage])
    if (length(p) != nrow(grid) || any(!is.finite(p) | p < 0) ||
        abs(sum(p) - 1) > 1e-7 || sum(p[!land]) > 1e-7) {
      stop("Invalid land-only Case05 location density for ", tip,
           " at ", time_ma, " Ma.")
    }
    posterior_edge <- p[from] * move_probability * move$directional_weight
    if (any(!is.finite(posterior_edge) | posterior_edge < 0)) {
      stop("Invalid forward-predictive edge mass.")
    }
    origin <- numeric(nrow(grid)); target <- origin
    if (length(posterior_edge)) {
      from_sum <- rowsum(posterior_edge, from, reorder = FALSE)
      to_sum <- rowsum(posterior_edge, to, reorder = FALSE)
      origin[as.integer(rownames(from_sum))] <- from_sum[, 1L]
      target[as.integer(rownames(to_sum))] <- to_sum[, 1L]
    }
    if (abs(sum(origin) - sum(target)) > 1e-10) {
      stop("Movement source/target accounting mismatch.")
    }
    max_value <- max(max_value, origin, target)
    member <- members[[tip]]
    candidate <- candidates[[tip]]
    selected_patch <- candidate$to_patch_id[candidate$time_ma == time_ma]
    patch_cell <- member[member$time_ma == time_ma &
      member$patch_id %in% selected_patch, c("cell_id", "patch_id"), drop = FALSE]
    refuge <- integer(nrow(grid))
    if (length(selected_patch)) {
      high <- candidate$to_patch_id[candidate$time_ma == time_ma &
        candidate$geographic_evidence == "location_mass_above_threshold"]
      idx <- match(patch_cell$cell_id, grid$cell_id)
      if (anyNA(idx)) stop("Candidate refugium contains a missing grid cell.")
      refuge[idx] <- ifelse(patch_cell$patch_id %in% high, 2L, 1L)
    }
    if (time_ma == max(as.numeric(names(field$lineage_location_density))) ||
        time_ma == min(as.numeric(names(field$lineage_location_density)))) {
      refuge[] <- NA_integer_
      map_status <- "not_evaluable_boundary_slice"
    } else map_status <- "evaluated"
    if (any(!is.na(refuge[!land]) & refuge[!land] > 0) ||
        any(origin[!land] > 0) || any(target[!land] > 0)) {
      stop("Map has biological value in marine cells.")
    }
    rank <- order(posterior_edge, decreasing = TRUE)
    rank <- head(rank[posterior_edge[rank] > 0 &
      abs(grid$lon[from[rank]] - grid$lon[to[rank]]) <= 180], 70L)
    top_edges <- data.frame(lon = grid$lon[from[rank]], lat = grid$lat[from[rank]],
      lon_to = grid$lon[to[rank]], lat_to = grid$lat[to[rank]],
      mass = posterior_edge[rank])
    fields[[tip]][[key]] <- list(land = land, refuge = refuge,
      origin = origin, target = target, top_edges = top_edges,
      lineage = lineage, map_status = map_status)
    diagnostics[[length(diagnostics) + 1L]] <- data.frame(focal_tip = tip,
      time_ma = time_ma, active_ancestor = lineage,
      location_mass_on_land = sum(p[land]),
      candidate_patches = length(unique(selected_patch)),
      candidate_land_cells = sum(refuge > 0, na.rm = TRUE),
      expected_active_movement_mass_per_myr = sum(origin),
      source_target_error = abs(sum(origin) - sum(target)),
      marine_source_mass = sum(origin[!land]),
      marine_target_mass = sum(target[!land]),
      source_model = "Case05_v8_terrain_location_marginal",
      movement_model = "native_snapshot_terrain_spherical_kernel_1Myr",
      edge_interpretation = "forward_predictive_not_joint_posterior",
      refuge_status = map_status)
  }
  message("Case07 real map inputs: ", k, "/", length(times),
          " native slices; ", time_ma, " Ma")
}

if (!is.finite(max_value) || max_value <= 0) {
  stop("All forward-predictive movement fields are zero.")
}
diagnostic_table <- do.call(rbind, diagnostics)
write.csv(diagnostic_table, file.path(map_dir, "case07_real_time_diagnostics.csv"),
  row.names = FALSE, na = "")

save_geotiff <- function(values, path) {
  raster <- terra::rast(ncols = 90, nrows = 45, xmin = -180,
    xmax = 180, ymin = -90, ymax = 90, crs = "EPSG:4326")
  cell <- terra::cellFromXY(raster, as.matrix(grid[, c("lon", "lat")]))
  full <- rep(NA_real_, terra::ncell(raster))
  full[cell] <- values
  terra::values(raster) <- full
  terra::writeRaster(raster, path, overwrite = TRUE,
    gdal = c("COMPRESS=DEFLATE"))
}
base_theme <- function() ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold", colour = "#173552"),
    plot.caption = ggplot2::element_text(colour = "#4f6070"))
map_rows <- list()
for (tip in focal) {
  for (time_ma in times) {
    key <- as.character(time_ma)
    z <- fields[[tip]][[key]]
    current <- grid[, c("lon", "lat")]
    for (metric in c("candidate_refugia", "potential_movement_origins",
                     "potential_corridor_pressure")) {
      folder <- file.path(map_dir, tip, metric)
      dir.create(folder, recursive = TRUE, showWarnings = FALSE)
      filename <- sprintf("%s_%03dMa", metric, time_ma)
      png <- file.path(folder, paste0(filename, ".png"))
      tif <- file.path(folder, paste0(filename, ".tif"))
      current$value <- switch(metric,
        candidate_refugia = as.numeric(z$refuge),
        potential_movement_origins = z$origin,
        potential_corridor_pressure = z$target)
      current$value[!z$land] <- NA_real_
      if (metric == "candidate_refugia") {
        current$class <- factor(current$value, levels = 0:2,
          labels = c("Other land", "Candidate", "Candidate + location evidence"))
        plot <- ggplot2::ggplot(current, ggplot2::aes(lon, lat)) +
          ggplot2::geom_tile(data = current[z$land, , drop = FALSE],
                            fill = "#e0e3e2", colour = NA) +
          ggplot2::geom_tile(ggplot2::aes(fill = class), colour = NA) +
          ggplot2::scale_fill_manual(values = c("Other land" = "#e0e3e2",
            Candidate = "#2f9174", "Candidate + location evidence" = "#e2a13d"),
            drop = FALSE, na.translate = FALSE, na.value = NA,
            name = "Evidence class") +
          ggplot2::labs(title = paste(tip, "| ecological refugium candidates |",
                                      time_ma, "Ma"),
            subtitle = if (z$map_status == "evaluated")
              "Contraction-surviving supported patches; real corrected Case05 inputs"
              else "Boundary time: strict incoming-and-next-link test is not evaluable",
            caption = paste("Candidate environmental continuity, not verified population survival.",
              "Ocean is white; boundary-time land is NA."))
        evidence <- "corrected_plant200_environmental_opportunity"
        unit <- "class_0_1_2"
      } else {
        current$display <- current$value
        current$display[!is.na(current$display) & current$display <= 0] <- NA_real_
        plot <- ggplot2::ggplot(current, ggplot2::aes(lon, lat)) +
          ggplot2::geom_tile(data = current[z$land, , drop = FALSE],
                            fill = "#e1e4e2", colour = NA) +
          ggplot2::geom_tile(ggplot2::aes(fill = display), colour = NA) +
          ggplot2::scale_fill_viridis_c(option = "C", trans = "log10",
            limits = c(1e-10, max_value), oob = scales::squish,
            na.value = NA, name = "Predictive mass\nper 1 Myr")
        if (metric == "potential_corridor_pressure" &&
            nrow(z$top_edges)) {
          plot <- plot + ggplot2::geom_segment(data = z$top_edges,
            ggplot2::aes(x = lon, y = lat, xend = lon_to, yend = lat_to),
            inherit.aes = FALSE, colour = "#9b3d35", alpha = .55,
            linewidth = .35, arrow = grid::arrow(length = grid::unit(.07, "cm")))
        }
        plot <- plot + ggplot2::labs(
          title = paste(tip, "|", if (metric == "potential_movement_origins")
            "potential movement origins" else "potential corridor pressure |",
            time_ma, "Ma"),
          subtitle = paste("Case05 location marginal + 1 Myr terrain movement kernel;",
                           "active ancestor", z$lineage),
          caption = paste("Forward-predictive scenario, NOT a posterior-used corridor or",
            "demographic source. Grey land has zero displayed mass; ocean is white."))
        evidence <- "real_paleogeography_forward_predictive_scenario"
        unit <- "expected_location_transition_mass_per_1Myr"
      }
      plot <- plot + ggplot2::coord_equal(expand = FALSE) +
        ggplot2::labs(x = "Palaeolongitude", y = "Palaeolatitude") + base_theme()
      ggplot2::ggsave(png, plot, width = 11.5, height = 6.1, dpi = 120)
      save_geotiff(current$value, tif)
      map_rows[[length(map_rows) + 1L]] <- data.frame(
        focal_tip = tip, time_ma = time_ma, active_ancestor = z$lineage,
        metric = metric, evidence_class = evidence, status = z$map_status,
        unit = unit, scale_min = if (metric == "candidate_refugia") 0 else 1e-10,
        scale_max = if (metric == "candidate_refugia") 2 else max_value,
        png = normalizePath(png, winslash = "/", mustWork = TRUE),
        geotiff = normalizePath(tif, winslash = "/", mustWork = TRUE))
    }
  }
  message("Case07 maps finished for ", tip)
}
map_index <- do.call(rbind, map_rows)
write.csv(map_index, file.path(map_dir, "case07_real_time_map_index.csv"),
  row.names = FALSE, na = "")
quality <- data.frame(check = c("all_native_times", "all_focal_paths",
  "all_three_map_types", "all_png_tif_files", "same_scale_within_metric",
  "source_target_mass_conserved", "marine_mass_zero",
  "real_corridor_not_mislabeled_posterior"),
  pass = c(length(times) == 66L,
    setequal(unique(map_index$focal_tip), focal),
    setequal(unique(map_index$metric), c("candidate_refugia",
      "potential_movement_origins", "potential_corridor_pressure")),
    all(file.exists(map_index$png)) && all(file.exists(map_index$geotiff)),
    length(unique(map_index$scale_max[map_index$metric != "candidate_refugia"])) == 1L,
    all(diagnostic_table$source_target_error < 1e-10),
    all(diagnostic_table$marine_source_mass == 0 &
        diagnostic_table$marine_target_mass == 0),
    all(map_index$evidence_class[map_index$metric != "candidate_refugia"] ==
      "real_paleogeography_forward_predictive_scenario")),
  stringsAsFactors = FALSE)
write.csv(quality, file.path(map_dir, "case07_real_time_quality_gates.csv"),
  row.names = FALSE)
if (!all(quality$pass)) {
  stop("Case07 native-time map atlas failed: ",
    paste(quality$check[!quality$pass], collapse = ", "), call. = FALSE)
}
} else {
  map_index <- utils::read.csv(file.path(map_dir,
    "case07_real_time_map_index.csv"), stringsAsFactors = FALSE)
  quality <- utils::read.csv(file.path(map_dir,
    "case07_real_time_quality_gates.csv"), stringsAsFactors = FALSE)
  if (!all(quality$pass) || !all(file.exists(map_index$png)) ||
      !all(file.exists(map_index$geotiff))) {
    stop("Existing atlas failed quality gates; rerun the complete script.")
  }
}

relative <- function(path) substring(path, nchar(map_dir) + 2L)
entries <- vapply(seq_len(nrow(map_index)), function(i) {
  x <- map_index[i, ]
  sprintf("{tip:'%s',time:%s,metric:'%s',path:'%s',status:'%s'}",
    x$focal_tip, x$time_ma, x$metric,
    gsub("\\\\", "/", relative(x$png)), x$status)
}, character(1))
html <- c('<!doctype html><html lang="zh"><head><meta charset="utf-8">',
  '<meta name="viewport" content="width=device-width,initial-scale=1">',
  '<title>Case07 native-time maps</title><style>',
  'body{font:16px system-ui,Arial,sans-serif;color:#173552;margin:0;background:#f5f7f6}',
  'header{background:#fff;border-bottom:1px solid #cbd5d2;padding:18px 24px}',
  'h1{font-size:24px;margin:0 0 6px}p{margin:6px 0;color:#435d66}',
  'main{max-width:1280px;margin:auto;padding:18px 24px}',
  '.controls{display:flex;gap:14px;flex-wrap:wrap;margin-bottom:14px}',
  'label{display:flex;flex-direction:column;gap:5px;font-weight:600}',
  'select{font:inherit;padding:8px;border:1px solid #92ada9;background:white;min-width:180px}',
  '.frame{background:white;border:1px solid #cbd5d2;padding:12px}',
  'img{width:100%;height:auto;display:block}#status{font-weight:600;margin:10px 0}',
  '</style></head><body><header><h1>Case07 | 逐期真实古地理地图</h1>',
  '<p>候选避难所来自更正的 Plant200 结果；移动起点和通道是前向预测情景，不是后验已使用路径。</p>',
  '</header><main><div class="controls">',
  '<label>祖先路径<select id="tip"><option>Ilex_aquifolium</option><option>Camellia_sinensis</option><option>Jasminum_officinale</option></select></label>',
  '<label>指标<select id="metric"><option value="candidate_refugia">候选生态避难所</option><option value="potential_movement_origins">潜在移动起点</option><option value="potential_corridor_pressure">潜在通道压力</option></select></label>',
  paste0('<label>时间 (Ma)<select id="time">',
    paste(sprintf('<option value="%s">%s Ma</option>', times, times), collapse = ''),
    '</select></label></div><p id="status"></p><div class="frame">'),
  '<img id="map" alt="Selected palaeogeographic map"></div>',
  '<p>每张地图的 PNG 与 GeoTIFF 路径在 case07_real_time_map_index.csv；所有时间使用相同指标色标。边界时间的严格避难所测试为 NA。</p>',
  '<p><a href="../07_dispersal_source_corridor_maps/index.html">Plant200 全体活动谱系的源区/通道地图</a> · <a href="../Case07_patch_evidence_ladder.html">阅读六函数完整教程</a></p>',
  '<script>const data=[', paste(entries, collapse = ',\n'), '];',
  'const tip=document.getElementById("tip"),metric=document.getElementById("metric"),time=document.getElementById("time");',
  'function update(){const row=data.find(x=>x.tip===tip.value&&x.metric===metric.value&&x.time===Number(time.value));',
  'document.getElementById("map").src=row.path;',
  'document.getElementById("status").textContent=`${tip.value} · ${time.value} Ma · ${metric.options[metric.selectedIndex].text} · ${row.status}`;}',
  '[tip,metric,time].forEach(x=>x.addEventListener("change",update));update();</script></main></body></html>')
writeLines(html, file.path(map_dir, "index.html"), useBytes = TRUE)
message("Case07 real-time atlas complete: ", nrow(map_index),
        " PNG + GeoTIFF pairs; ", length(times), " times; ", map_dir)
