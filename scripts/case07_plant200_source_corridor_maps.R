#!/usr/bin/env Rscript

# Forward-predictive local movement maps from real Case05 Plant200 marginals.
# These are not posterior-used corridors or demographic source populations.

`%||%` <- function(x, y) if (is.null(x)) y else x
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (!length(hit)) default else substring(hit[[1L]], nchar(name) + 4L)
}
script <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
if (!length(script)) stop("Run with Rscript.", call. = FALSE)
pkg_root <- normalizePath(file.path(dirname(script[[1L]]), ".."), winslash = "/")
root <- normalizePath(file.path(pkg_root, ".."), winslash = "/")
base_out <- file.path(root, "outputs", "HmscEcoEvo",
  "case07_patch_evidence_ladder_20260926")
out <- normalizePath(get_arg("output", file.path(base_out,
  "07_dispersal_source_corridor_maps")), winslash = "/", mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (package in c("HmscEcoEvo", "ggplot2", "terra", "scales")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Missing R package: ", package, call. = FALSE)
  }
}
if (utils::packageVersion("HmscEcoEvo") < package_version("1.0.4") ||
    sum(HmscEcoEvo::hee_core_process_catalog()$layer == "core_process") != 6L) {
  stop("Case07 maps require HmscEcoEvo >= 1.0.4 with the six-process API.",
       call. = FALSE)
}
case05 <- file.path(root, "outputs", "HmscEcoEvo",
  "case05_v8_full200_terrain_resistance_20260926")
cache <- file.path(pkg_root, "derived_inputs",
  "case05_v2_paleomap_h3_carrier_325_0Ma_1deg",
  "movement_cache_arias_topographic_rate0p4")
field <- readRDS(file.path(case05, "05_shared_scale_maps",
  "case05_v8_map_values.rds"))
times <- as.numeric(names(field$lineage_location_density))
if (length(times) != 66L || !identical(times, sort(times, decreasing = TRUE))) {
  stop("Expected the 66 native Case05 time slices (325 to 0 Ma).")
}
index <- utils::read.csv(file.path(cache,
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
edge <- HmscEcoEvo::hee_dispersal_connectivity_graph(base_cells,
  resolution_deg = 4, neighbours = 8, omega_topo = 0
)$edges[, c("from_cell_id", "to_cell_id")]
from0 <- match(edge$from_cell_id, grid$cell_id)
to0 <- match(edge$to_cell_id, grid$cell_id)
if (anyNA(from0) || anyNA(to0)) stop("Invalid 4-degree graph endpoints.")
ix <- rep(0:89, times = 45); iy <- rep(0:44, each = 90)
dx <- pmin(abs(ix[from0] - ix[to0]), 90 - abs(ix[from0] - ix[to0]))
dy <- abs(iy[from0] - iy[to0])
if (any(dx > 1 | dy > 1)) stop("Local corridor graph has nonlocal edges.")
diagonal <- dx == 1 & dy == 1
via_a <- iy[from0] * 90L + ix[to0] + 1L
via_b <- iy[to0] * 90L + ix[from0] + 1L
distance <- HmscEcoEvo::hee_great_circle_distance_km(
  grid$lon[from0], grid$lat[from0], grid$lon[to0], grid$lat[to0])

get_carrier <- function(time_ma) {
  path <- index$file[match(time_ma, index$time_ma)]
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("Missing palaeogeographic carrier for ", time_ma, " Ma.")
  }
  readRDS(path)$carrier
}
get_permeability <- function(carrier) {
  valid <- !is.na(carrier$active_land) & carrier$active_land &
    is.finite(carrier$land_area_km2) & carrier$land_area_km2 > 0 &
    is.finite(carrier$paleo_lon) & is.finite(carrier$paleo_lat) &
    is.finite(carrier$topographic_permeability)
  z <- carrier[valid, , drop = FALSE]
  b <- bin4(z$paleo_lon, z$paleo_lat)
  num <- rowsum(z$land_area_km2 * z$topographic_permeability,
    b, reorder = FALSE)
  den <- rowsum(z$land_area_km2, b, reorder = FALSE)
  result <- numeric(nrow(grid)); at <- as.integer(rownames(num))
  result[at] <- num[, 1L] / den[match(rownames(num), rownames(den)), 1L]
  pmin(1, pmax(0, result))
}
sum_by_cell <- function(value, index) {
  result <- numeric(nrow(grid))
  if (length(value)) {
    sums <- rowsum(value, index, reorder = FALSE)
    result[as.integer(rownames(sums))] <- sums[, 1L]
  }
  result
}
save_tif <- function(value, path) {
  raster <- terra::rast(ncols = 90, nrows = 45, xmin = -180,
    xmax = 180, ymin = -90, ymax = 90, crs = "EPSG:4326")
  cell <- terra::cellFromXY(raster, as.matrix(grid[, c("lon", "lat")]))
  values <- rep(NA_real_, terra::ncell(raster)); values[cell] <- value
  terra::values(raster) <- values
  terra::writeRaster(raster, path, overwrite = TRUE, gdal = "COMPRESS=DEFLATE")
}

fields <- vector("list", length(times)); names(fields) <- as.character(times)
diagnostics <- vector("list", length(times))
max_source <- 0; max_corridor <- 0
edge_dir <- file.path(out, "edge_flux_tables")
dir.create(edge_dir, recursive = TRUE, showWarnings = FALSE)
for (k in seq_along(times)) {
  t <- times[[k]]; key <- as.character(t)
  land_fraction <- field$landscape_weight[[k]]
  land <- land_fraction > 0
  if (length(land) != nrow(grid) || anyNA(land_fraction)) {
    stop("Invalid land mask at ", t, " Ma.")
  }
  density <- field$lineage_location_density[[key]]
  if (nrow(density) != nrow(grid) || !ncol(density) ||
      any(!is.finite(density) | density < 0) ||
      any(abs(colSums(density) - 1) > 1e-7) ||
      any(density[!land, , drop = FALSE] > 1e-7)) {
    stop("Invalid active-lineage location marginals at ", t, " Ma.")
  }
  p <- rowMeans(density)
  keep <- land[from0] & land[to0] &
    (!diagonal | land[via_a] | land[via_b])
  current <- edge[keep, , drop = FALSE]
  from <- from0[keep]; to <- to0[keep]
  permeability <- get_permeability(get_carrier(t))
  current$time_ma <- t
  current$distance_km <- distance[keep]
  current$effective_cost_km <- distance[keep] *
    (1 + 0.5 * (1 - (permeability[from] + permeability[to]) / 2))
  landscape_cells <- data.frame(cell_id = grid$cell_id, time_ma = t,
    lon = grid$lon, lat = grid$lat, movement_weight = land_fraction,
    cell_area_km2 = grid$cell_area_km2)
  landscape <- HmscEcoEvo::hee_dispersal_landscape_weights(
    landscape_cells, movement_weight_col = "movement_weight",
    cell_area_col = "cell_area_km2")
  kernel <- HmscEcoEvo::hee_dispersal_spherical_kernel(
    cells = landscape, edges = current, delta_t_myr = 1,
    diffusion_variance_rad2_per_myr = 0.0002,
    source_emigration_rate_per_myr = 0.4, time_ma = t,
    edge_cost_col = "effective_cost_km", spatial_domain = "global_explicit")
  move <- kernel$edges
  from <- match(move$from_cell_id, grid$cell_id)
  to <- match(move$to_cell_id, grid$cell_id)
  if (anyNA(from) || anyNA(to) || any(!land[from] | !land[to])) {
    stop("Marine or missing endpoint in local movement kernel at ", t, " Ma.")
  }
  flux <- p[from] * (1 - exp(-0.4)) * move$directional_weight
  if (any(!is.finite(flux) | flux < 0)) stop("Invalid predicted edge flux.")
  origin <- sum_by_cell(flux, from)
  target <- sum_by_cell(flux, to)
  corridor <- (origin + target) / 2
  if (abs(sum(origin) - sum(target)) > 1e-10 ||
      any(origin[!land] > 0) || any(corridor[!land] > 0)) {
    stop("Movement mass or marine-mask check failed at ", t, " Ma.")
  }
  positive <- origin[land & origin > 0]
  threshold <- if (length(positive))
    unname(stats::quantile(positive, 0.9, names = FALSE)) else Inf
  source_area <- land & origin >= threshold & origin > 0
  pair <- paste(pmin(from, to), pmax(from, to), sep = "_")
  pair_sum <- rowsum(flux, pair, reorder = FALSE)
  pair_first <- match(rownames(pair_sum), pair)
  rank <- head(order(pair_sum[, 1L], decreasing = TRUE), 80L)
  pair_from <- pmin(from[pair_first[rank]], to[pair_first[rank]])
  pair_to <- pmax(from[pair_first[rank]], to[pair_first[rank]])
  display <- abs(grid$lon[pair_from] - grid$lon[pair_to]) <= 180 &
    pair_sum[rank, 1L] > 0
  top <- data.frame(lon = grid$lon[pair_from[display]],
    lat = grid$lat[pair_from[display]],
    lon_to = grid$lon[pair_to[display]],
    lat_to = grid$lat[pair_to[display]],
    flux = pair_sum[rank[display], 1L])
  edge_table <- data.frame(time_ma = t, from_cell_id = move$from_cell_id,
    to_cell_id = move$to_cell_id, predicted_mass_per_1Myr = flux,
    distance_km = move$distance_km,
    effective_cost_km = move$effective_movement_distance_km)
  edge_file <- file.path(edge_dir, sprintf("predicted_local_edges_%03dMa.csv", t))
  utils::write.csv(edge_table, edge_file, row.names = FALSE)
  max_source <- max(max_source, origin)
  max_corridor <- max(max_corridor, corridor)
  fields[[key]] <- list(land = land, origin = origin,
    corridor = corridor, source_area = source_area, top_edges = top,
    active_lineages = ncol(density), edge_file = edge_file)
  diagnostics[[k]] <- data.frame(time_ma = t, active_lineages = ncol(density),
    n_land_cells = sum(land), n_candidate_source_cells = sum(source_area),
    n_active_edges = nrow(edge_table), n_drawn_corridor_links = nrow(top),
    expected_movement_mass_per_mean_lineage = sum(origin),
    conservation_error = abs(sum(origin) - sum(target)),
    marine_source_mass = sum(origin[!land]),
    marine_corridor_mass = sum(corridor[!land]),
    corner_cut_edges = sum(diagonal[keep] &
      !land[via_a[keep]] & !land[via_b[keep]]),
    stringsAsFactors = FALSE)
  message("Plant200 source/corridor inputs ", k, "/", length(times),
    ": ", t, " Ma; ", ncol(density), " active lineages")
}
if (max_source <= 0 || max_corridor <= 0) stop("Predicted movement is zero.")
diagnostics <- do.call(rbind, diagnostics)
utils::write.csv(diagnostics, file.path(out, "plant200_movement_diagnostics.csv"),
  row.names = FALSE)

map_rows <- list()
theme_map <- ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold", colour = "#173552"),
    plot.caption = ggplot2::element_text(colour = "#4f6070"))
for (t in times) {
  z <- fields[[as.character(t)]]
  for (metric in c("potential_dispersal_source_areas",
                   "potential_land_corridors")) {
    folder <- file.path(out, metric)
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    base <- sprintf("%s_%03dMa", metric, t)
    png <- file.path(folder, paste0(base, ".png"))
    tif <- file.path(folder, paste0(base, ".tif"))
    value <- if (metric == "potential_dispersal_source_areas")
      z$origin else z$corridor
    value[!z$land] <- NA_real_
    frame <- grid[, c("lon", "lat")]
    frame$value <- value
    frame$display <- value
    frame$display[!is.na(frame$display) & frame$display <= 0] <- NA_real_
    plot <- ggplot2::ggplot(frame, ggplot2::aes(lon, lat)) +
      ggplot2::geom_tile(data = frame[z$land, , drop = FALSE],
        fill = "#dfe4e1", colour = NA) +
      ggplot2::geom_tile(ggplot2::aes(fill = display), colour = NA) +
      ggplot2::scale_fill_viridis_c(option = if (metric ==
        "potential_dispersal_source_areas") "C" else "D",
        trans = "log10", limits = c(1e-10, if (metric ==
        "potential_dispersal_source_areas") max_source else max_corridor),
        oob = scales::squish, na.value = NA,
        name = "Predicted mass\nper 1 Myr")
    if (metric == "potential_dispersal_source_areas") {
      top <- frame[z$source_area, , drop = FALSE]
      plot <- plot + ggplot2::geom_tile(data = top,
        ggplot2::aes(x = lon, y = lat), inherit.aes = FALSE,
        fill = NA, colour = "#174f47", linewidth = 0.35) +
        ggplot2::labs(title = paste("Plant200 | potential dispersal source areas |",
          t, "Ma"),
          subtitle = paste(z$active_lineages, "active sampled-surviving lineages;",
            "outlined = upper 10% positive land exporters within this time"),
          caption = paste("Forward-predictive 1 Myr local emigration, averaged over active lineages.",
            "Not demographic source populations; ocean is blank."))
    } else {
      plot <- plot + ggplot2::geom_segment(data = z$top_edges,
        ggplot2::aes(x = lon, y = lat, xend = lon_to, yend = lat_to),
        inherit.aes = FALSE, colour = "white", linewidth = 2.1,
        alpha = 0.9) +
        ggplot2::geom_segment(data = z$top_edges,
        ggplot2::aes(x = lon, y = lat, xend = lon_to, yend = lat_to,
          linewidth = flux), inherit.aes = FALSE,
        colour = "#9e2e20", alpha = 0.94) +
        ggplot2::scale_linewidth(range = c(0.7, 1.5), guide = "none") +
        ggplot2::labs(title = paste("Plant200 | potential land corridors |",
          t, "Ma"),
          subtitle = paste(z$active_lineages, "active sampled-surviving lineages;",
            "outlined red = strongest 80 bidirectional local links"),
          caption = paste("Raster: mean endpoint throughput on the land graph.",
            "Predicted links are not historical paths or posterior-used movement."))
    }
    plot <- plot + ggplot2::coord_equal(expand = FALSE) +
      ggplot2::labs(x = "Palaeolongitude", y = "Palaeolatitude") + theme_map
    ggplot2::ggsave(png, plot, width = 11.5, height = 6.1, dpi = 140)
    save_tif(value, tif)
    map_rows[[length(map_rows) + 1L]] <- data.frame(
      time_ma = t, metric = metric, active_lineages = z$active_lineages,
      evidence_class = "real_inputs_forward_predictive_not_posterior",
      scale_min = 1e-10,
      scale_max = if (metric == "potential_dispersal_source_areas")
        max_source else max_corridor,
      png = normalizePath(png, winslash = "/", mustWork = TRUE),
      geotiff = normalizePath(tif, winslash = "/", mustWork = TRUE),
      edge_table = normalizePath(z$edge_file, winslash = "/", mustWork = TRUE))
  }
  message("Plant200 source/corridor maps complete: ", t, " Ma")
}
map_index <- do.call(rbind, map_rows)
utils::write.csv(map_index, file.path(out, "plant200_source_corridor_map_index.csv"),
  row.names = FALSE)
quality <- data.frame(check = c("all_66_times", "both_map_types",
  "all_png_tif_edge_tables", "land_only_movement",
  "source_target_mass_conserved", "no_diagonal_ocean_corner_cut",
  "shared_scale_within_each_metric", "tree_tip_count_at_zero",
  "root_lineage_count", "predictions_not_labeled_posterior"),
  pass = c(length(unique(map_index$time_ma)) == 66L,
    setequal(unique(map_index$metric), c("potential_dispersal_source_areas",
      "potential_land_corridors")),
    all(file.exists(map_index$png)) && all(file.exists(map_index$geotiff)) &&
      all(file.exists(map_index$edge_table)),
    all(diagnostics$marine_source_mass == 0 &
      diagnostics$marine_corridor_mass == 0),
    all(diagnostics$conservation_error < 1e-10),
    all(diagnostics$corner_cut_edges == 0),
    all(vapply(split(map_index$scale_max, map_index$metric),
      function(x) length(unique(x)) == 1L, logical(1))),
    diagnostics$active_lineages[diagnostics$time_ma == 0] == 200L,
    diagnostics$active_lineages[diagnostics$time_ma == 325] == 1L,
    all(map_index$evidence_class ==
      "real_inputs_forward_predictive_not_posterior")),
  stringsAsFactors = FALSE)
utils::write.csv(quality, file.path(out, "plant200_source_corridor_quality.csv"),
  row.names = FALSE)
if (!all(quality$pass)) stop("Map quality failed: ",
  paste(quality$check[!quality$pass], collapse = ", "), call. = FALSE)

relative <- function(path) substring(path, nchar(out) + 2L)
entries <- vapply(seq_len(nrow(map_index)), function(i) {
  x <- map_index[i, ]
  sprintf("{time:%s,metric:'%s',path:'%s',lineages:%s}",
    x$time_ma, x$metric, relative(x$png), x$active_lineages)
}, character(1))
html <- c('<!doctype html><html lang="zh"><head><meta charset="utf-8">',
  '<meta name="viewport" content="width=device-width,initial-scale=1">',
  '<title>Plant200 source areas and land corridors</title><style>',
  'body{font:16px system-ui,Arial,sans-serif;color:#173552;margin:0;background:#f5f7f6}',
  'header{background:white;border-bottom:1px solid #cbd5d2;padding:18px 24px}',
  'h1{font-size:24px;margin:0 0 6px}p{margin:6px 0;color:#435d66}',
  'main{max-width:1280px;margin:auto;padding:18px 24px}',
  '.controls{display:flex;gap:14px;flex-wrap:wrap;margin-bottom:14px}',
  'label{display:flex;flex-direction:column;gap:5px;font-weight:600}',
  'select{font:inherit;padding:8px;border:1px solid #92ada9;background:white;min-width:185px}',
  '.frame{background:white;border:1px solid #cbd5d2;padding:12px}',
  'img{width:100%;height:auto;display:block}#status{font-weight:600;margin:10px 0}',
  '</style></head><body><header><h1>Plant200 逐期扩散源区与陆地通道</h1>',
  '<p>真实古陆地与当期活动谱系的 Case05 位置边际 + 一步局地移动核。',
  '这是前向预测情景，不是后验已使用路径或种群人口统计学源地。</p>',
  '</header><main><div class="controls">',
  '<label>地图<select id="metric"><option value="potential_dispersal_source_areas">',
  '潜在扩散源区</option><option value="potential_land_corridors">潜在陆地扩散通道</option>',
  '</select></label>',
  paste0('<label>时间<select id="time">',
    paste(sprintf('<option value="%s">%s Ma</option>', times, times),
      collapse = ''), '</select></label></div>'),
  '<p id="status"></p><div class="frame"><img id="map" alt="Plant200 palaeogeographic movement map"></div>',
  '<p>深色描边的源区是当期正输出网格中的前 10%；通道图的红线是预测质量最大的 80 条局地陆地边。',
  '底图 4°，海洋留白；同一指标的全部时间使用同一色标。',
  '逐期 PNG、GeoTIFF 和完整有向边表列于地图索引。</p>',
  '<p><a href="../06_real_time_maps/index.html">三条祖先路径及候选避难所</a> · ',
  '<a href="../Case07_patch_evidence_ladder.html">案例教程</a></p>',
  '<script>const data=[', paste(entries, collapse = ',\n'), '];',
  'const metric=document.getElementById("metric"),time=document.getElementById("time");',
  'function update(){const row=data.find(x=>x.metric===metric.value&&x.time===Number(time.value));',
  'document.getElementById("map").src=row.path;',
  'document.getElementById("status").textContent=`${time.value} Ma · ${row.lineages} active lineages · ${metric.options[metric.selectedIndex].text}`;}',
  '[metric,time].forEach(x=>x.addEventListener("change",update));update();</script></main></body></html>')
writeLines(html, file.path(out, "index.html"), useBytes = TRUE)
message("Plant200 source/corridor atlas complete: ", nrow(map_index),
  " PNG and GeoTIFF pairs across ", length(times), " times; ", out)
