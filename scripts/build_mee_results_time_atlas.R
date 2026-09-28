#!/usr/bin/env Rscript

# Browser for all saved Case05/Case07 times. Re-render the archived support
# diversity field after excluding intensities below numerical precision.
stopifnot(all(vapply(c("ggplot2", "ragg", "terra", "scales", "jsonlite"),
                     requireNamespace, logical(1), quietly = TRUE)))

root <- normalizePath(file.path(getwd(), ".."), winslash = "/")
archive <- file.path(root, "outputs", "HmscEcoEvo")
out <- file.path(archive, "manuscript_results_six_process_20260927", "time_map_atlas")
div_dir <- file.path(out, "diversity_corrected")
dir.create(div_dir, recursive = TRUE, showWarnings = FALSE)
field_path <- file.path(archive, "case05_v8_full200_terrain_resistance_20260926",
                        "05_shared_scale_maps", "case05_v8_map_values.rds")
shared_path <- file.path(archive, "case05_v8_four_scheme_comparison_20260926",
                         "07_cross_scheme_shared_scale_maps",
                         "cross_scheme_shared_scale_map_index.csv")
movement_path <- file.path(archive, "case07_patch_evidence_ladder_20260926",
                           "07_dispersal_source_corridor_maps",
                           "plant200_source_corridor_map_index.csv")
refuge_path <- file.path(archive, "case07_patch_evidence_ladder_20260926",
                         "06_real_time_maps", "case07_real_time_map_index.csv")
stopifnot(all(file.exists(c(field_path, shared_path, movement_path, refuge_path))))

field <- readRDS(field_path)
shared <- utils::read.csv(shared_path, stringsAsFactors = FALSE)
movement <- utils::read.csv(movement_path, stringsAsFactors = FALSE)
refuges <- utils::read.csv(refuge_path, stringsAsFactors = FALSE)
grid <- expand.grid(lon = seq(-178, 178, 4), lat = seq(-88, 88, 4))
stopifnot(nrow(grid) == 4050L)
times <- sort(as.numeric(names(field$location_density)), decreasing = TRUE)
scale_rows <- shared[shared$metric == "07_support_effective", , drop = FALSE]
stopifnot(length(unique(scale_rows$scale_min)) == 1L,
          length(unique(scale_rows$scale_max)) == 1L,
          length(times) == 66L)
maximum <- unique(scale_rows$scale_max)
template <- terra::rast(ncols = 90, nrows = 45, xmin = -180, xmax = 180,
                        ymin = -90, ymax = 90, crs = "EPSG:4326")
cells <- terra::cellFromXY(template, as.matrix(grid[, c("lon", "lat")]))
stopifnot(!anyNA(cells), !anyDuplicated(cells))

div_index <- do.call(rbind, lapply(times, function(time_ma) {
  key <- as.character(time_ma)
  land_at <- match(key, names(field$location_density))
  land <- as.numeric(field$landscape_weight[[land_at]]) > 0
  intensity <- as.numeric(field$terminal_calibrated_lineage_support_intensity[[key]])
  effective <- as.numeric(field$support_mixture_effective_lineages[[key]])
  effective[!land | intensity <= .Machine$double.eps] <- NA_real_
  active <- ncol(field$lineage_location_density[[key]])
  if (max(effective, na.rm = TRUE) > active + 1e-9) {
    stop("Effective support exceeds active lineages at ", time_ma, " Ma")
  }

  tif <- file.path(div_dir, sprintf("effective_support_%sMa.tif", key))
  png <- file.path(div_dir, sprintf("effective_support_%sMa.png", key))
  raster <- template
  values <- rep(NA_real_, terra::ncell(raster))
  values[cells] <- effective
  terra::values(raster) <- values
  terra::writeRaster(raster, tif, overwrite = TRUE, datatype = "FLT4S")

  tab <- grid
  tab$land <- land
  tab$value <- effective
  p <- ggplot2::ggplot(tab, ggplot2::aes(lon, lat)) +
    ggplot2::geom_tile(data = tab[tab$land, , drop = FALSE],
                       fill = "#E7ECE9", width = 4, height = 4) +
    ggplot2::geom_tile(ggplot2::aes(fill = value), width = 4, height = 4) +
    ggplot2::scale_fill_viridis_c(option = "C", limits = c(1, maximum),
      oob = scales::squish, na.value = NA,
      name = "Effective support\nlineages") +
    ggplot2::coord_fixed(expand = FALSE) +
    ggplot2::scale_x_continuous(limits = c(-180, 180), breaks = c(-120, 0, 120)) +
    ggplot2::scale_y_continuous(limits = c(-90, 90), breaks = c(-60, 0, 60)) +
    ggplot2::labs(title = sprintf("Lineage-support mixture | %s Ma", key),
      subtitle = sprintf("Plant200 terrain-resistance setting | %s active sampled lineages", active),
      x = "Palaeolongitude", y = "Palaeolatitude",
      caption = "exp(Shannon) of support weights, not historical species richness. Ocean white; undefined land grey.") +
    ggplot2::theme_minimal(base_size = 12, base_family = "Arial") +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.position = "right", plot.title = ggplot2::element_text(face = "bold",
      colour = "#243B50"), plot.subtitle = ggplot2::element_text(colour = "#627681"),
      plot.caption = ggplot2::element_text(colour = "#53656F", size = 8))
  ggplot2::ggsave(png, p, width = 9, height = 5.4, dpi = 220,
                  device = ragg::agg_png, bg = "white")
  data.frame(time_ma = time_ma, active_lineages = active,
    land_cells = sum(land), defined_cells = sum(is.finite(effective)),
    max_effective_support = max(effective, na.rm = TRUE), png = png,
    geotiff = tif, stringsAsFactors = FALSE)
}))
utils::write.csv(div_index, file.path(out, "corrected_diversity_map_index.csv"),
                 row.names = FALSE)

file_url <- function(path) paste0("file:///", normalizePath(path, winslash = "/"))
common <- function(metric, scheme, focal, time_ma, status, png, geotiff) {
  data.frame(metric = metric, scheme = scheme, focal = focal,
    time_ma = as.numeric(time_ma), status = status,
    png = vapply(png, file_url, character(1)),
    geotiff = vapply(geotiff, file_url, character(1)),
    stringsAsFactors = FALSE)
}
loc <- shared[shared$metric == "02_log_density", , drop = FALSE]
move <- movement[movement$metric == "potential_land_corridors", , drop = FALSE]
ref <- refuges[refuges$metric == "candidate_refugia", , drop = FALSE]
records <- rbind(
  common("location", loc$scheme, "", loc$time_ma, "evaluated", loc$png, loc$geotiff),
  common("movement", "terrain_resistance", "", move$time_ma, "evaluated",
         move$png, move$geotiff),
  common("diversity", "terrain_resistance", "", div_index$time_ma, "evaluated",
         div_index$png, div_index$geotiff),
  common("refugia", "terrain_resistance", ref$focal_tip, ref$time_ma, ref$status,
         ref$png, ref$geotiff)
)
stopifnot(nrow(records) == 66L * (4L + 1L + 1L + 3L),
          all(file.exists(sub("^file:///", "", records$png))))
payload <- jsonlite::toJSON(records, dataframe = "rows", auto_unbox = TRUE,
                            na = "null", digits = NA)

html <- paste(c(
  '<!doctype html><html lang="en"><head><meta charset="utf-8">',
  '<meta name="viewport" content="width=device-width,initial-scale=1">',
  '<title>HmscEE Plant200 time-map atlas</title>',
  '<style>',
  ':root{font-family:Arial,sans-serif;color:#243746;background:#fff}',
  'body{margin:0}[hidden]{display:none!important}header{padding:20px max(20px,4vw);border-bottom:1px solid #d9e1e4}',
  'h1{font-size:24px;margin:0 0 6px}p{margin:4px 0;line-height:1.45;color:#53656f}',
  '.controls{display:flex;flex-wrap:wrap;gap:14px;align-items:end;padding:18px max(20px,4vw);background:#f4f8f7}',
  'label{display:grid;gap:6px;font-size:13px;font-weight:bold}select,input,button{font:inherit;padding:8px;border:1px solid #a9b9bd;border-radius:4px;background:#fff;color:#243746}',
  'button{cursor:pointer;min-width:38px}button:hover{background:#dbece6}',
  '.time{display:flex;gap:8px;align-items:center}.time input{width:min(40vw,280px)}',
  'main{padding:20px max(20px,4vw);max-width:1350px;margin:auto}',
  '#map{display:block;width:100%;max-height:74vh;object-fit:contain;background:#fff}',
  '.frame{border:1px solid #d9e1e4;border-radius:4px;padding:12px}',
  '.detail{display:flex;justify-content:space-between;gap:18px;flex-wrap:wrap;margin:12px 0}',
  '#status{font-weight:bold;color:#166f5a}a{color:#176a77}',
  '.note{border-left:4px solid #d38b27;padding:4px 14px;margin:16px 0}',
  '</style></head><body>',
  '<header><h1>Plant200 time-map atlas</h1>',
  '<p>All 66 saved palaeogeographic slices. The comparable modelled sequence ends at 5 Ma; 0 Ma is a separate record-conditioned endpoint diagnostic.</p></header>',
  '<div class="controls">',
  '<label>Result<select id="metric"><option value="location">Lineage location</option><option value="movement">Predictive movement</option><option value="diversity">Support-mixture diversity</option><option value="refugia">Candidate refugia</option></select></label>',
  '<label id="schemeBox">Movement setting<select id="scheme"><option value="terrain_resistance">Terrain resistance</option><option value="spherical_no_ldd">Spherical, no rare jumps</option><option value="spherical_land">Spherical land + rare jumps</option><option value="particle_topography">Finite propagule</option></select></label>',
  '<label id="focalBox" hidden>Ancestor path<select id="focal"><option value="Camellia_sinensis">Camellia sinensis</option><option value="Ilex_aquifolium">Ilex aquifolium</option><option value="Jasminum_officinale">Jasminum officinale</option></select></label>',
  '<label>Time<select id="time"></select></label>',
  '<div class="time"><button id="older" title="Older slice" aria-label="Older slice">←</button><input id="slider" type="range" min="0" max="65" step="1"><button id="younger" title="Younger slice" aria-label="Younger slice">→</button></div>',
  '</div><main><div class="detail"><span id="status"></span><span><a id="pngLink">Open PNG</a> · <a id="tifLink">Open GeoTIFF</a></span></div>',
  '<div class="frame"><img id="map" alt="Selected time-resolved Plant200 map"></div>',
  '<p class="note" id="interpretation"></p>',
  '<p>These maps are distinct evidence classes. White marine space is not a zero-diversity estimate. Candidate refugia are not verified population persistence.</p>',
  '</main><script>const records=__DATA__;',
  'const times=[...new Set(records.map(r=>r.time_ma))].sort((a,b)=>b-a);',
  'const $=id=>document.getElementById(id);',
  'const notes={location:"Tree-conditioned location density on land, not historical occupancy or richness. Log display is fixed across schemes and times.",movement:"Forward-predictive local-link throughput under a specified kernel, not posterior-used migration routes. One scale across all 66 times.",diversity:"Corrected effective number of lineage-support fields: exp(Shannon). Cells below numerical support precision are undefined, not absent. Not ancient species richness.",refugia:"Patch-history ecological candidates on the selected ancestral path. Orange additionally meets a modelled location-mass threshold. Neither class proves survival."};',
  'times.forEach(t=>{let o=document.createElement("option");o.value=t;o.textContent=t===0?"0 Ma | record endpoint":t+" Ma";$("time").append(o)});',
  '$("time").value="5";$("slider").value=times.indexOf(5);',
  'function update(){let metric=$("metric").value,time=Number($("time").value);',
  '$("schemeBox").hidden=!(metric==="location");$("focalBox").hidden=!(metric==="refugia");',
  'let scheme=metric==="location"?$("scheme").value:"terrain_resistance";',
  'let focal=metric==="refugia"?$("focal").value:"";',
  'let row=records.find(r=>r.metric===metric&&r.scheme===scheme&&r.focal===focal&&r.time_ma===time);',
  '$("slider").value=times.indexOf(time);$("interpretation").textContent=notes[metric]+(time===0?" At 0 Ma, exact conditioning on recorded terminal bins makes unrecorded land appear blank; this is not confirmed absence or a biological collapse.":"");',
  'if(!row||row.status!=="evaluated"){$("map").hidden=true;$("status").textContent="Not evaluable at this boundary slice";$("pngLink").hidden=true;$("tifLink").hidden=true;return}',
  '$("map").hidden=false;$("map").src=row.png;$("status").textContent=time+" Ma · "+metric.replaceAll("_"," ")+" · "+(focal||scheme);',
  '$("pngLink").hidden=false;$("tifLink").hidden=false;$("pngLink").href=row.png;$("tifLink").href=row.geotiff}',
  '$("metric").onchange=update;$("scheme").onchange=update;$("focal").onchange=update;$("time").onchange=update;',
  '$("slider").oninput=()=>{$("time").value=times[Number($("slider").value)];update()};',
  '$("older").onclick=()=>{$("time").value=times[Math.max(0,Number($("slider").value)-1)];update()};',
  '$("younger").onclick=()=>{$("time").value=times[Math.min(times.length-1,Number($("slider").value)+1)];update()};',
  'update();</script></body></html>'
), collapse = "\n")
html <- sub("__DATA__", payload, html, fixed = TRUE)
writeLines(html, file.path(out, "index.html"), useBytes = TRUE)
message("Saved corrected 66-time diversity atlas and browser: ", out)
