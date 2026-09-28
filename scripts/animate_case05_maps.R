#!/usr/bin/env Rscript

# Render time-ordered Case05 rasters as MP4s with a matched spatial summary.
# Run from the HmscEcoEvo package root. Set HMSCEE_ANIM_TEST=1 for a short check.

suppressPackageStartupMessages({
  library(terra)
  library(av)
})

default_root <- file.path(getwd(), "..", "outputs", "HmscEcoEvo")
root <- normalizePath(Sys.getenv("HMSCEE_OUTPUT_ROOT", default_root),
                      winslash = "/", mustWork = TRUE)
case <- file.path(root, "case05_v8_four_scheme_comparison_20260926")
main_index_path <- file.path(case, "07_cross_scheme_shared_scale_maps",
                             "cross_scheme_shared_scale_map_index.csv")
selected_index_path <- file.path(case, "10_patch_history_results_corrected",
                                 "case05_patch_map_index.csv")
difference_index_path <- file.path(case, "four_scheme_map_index.csv")
corrected_dir <- file.path(root, "manuscript_results_six_process_20260927",
                           "time_map_atlas", "diversity_corrected")
output <- file.path(root, "case05_plant200_time_animations_20260928")
test_mode <- identical(Sys.getenv("HMSCEE_ANIM_TEST"), "1")
if (test_mode) output <- paste0(output, "_preflight")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output, "videos"), showWarnings = FALSE)
dir.create(file.path(output, "posters"), showWarnings = FALSE)

read_index <- function(path) {
  if (!file.exists(path)) stop("Missing index: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

main <- read_index(main_index_path)
main$category <- "full_66"
main$group_id <- paste(main$category, main$scheme, main$metric, sep = "__")
main$path <- main$geotiff
main$kind <- "continuous"

selected <- read_index(selected_index_path)
selected$scheme <- ifelse(!is.na(selected$focal_tip) & nzchar(selected$focal_tip),
                          selected$focal_tip, "all_lineages")
selected$metric <- selected$metric
selected$category <- "sampled_3_or_6"
selected$group_id <- paste(selected$category, selected$scheme,
                           selected$metric, sep = "__")
selected$path <- selected$geotiff
selected$kind <- ifelse(selected$metric == "candidate_ecological_refugia_evidence_class",
                        "refugia", ifelse(selected$metric == "modern_training_noanalog_flag",
                                           "binary", "continuous"))
selected$scale_min <- NA_real_
selected$scale_max <- NA_real_

differences <- read_index(difference_index_path)
differences$scheme <- differences$comparison
differences$metric <- "location_density_difference"
differences$category <- "sampled_scheme_difference"
differences$group_id <- paste(differences$category, differences$scheme,
                              differences$metric, sep = "__")
differences$path <- differences$geotiff
differences$kind <- "difference"

corrected_files <- list.files(corrected_dir, "^effective_support_[0-9]+Ma[.]tif$",
                              full.names = TRUE)
corrected <- data.frame(
  category = "corrected_diversity_66",
  scheme = "manuscript_corrected",
  metric = "effective_support",
  time_ma = as.numeric(sub("^effective_support_([0-9]+)Ma[.]tif$",
                           "\\1", basename(corrected_files))),
  path = corrected_files,
  scale_min = NA_real_,
  scale_max = NA_real_,
  kind = "continuous",
  stringsAsFactors = FALSE
)
corrected$group_id <- paste(corrected$category, corrected$scheme,
                            corrected$metric, sep = "__")

cols <- c("category", "scheme", "metric", "time_ma", "path",
          "scale_min", "scale_max", "kind", "group_id")
rows <- rbind(main[, cols], selected[, cols], differences[, cols],
              corrected[, cols])
rows$time_ma <- as.numeric(rows$time_ma)
rows$scale_min <- as.numeric(rows$scale_min)
rows$scale_max <- as.numeric(rows$scale_max)
rows <- rows[order(rows$group_id, -rows$time_ma), ]
if (any(!file.exists(rows$path))) {
  stop("A listed GeoTIFF is missing: ", rows$path[which(!file.exists(rows$path))[1]])
}

metric_text <- list(
  "01_density" = c("Lineage location density", "Expected selected-lineage density",
                   "A lineage-location summary, not observed species richness."),
  "02_log_density" = c("Log lineage location density", "log10 density",
                       "The logarithmic view reveals low-density spatial support."),
  "03_ancestral_eta" = c("Ancestral environmental response", "HMSC linear predictor",
                        "Posterior-location-weighted response, not occupancy."),
  "04_land_fraction" = c("Available land fraction", "Land fraction",
                         "An external palaeogeographic condition, not dispersal."),
  "05_support_intensity" = c("Lineage-support intensity", "Support intensity",
                             "Terminal-calibrated model support, not a census."),
  "06_support_shannon" = c("Support-mixture Shannon index", "Shannon index",
                           "Diversity of modelled lineage support, not sampled community Shannon."),
  "07_support_effective" = c("Effective support lineages", "Effective lineages",
                             "exp(Shannon) of lineage-support mixture, not species richness."),
  "effective_support" = c("Corrected effective support diversity", "Effective lineages",
                          "Manuscript-corrected support-mixture diversity."),
  "active_lineage_environmental_support_count" =
    c("Environmentally supported active lineages", "Lineage count",
      "Count above the stated environmental-support threshold."),
  "active_lineage_environmental_support_fraction" =
    c("Fraction of active lineages supported", "Fraction",
      "Environmental opportunity, not realised occupancy."),
  "modern_training_noanalog_flag" =
    c("No-analog environments", "Flag (0-1)",
      "Outside the modern training envelope; not biological absence."),
  "ancestral_environmental_support" =
    c("Ancestral environmental support", "Support (0-1)",
      "Ancestral response matched to palaeoenvironment, not occupancy."),
  "candidate_ecological_refugia_evidence_class" =
    c("Candidate ecological refugia", "Highlighted land (%)",
      "Candidate patches, not confirmed population persistence."),
  "location_density_difference" =
    c("Dispersal-scheme location difference", "Density difference",
      "Signed comparison of modelled location density.")
)

scheme_text <- c(
  spherical_no_ldd = "Spherical movement, no rare long-distance jumps",
  spherical_land = "Spherical movement with rare long-distance jumps",
  terrain_resistance = "Terrain-resistance movement",
  particle_topography = "Finite-propagule terrain movement",
  manuscript_corrected = "Manuscript-corrected Plant200 result",
  all_lineages = "All active lineages"
)

make_palette <- function(metric) {
  if (metric %in% c("03_ancestral_eta", "location_density_difference")) {
    return(colorRampPalette(c("#255F8F", "#85B6D4", "#F7F7F3",
                              "#EAA368", "#AA3B3B"))(256))
  }
  if (metric == "04_land_fraction") {
    return(colorRampPalette(c("#DCEAD9", "#79B79B", "#1D755A"))(256))
  }
  if (metric %in% c("01_density", "02_log_density")) {
    return(colorRampPalette(c("#182C60", "#245C96", "#23A09A",
                              "#B9D97C", "#F5CF58"))(256))
  }
  if (metric %in% c("active_lineage_environmental_support_fraction",
                    "ancestral_environmental_support")) {
    return(colorRampPalette(c("#213F69", "#3C8FA7", "#86CBA3",
                              "#E9D87A"))(256))
  }
  colorRampPalette(c("#1D176B", "#5A1F94", "#B64E98",
                     "#F49757", "#F5D95C"))(256)
}

land_paths <- main[main$scheme == "spherical_land" &
                     main$metric == "04_land_fraction",
                   c("time_ma", "path")]
land_lookup <- setNames(land_paths$path, as.character(land_paths$time_ma))
land_cache <- new.env(parent = emptyenv())
get_land <- function(time_ma) {
  key <- as.character(time_ma)
  if (exists(key, envir = land_cache, inherits = FALSE)) {
    return(get(key, envir = land_cache, inherits = FALSE))
  }
  path <- land_lookup[[key]]
  if (is.null(path)) stop("No land mask for ", time_ma, " Ma")
  mask <- as.matrix(rast(path), wide = TRUE)
  mask <- ifelse(is.finite(mask) & mask > 0, 1, NA_real_)
  assign(key, mask, envir = land_cache)
  mask
}

read_grid <- function(path, time_ma) {
  grid <- as.matrix(rast(path), wide = TRUE)
  land <- get_land(time_ma)
  if (!identical(dim(grid), dim(land))) stop("Grid mismatch: ", path)
  grid[is.na(land)] <- NA_real_
  grid
}

latitudes <- seq(88, -88, by = -4)
cell_weight <- matrix(rep(pmax(cos(latitudes * pi / 180), 0), 90),
                      nrow = 45, ncol = 90)

weighted_quantile <- function(x, w, probs) {
  order_x <- order(x)
  x <- x[order_x]
  w <- w[order_x]
  cum <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cum >= p)[1]], numeric(1))
}

summarise_grid <- function(grid, kind) {
  keep <- is.finite(grid) & is.finite(cell_weight) & cell_weight > 0
  if (!any(keep)) return(c(mean = NA_real_, q25 = NA_real_,
                           q75 = NA_real_, secondary = NA_real_))
  x <- grid[keep]
  w <- cell_weight[keep]
  if (kind == "refugia") {
    return(c(mean = 100 * sum(w[x == 2]) / sum(w), q25 = NA_real_,
             q75 = NA_real_, secondary = 100 * sum(w[x == 1]) / sum(w)))
  }
  if (kind == "binary") x <- 100 * as.numeric(x > 0)
  if (kind == "difference") x <- abs(x)
  q <- weighted_quantile(x, w, c(0.25, 0.75))
  c(mean = sum(x * w) / sum(w), q25 = q[1], q75 = q[2],
    secondary = NA_real_)
}

summary_path <- file.path(output, "time_summary.csv")
if (file.exists(summary_path) && !test_mode) {
  summaries <- read.csv(summary_path, stringsAsFactors = FALSE)
} else {
  message("Summarising ", nrow(rows), " GeoTIFFs...")
  summary_rows <- vector("list", nrow(rows))
  for (i in seq_len(nrow(rows))) {
    stat <- summarise_grid(read_grid(rows$path[i], rows$time_ma[i]),
                           rows$kind[i])
    summary_rows[[i]] <- cbind(rows[i, c("category", "scheme", "metric",
                                        "time_ma", "group_id")],
                                as.data.frame(as.list(stat)))
  }
  summaries <- do.call(rbind, summary_rows)
  write.csv(summaries, summary_path, row.names = FALSE)
}

group_ids <- unique(rows$group_id)
if (test_mode) {
  group_ids <- c("full_66__spherical_land__07_support_effective",
                 "sampled_3_or_6__Camellia_sinensis__candidate_ecological_refugia_evidence_class",
                 "sampled_scheme_difference__terrain_resistance_minus_spherical_land__location_density_difference")
}
force_refugia <- identical(Sys.getenv("HMSCEE_ANIM_FORCE_REFUGIA"), "1")
force_all <- identical(Sys.getenv("HMSCEE_ANIM_FORCE_ALL"), "1")
force_pattern <- Sys.getenv("HMSCEE_ANIM_FORCE_PATTERN")
if (force_refugia && !test_mode) {
  group_ids <- group_ids[grepl("candidate_ecological_refugia_evidence_class",
                                group_ids, fixed = TRUE)]
}
if (nzchar(force_pattern) && !test_mode) {
  group_ids <- group_ids[grepl(force_pattern, group_ids, fixed = TRUE)]
  if (!length(group_ids)) stop("No animation group matches force pattern")
}

scale_for <- function(group) {
  low <- group$scale_min[is.finite(group$scale_min)]
  high <- group$scale_max[is.finite(group$scale_max)]
  if (length(low) && length(high)) return(c(min(low), max(high)))
  if (group$kind[1] == "refugia") return(c(0, 2))
  if (group$kind[1] == "binary") return(c(0, 1))
  if (group$metric[1] %in% c("ancestral_environmental_support",
                             "active_lineage_environmental_support_fraction")) {
    return(c(0, 1))
  }
  paths <- rows$path[rows$metric == group$metric[1] &
                       rows$category == group$category[1]]
  ranges <- vapply(paths, function(path) {
    value <- terra::global(rast(path), c("min", "max"), na.rm = TRUE)
    c(value[1, 1], value[1, 2])
  }, numeric(2))
  range(ranges, finite = TRUE)
}

draw_map <- function(grid, land, limits, palette, kind, unit) {
  par(mar = c(2.2, 3.2, 0.2, 1.2), xaxs = "i", yaxs = "i")
  x <- seq(-178, 178, by = 4)
  y <- seq(-88, 88, by = 4)
  oriented <- function(m) t(m[nrow(m):1, , drop = FALSE])
  image(x, y, oriented(land), col = "#E8EDEB", zlim = c(0, 1),
        xlim = c(-180, 244), ylim = c(-90, 90), asp = 1,
        axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
  if (kind == "refugia") {
    image(x, y, oriented(grid), col = c("#E8EDEB", "#198B75", "#DA8B29"),
          zlim = c(0, 2), add = TRUE, useRaster = TRUE)
  } else {
    image(x, y, oriented(grid), col = palette, zlim = limits,
          add = TRUE, useRaster = TRUE)
  }
  axis(1, at = c(-120, 0, 120), labels = c("120 W", "0", "120 E"),
       col = NA, col.axis = "#607587", cex.axis = 0.9)
  axis(2, at = c(-60, 0, 60), labels = c("60 S", "0", "60 N"),
       las = 1, col = NA, col.axis = "#607587", cex.axis = 0.9)
  if (kind == "refugia") {
    rect(191, 47, 203, 55, col = "#198B75", border = NA)
    text(210, 51, "Ecological candidate", pos = 4, cex = 0.9,
         col = "#22384C")
    rect(191, 25, 203, 33, col = "#DA8B29", border = NA)
    text(210, 29, "Candidate + location mass", pos = 4,
         cex = 0.9, col = "#22384C")
  } else {
    ymin <- -43
    ymax <- 55
    n <- length(palette)
    for (i in seq_len(n)) {
      a <- ymin + (i - 1) * (ymax - ymin) / n
      b <- ymin + i * (ymax - ymin) / n
      rect(193, a, 206, b, col = palette[i], border = NA)
    }
    for (j in 0:4) {
      v <- limits[1] + j * diff(limits) / 4
      yy <- ymin + j * (ymax - ymin) / 4
      segments(206, yy, 210, yy, col = "#607587")
      text(212, yy, formatC(v, digits = 2, format = "fg"),
           pos = 4, cex = 0.85, col = "#22384C")
    }
    text(194, 66, unit, adj = c(0, 0.5), cex = 0.85, col = "#22384C")
  }
  text(193, -67, "White: ocean", adj = c(0, 0.5),
       cex = 0.8, col = "#5E7180")
  grey_label <- if (kind == "refugia") {
    "Grey: not flagged as candidate"
  } else "Grey: land without value"
  text(193, -77, grey_label, adj = c(0, 0.5),
       cex = 0.8, col = "#5E7180")
}

draw_chart <- function(group, summaries, current_time, current_scheme) {
  par(mar = c(3.9, 4.8, 1.1, 1.7), xaxs = "i")
  own <- summaries[summaries$group_id == group$group_id[1], ]
  own <- own[order(-own$time_ma), ]
  is_full <- group$category[1] == "full_66"
  compared <- if (is_full) {
    summaries[summaries$category == "full_66" &
                summaries$metric == group$metric[1], ]
  } else own
  vals <- c(compared$mean, own$q25, own$q75, own$secondary)
  ylim <- range(vals, finite = TRUE)
  if (group$kind[1] %in% c("refugia", "binary")) ylim <- c(0, max(10, ylim[2]))
  if (!all(is.finite(ylim))) ylim <- c(0, 1)
  if (diff(ylim) < 1e-9) ylim <- ylim + c(-0.5, 0.5)
  ylim <- ylim + c(-1, 1) * diff(ylim) * 0.12
  plot(NA, xlim = c(325, 0), ylim = ylim, axes = FALSE,
       xlab = "", ylab = "", bty = "n")
  abline(h = pretty(ylim, 4), col = "#E7ECEF", lwd = 1)
  abline(v = c(300, 250, 200, 150, 100, 50, 0),
         col = "#EDF1F3", lwd = 1)
  axis(1, at = c(325, 300, 250, 200, 150, 100, 50, 0),
       col = NA, col.axis = "#51677B", cex.axis = 0.86)
  axis(2, las = 1, col = NA, col.axis = "#51677B", cex.axis = 0.86)
  mtext("Time before present (Ma)", side = 1, line = 2.6,
        cex = 0.86, col = "#344E62")
  ylabel <- if (group$kind[1] == "refugia") {
    "Land share in highlighted class (%)"
  } else if (group$kind[1] == "binary") {
    "No-analog land (%)"
  } else if (group$kind[1] == "difference") {
    "Land-area-weighted mean absolute difference"
  } else "Land-area-weighted spatial mean"
  mtext(ylabel, side = 2, line = 3.4, cex = 0.82, col = "#344E62")
  if (is_full) {
    scheme_colors <- c(spherical_no_ldd = "#3976AA",
                       spherical_land = "#D1902D",
                       terrain_resistance = "#219178",
                       particle_topography = "#A05798")
    for (scheme in names(scheme_colors)) {
      series <- compared[compared$scheme == scheme, ]
      series <- series[order(-series$time_ma), ]
      lines(series$time_ma, series$mean, col = adjustcolor(
        scheme_colors[[scheme]], alpha.f = if (scheme == current_scheme) 1 else 0.35),
        lwd = if (scheme == current_scheme) 3.5 else 1.8)
    }
  } else {
    main_color <- if (group$kind[1] == "refugia") "#DA8B29" else "#128F88"
    lines(own$time_ma, own$mean, col = main_color, lwd = 3,
          lty = if (nrow(own) < 20) 2 else 1)
    points(own$time_ma, own$mean, pch = 16, cex = 0.8,
           col = main_color)
    if (group$kind[1] == "refugia" && any(is.finite(own$secondary))) {
      lines(own$time_ma, own$secondary, col = "#198B75", lwd = 2.6, lty = 2)
      points(own$time_ma, own$secondary, col = "#198B75", pch = 16)
    }
  }
  abline(v = current_time, col = "#E2A531", lwd = 2)
  value <- own$mean[match(current_time, own$time_ma)]
  points(current_time, value, pch = 21, bg = "#E2A531",
         col = "white", lwd = 1.5, cex = 1.6)
  if (is_full) {
    legend("topright", legend = c("No LDD", "Rare LDD", "Terrain", "Particles"),
           col = c("#3976AA", "#D1902D", "#219178", "#A05798"),
           lwd = 2.5, horiz = TRUE, bty = "n", cex = 0.75)
  } else if (group$kind[1] == "refugia") {
    legend("topright", legend = c("Candidate + location mass",
                                  "Ecological candidate"),
           col = c("#DA8B29", "#198B75"), lwd = 2.5,
           horiz = TRUE, bty = "n", cex = 0.76)
  }
}

draw_frame <- function(group, grid, land, limits, palette, summaries,
                       current_time) {
  layout(matrix(1:3, ncol = 1), heights = c(1.05, 6.3, 2.65))
  par(mar = c(0, 1, 0, 1))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  label <- metric_text[[group$metric[1]]]
  title <- if (is.null(label)) gsub("_", " ", group$metric[1]) else label[1]
  subtitle <- unname(scheme_text[group$scheme[1]])
  if (is.na(subtitle)) subtitle <- gsub("_", " ", group$scheme[1])
  subtitle <- paste(subtitle, "|", nrow(group), "stored time slices")
  if (current_time == 0) {
    subtitle <- paste(subtitle, "| 0 Ma uses modern endpoint conditioning")
  }
  text(0.02, 0.76, title, adj = c(0, 0.5), cex = 1.8,
       font = 2, col = "#17344B")
  text(0.02, 0.29, subtitle, adj = c(0, 0.5), cex = 1.08,
       col = "#506779")
  rect(0.88, 0.18, 0.985, 0.84, col = "#E6F2EE", border = NA)
  text(0.932, 0.51, paste0(current_time, " Ma"), cex = 1.25,
       font = 2, col = "#187C6C")
  palette <- make_palette(group$metric[1])
  unit <- if (is.null(label)) "Value" else label[2]
  draw_map(grid, land, limits, palette, group$kind[1], unit)
  draw_chart(group, summaries, current_time, group$scheme[1])
  layout(1)
}

render_group <- function(group_id) {
  group <- rows[rows$group_id == group_id, ]
  group <- group[order(-group$time_ma), ]
  if (test_mode) group <- head(group, 3)
  if (anyDuplicated(group$time_ma)) stop("Duplicate time in ", group_id)
  limits <- scale_for(group)
  if (!all(is.finite(limits)) || diff(limits) <= 0) {
    stop("Invalid colour scale in ", group_id)
  }
  frame_dir <- file.path(tempdir(), gsub("[^A-Za-z0-9_-]", "_", group_id))
  dir.create(frame_dir, recursive = TRUE, showWarnings = FALSE)
  frame_paths <- file.path(frame_dir, sprintf("frame_%03d.png", seq_len(nrow(group))))
  for (i in seq_len(nrow(group))) {
    time_ma <- group$time_ma[i]
    ragg::agg_png(frame_paths[i], width = 1600, height = 900,
                  res = 125, background = "white")
    draw_frame(group, read_grid(group$path[i], time_ma),
               get_land(time_ma), limits, make_palette(group$metric[1]),
               summaries, time_ma)
    dev.off()
  }
  poster_time <- if (nrow(group) > 20) 100 else 65
  poster_i <- which.min(abs(group$time_ma - poster_time))
  poster_path <- file.path(output, "posters", paste0(group_id, ".png"))
  file.copy(frame_paths[poster_i], poster_path, overwrite = TRUE)
  video_path <- file.path(output, "videos", paste0(group_id, ".mp4"))
  frame_repeat <- if (nrow(group) > 20) 2L else if (nrow(group) > 3) 12L else 24L
  av::av_encode_video(rep(frame_paths, each = frame_repeat),
                      output = video_path, framerate = 10,
                      codec = "libx264", vfilter = "format=yuv420p",
                      verbose = FALSE)
  if (!file.exists(video_path) || file.info(video_path)$size < 10000) {
    stop("Video missing or too small: ", video_path)
  }
  safe_root <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
  safe_frames <- normalizePath(frame_dir, winslash = "/", mustWork = TRUE)
  if (!startsWith(safe_frames, paste0(safe_root, "/"))) {
    stop("Unsafe temporary frame directory")
  }
  unlink(frame_dir, recursive = TRUE)
  data.frame(group_id = group_id, category = group$category[1],
             scheme = group$scheme[1], metric = group$metric[1],
             n_time_slices = nrow(group),
             oldest_ma = max(group$time_ma),
             youngest_ma = min(group$time_ma),
             scale_min = limits[1], scale_max = limits[2],
             video = video_path, poster = poster_path,
             stringsAsFactors = FALSE)
}

index_path <- file.path(output, "animation_index.csv")
existing <- if (file.exists(index_path)) read.csv(index_path, stringsAsFactors = FALSE) else NULL
finished <- if (is.null(existing)) character() else
  existing$group_id[file.exists(existing$video)]
if (force_refugia) finished <- setdiff(finished, group_ids)
if (force_all) finished <- character()
if (nzchar(force_pattern)) finished <- setdiff(finished, group_ids)
for (id in group_ids) {
  if (id %in% finished) next
  message("Rendering ", id, " (", match(id, group_ids), "/", length(group_ids), ")")
  result <- render_group(id)
  if (!is.null(existing)) existing <- existing[existing$group_id != id, ]
  existing <- if (is.null(existing)) result else rbind(existing, result)
  write.csv(existing, index_path, row.names = FALSE)
}
message("Done: ", nrow(existing), " videos in ", output)
