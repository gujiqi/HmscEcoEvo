suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}

ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(ofile) || length(ofile) == 0L || is.na(ofile) || !nzchar(ofile)) {
  ofile <- "scripts/rebuild_case01_complete_report.R"
}
script_path <- normalizePath(ofile, winslash = "/", mustWork = FALSE)
pkg_root <- normalizePath(file.path(dirname(script_path), ".."),
                          winslash = "/", mustWork = FALSE)
if (file.exists(file.path(pkg_root, "DESCRIPTION")) &&
    requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(pkg_root, quiet = TRUE)
} else {
  library(HmscEcoEvo)
}

output_root <- normalizePath(
  arg_value("output", file.path(pkg_root, "outputs/case01_simulated_full_540Ma")),
  winslash = "/", mustWork = TRUE
)
fig_dir <- file.path(output_root, "figures")
tab_dir <- file.path(output_root, "tables")
rds_dir <- file.path(output_root, "rds")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, name, width = 7.2, height = 4.8) {
  hee_save_plot(plot, file.path(fig_dir, tools::file_path_sans_ext(name)),
                width = width, height = height, dpi = 320)
}
write_table <- function(x, name) {
  path <- file.path(tab_dir, name)
  utils::write.csv(as.data.frame(x), path, row.names = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

legacy_times <- as.numeric(strsplit(arg_value("times", "420,240,180,120,60"),
                                    ",", fixed = TRUE)[[1]])
p_path <- file.path(rds_dir, "M5_env_phylo_bgb_dynamic_dispersal_probability_full.rds")
r_path <- file.path(rds_dir, "M5_env_phylo_bgb_dynamic_dispersal_richness_full.rds")
e_path <- file.path(rds_dir, "extrapolation_full.rds")
for (p in c(p_path, r_path, e_path)) {
  if (!file.exists(p)) stop("Missing required saved output: ", p, call. = FALSE)
}

message("Loading saved full projection RDS files")
P <- readRDS(p_path)
rich <- readRDS(r_path)
extrap <- readRDS(e_path)
species <- sort(unique(as.character(P$species)))
if (length(species) == 0L) stop("Saved projection has no species.", call. = FALSE)
focal_species <- species[1]
pal <- c("#163A5F", "#2F6DB5", "#2D7F7B", "#2D7F5E", "#D8A03D",
         "#C76D2A", "#B23A3A", "#6B4FA3", "#6E7781")

notes <- list()
for (tm in legacy_times) {
  if (!tm %in% unique(P$time_ma)) next
  message("Rebuilding legacy representative maps for ", tm, " Ma")
  one_sp <- P[P$time_ma == tm & P$species == focal_species, , drop = FALSE]
  one_rich <- rich[rich$time_ma == tm, , drop = FALSE]
  one_ex <- extrap[extrap$time_ma == tm, , drop = FALSE]
  acc_tm <- stats::aggregate(accessibility ~ region + time_ma,
                             P[P$time_ma == tm, , drop = FALSE],
                             mean, na.rm = TRUE)
  phy_tm <- stats::aggregate(phylo_existence ~ species,
                             P[P$time_ma == tm, , drop = FALSE],
                             max, na.rm = TRUE)
  names(phy_tm)[2] <- "valid"

  save_plot(ggplot(one_sp, aes(x = lon, y = lat, color = suitability)) +
              geom_point(size = 0.3) + coord_equal() +
              scale_color_gradient(low = "white", high = pal[2]) +
              labs(title = paste0("Suitability map, ", tm, " Ma"),
                   x = "Longitude", y = "Latitude", color = "Suitability"),
            sprintf("deep_suitability_map_%sMa.png", tm))
  save_plot(ggplot(acc_tm, aes(x = region, y = accessibility, fill = region)) +
              geom_col(show.legend = FALSE) +
              labs(title = paste0("Accessibility summary, ", tm, " Ma"),
                   x = "Region", y = "Accessibility"),
            sprintf("deep_accessibility_map_%sMa.png", tm))
  save_plot(ggplot(phy_tm, aes(x = species, y = valid)) +
              geom_col(fill = pal[4]) + coord_flip() +
              labs(title = paste0("Phylo mask summary, ", tm, " Ma"),
                   x = NULL, y = "E_phylo"),
            sprintf("deep_phylo_mask_summary_%sMa.png", tm))
  save_plot(ggplot(one_sp, aes(x = lon, y = lat, color = probability)) +
              geom_point(size = 0.3) + coord_equal() +
              scale_color_gradient(low = "white", high = pal[7]) +
              labs(title = paste0("Final probability example map, ", tm, " Ma"),
                   x = "Longitude", y = "Latitude", color = "Probability"),
            sprintf("deep_final_probability_example_map_%sMa.png", tm))
  save_plot(ggplot(one_rich, aes(x = lon, y = lat, color = expected_richness)) +
              geom_point(size = 0.3) + coord_equal() +
              scale_color_gradient(low = "white", high = pal[7]) +
              labs(title = paste0("Richness map, ", tm, " Ma"),
                   x = "Longitude", y = "Latitude", color = "Expected richness"),
            sprintf("deep_richness_map_%sMa.png", tm))

  unc <- stats::aggregate(probability ~ cell_id + time_ma + lon + lat,
                          P[P$time_ma == tm, , drop = FALSE],
                          function(p) sqrt(sum(pmax(p, 0) * pmax(1 - p, 0),
                                               na.rm = TRUE)))
  names(unc)[5] <- "uncertainty"
  save_plot(ggplot(unc, aes(x = lon, y = lat, color = uncertainty)) +
              geom_point(size = 0.3) + coord_equal() +
              scale_color_gradient(low = "white", high = pal[8]) +
              labs(title = paste0("Conditional richness uncertainty proxy, ",
                                  tm, " Ma"),
                   x = "Longitude", y = "Latitude", color = "Proxy"),
            sprintf("deep_uncertainty_map_%sMa.png", tm))
  save_plot(ggplot(one_ex, aes(x = lon, y = lat,
                               color = extrapolation_score)) +
              geom_point(size = 0.3) + coord_equal() +
              scale_color_gradient(low = "white", high = pal[6]) +
              labs(title = paste0("Extrapolation risk map, ", tm, " Ma"),
                   x = "Longitude", y = "Latitude",
                   color = "Extrapolation"),
            sprintf("deep_extrapolation_risk_map_%sMa.png", tm))
  notes[[length(notes) + 1L]] <- data.frame(
    time_ma = tm,
    focal_species = focal_species,
    restored_from = "saved full M5 projection/richness/extrapolation RDS",
    uncertainty_note = "deep_uncertainty_map uses conditional richness variance proxy from saved probabilities; future full script runs compute posterior-draw SD when beta_draws are available.",
    stringsAsFactors = FALSE
  )
}
if (length(notes) > 0L) {
  write_table(do.call(rbind, notes), "legacy_representative_map_restore_notes.csv")
}

fig_count <- length(list.files(fig_dir, pattern = "\\.png$"))
message("Rebuilding complete Word report from ", fig_count, " PNG figures")
report <- hee_write_full_workflow_report(
  output_root,
  output_file = file.path(output_root, "reports",
                          "hmscecoevo_full_540Ma_workflow.docx"),
  include_all_figures = TRUE,
  max_table_rows = 8
)
tab_count <- length(list.files(tab_dir, pattern = "\\.csv$"))
summary_path <- file.path(tab_dir, "case01_run_summary.csv")
if (file.exists(summary_path)) {
  s <- utils::read.csv(summary_path, check.names = FALSE)
  s$n_figures_png <- fig_count
  s$n_tables <- tab_count
  s$report_docx <- report
  utils::write.csv(s, summary_path, row.names = FALSE)
}
hee_write_ordered_outputs(output_root, include_rds = TRUE, overwrite = TRUE)
cat("Complete report rebuilt:\n", report, "\n", sep = "")
