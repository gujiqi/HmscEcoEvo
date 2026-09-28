#!/usr/bin/env Rscript

# Recompose the original Case07 diagnostic plots without their framed page layout.
# The saved metric objects, rather than the rendered PNGs, supply every panel.

suppressPackageStartupMessages({
  library(HmscEcoEvo)
  library(ggplot2)
  library(grid)
  library(ragg)
})

root <- "C:/Users/Google/Documents/HMSC-HIST"
case_dir <- file.path(root, "HmscEcoEvo", "inst", "extdata", "ecoevo_cases",
                      "outputs", "case07_validation_dashboard")
out_dir <- file.path(root, "outputs", "HmscEcoEvo",
                     "manuscript_results_six_process_20260927",
                     "hmscHist_supplement_clean_20260928")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

metrics <- readRDS(file.path(case_dir, "metrics.rds"))
spec <- list(
  list(number = 13L, key = "niche", method = "niche_summary",
       title = "Environmental response shape",
       keep = c(D = "Optimum and breadth", F = "Multivariate response direction")),
  list(number = 14L, key = "phylo", method = "phylo_signal",
       title = "Where phylogenetic signal occurs",
       keep = c(B = "Signal by phylogenetic distance",
                C = "Scale and clade conservatism",
                D = "Tip-level local signal", E = "Internal-node signal")),
  list(number = 15L, key = "evo", method = "evo_transition",
       title = "Evolutionary transitions in response history",
       keep = c(A = "Branch shift support", B = "Shift magnitude by response axis",
                C = "Shifts and rates on the tree",
                D = "Branch and clade rate diagnostics")),
  list(number = 16L, key = "trait", method = "trait_mediation",
       title = "Trait explanation and residual history",
       keep = c(B = "Phylogenetic signal decomposition",
                C = "Trait explanation and residual signal",
                D = "Trait phylogenetic redundancy",
                E = "Trait omission sensitivity")),
  list(number = 17L, key = "gamma", method = "gamma_evolution",
       title = "Variation in trait-response relationships",
       keep = c(A = "Clade-specific Gamma", B = "Regime-specific Gamma",
                C = "Gamma shifts on the tree",
                D = "Sign reversals and between-clade variance")),
  list(number = 18L, key = "population", method = "population_evolution",
       title = "Population-level response diagnostics",
       keep = c(B = "Within-species niche divergence",
                C = "Sources of response variation",
                D = "Local-adaptation to plasticity ratio",
                E = "Supplied genotype reaction norms")),
  list(number = 19L, key = "validation", method = "validation_dashboard",
       title = "Validation of derived hmscHist metrics",
       keep = c(A = "Recovery of derived metrics",
                C = "Tree uncertainty sensitivity",
                D = "Prior, environment and spatial sensitivity"))
)

panel_layout <- function(n) {
  if (n == 2L) {
    return(data.frame(x = c(0.045, 0.525), y = c(0.13, 0.13),
                      w = c(0.43, 0.43), h = c(0.70, 0.70)))
  }
  if (n == 3L) {
    return(data.frame(x = c(0.045, 0.045, 0.525),
                      y = c(0.47, 0.08, 0.08),
                      w = c(0.91, 0.43, 0.43),
                      h = c(0.35, 0.32, 0.32)))
  }
  stopifnot(n == 4L)
  data.frame(x = c(0.045, 0.525, 0.045, 0.525),
             y = c(0.47, 0.47, 0.08, 0.08),
             w = rep(0.43, 4), h = rep(0.33, 4))
}

draw_panel <- function(source_panel, label, title, placement) {
  pushViewport(viewport(x = placement$x, y = placement$y,
                        width = placement$w, height = placement$h,
                        just = c("left", "bottom")))
  grid.text(label, x = 0, y = 0.99, just = c("left", "top"),
            gp = gpar(fontfamily = "Arial", fontface = "bold",
                      fontsize = 11, col = "#19354A"))
  grid.text(title, x = 0.065, y = 0.99, just = c("left", "top"),
            gp = gpar(fontfamily = "Arial", fontface = "bold",
                      fontsize = 9.2, col = "#19354A"))
  pushViewport(viewport(x = 0, y = 0.005, width = 1, height = 0.87,
                        clip = "on",
                        just = c("left", "bottom")))
  plot <- source_panel$plot
  if (inherits(plot, "ggplot")) {
    print(plot, newpage = FALSE)
  } else if (inherits(plot, "grob")) {
    grid.draw(plot)
  } else {
    stop("Unsupported panel type: ", paste(class(plot), collapse = ", "))
  }
  popViewport(2)
}

draw_figure <- function(item, page) {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  grid.text(paste0("S", item$number, "  ", item$title),
            x = 0.045, y = 0.955, just = c("left", "top"),
            gp = gpar(fontfamily = "Arial", fontface = "bold",
                      fontsize = 13, col = "#19354A"))
  grid.text("Controlled software example; not a Plant200 empirical estimate",
            x = 0.045, y = 0.909, just = c("left", "top"),
            gp = gpar(fontfamily = "Arial", fontsize = 8.8,
                      col = "#586977"))
  positions <- panel_layout(length(item$keep))
  for (i in seq_along(item$keep)) {
    old_label <- names(item$keep)[i]
    source_panel <- page$panels[[match(old_label,
                                      vapply(page$panels, `[[`, character(1), "label"))]]
    if (is.null(source_panel)) stop("Missing source panel: ", old_label)
    draw_panel(source_panel, LETTERS[i], unname(item$keep[i]), positions[i, ])
  }
}

trace <- data.frame(supplement = character(), original_figure = character(),
                    original_panel = character(), revised_panel = character(),
                    decision = character(), reason = character(),
                    stringsAsFactors = FALSE)
removed_reasons <- c(
  "13A" = "Direct HMSC beta-by-species display",
  "13B" = "Direct HMSC beta posterior summary",
  "13C" = "Ordinary HMSC environmental response curves",
  "13E" = "Generic ordination of tip beta coefficients",
  "14A" = "Global HMSC rho and standard overall signal summaries",
  "14F" = "Direct residual-beta heatmap",
  "15E" = "Additional synthetic diversification composite omitted for focus",
  "15F" = "Crowded convergence composite omitted for readability",
  "16A" = "Ordinary HMSC trait-by-response Gamma matrix",
  "16F" = "Direct residual-beta heatmap",
  "17E" = "Redundant trait-function coefficient display",
  "17F" = "Additional turnover composite omitted for focus",
  "18A" = "Generic species-versus-population coefficient comparison",
  "18F" = "Supplied feedback score lacks independent causal identification",
  "19B" = "Uninformative collapsed uncertainty histogram in example fixture",
  "19E" = "Routine HMSC spatial block cross-validation",
  "19F" = "Routine HMSC posterior predictive and MCMC diagnostics"
)

for (item in spec) {
  plot_fn <- get(paste0("plot_", item$method), asNamespace("HmscEcoEvo"))
  page <- plot_fn(metrics[[item$key]], style = "nature")
  if (item$number == 16L) {
    decomposition <- metrics$trait$signal_decomposition
    decomposition$component <- factor(
      decomposition$component,
      levels = c("TMNS", "HMNS", "THMNS", "HPNS"),
      labels = c("Traits", "History", "Combined", "Residual signal"))
    page$panels[[match("B", vapply(page$panels, `[[`, character(1), "label"))]]$plot <-
      ggplot(decomposition, aes(axis, component, fill = value)) +
      geom_tile(colour = "white", linewidth = 0.7) +
      geom_text(aes(label = sprintf("%.2f", value)), size = 2.5) +
      scale_fill_gradient(low = "#F1F5F9", high = "#237F82", limits = c(0, 1),
                          name = "Index") +
      labs(x = "Environmental response axis", y = NULL) +
      theme_classic(base_size = 8) +
      theme(panel.grid = element_blank(), legend.position = "right")
  }
  if (item$number %in% c(14L, 15L, 17L)) {
    tree_label <- if (item$number == 14L) "E" else "C"
    tree_index <- match(tree_label,
                        vapply(page$panels, `[[`, character(1), "label"))
    page$panels[[tree_index]]$plot <-
      page$panels[[tree_index]]$plot + theme(legend.position = "none")
  }
  if (item$number == 15L) {
    index <- match("B", vapply(page$panels, `[[`, character(1), "label"))
    page$panels[[index]]$plot <-
      page$panels[[index]]$plot + theme(legend.position = "none")
  }
  if (item$number == 19L) {
    index <- match("A", vapply(page$panels, `[[`, character(1), "label"))
    recovery <- page$panels[[index]]$plot
    recovery$data <- recovery$data[recovery$data$metric != "gamma", , drop = FALSE]
    page$panels[[index]]$plot <- recovery
    index <- match("C", vapply(page$panels, `[[`, character(1), "label"))
    tree_sensitivity <- page$panels[[index]]$plot
    tree_sensitivity$data <- tree_sensitivity$data[
      tree_sensitivity$data$metric != "gamma", , drop = FALSE]
    page$panels[[index]]$plot <- tree_sensitivity
  }
  base <- file.path(fig_dir, paste0("Supplementary_Figure_S", item$number,
                                   "_hmscHist_clean"))
  ragg::agg_png(paste0(base, ".png"), width = 210, height = 145,
                units = "mm", res = 350, background = "white")
  draw_figure(item, page)
  dev.off()
  grDevices::cairo_pdf(paste0(base, ".pdf"), width = 210 / 25.4,
                       height = 145 / 25.4, family = "Arial")
  draw_figure(item, page)
  dev.off()
  labels <- vapply(page$panels, `[[`, character(1), "label")
  labels <- labels[nzchar(labels)]
  for (old in labels) {
    kept <- old %in% names(item$keep)
    trace <- rbind(trace, data.frame(
      supplement = paste0("S", item$number),
      original_figure = basename(file.path(case_dir,
        paste0(item$method, ".png"))),
      original_panel = old,
      revised_panel = if (kept) LETTERS[match(old, names(item$keep))] else "",
      decision = if (kept) "retained" else "removed",
      reason = if (kept) "hmscHist derived diagnostic"
               else unname(removed_reasons[paste0(item$number, old)]),
      stringsAsFactors = FALSE))
  }
  cat("Wrote S", item$number, ": ", paste(names(item$keep), collapse = ", "), "\n",
      sep = "")
}
write.csv(trace, file.path(out_dir, "S13-S19_original_panel_trace.csv"),
          row.names = FALSE)
