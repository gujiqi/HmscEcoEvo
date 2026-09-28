# Nature-style figure utilities ---------------------------------------------

.nature_palette <- function() {
  list(
    navy = "#0B2A55",
    ink = "#111827",
    muted = "#4B5563",
    line = "#CBD5E1",
    fill = "#F8FAFC",
    blue = "#1F5FBF",
    light_blue = "#DDEBFF",
    red = "#D9272E",
    orange = "#F26B21",
    green = "#1B8A3A",
    purple = "#7A3DB8",
    brown = "#8B4A20",
    grey = "#6B7280",
    clade = c(A = "#1F5FBF", B = "#1B8A3A", C = "#F26B21",
              D = "#7A3DB8", E = "#D9272E"),
    axis = c(temp = "#D9272E", precip = "#1F5FBF", forest = "#1B8A3A",
             soil = "#8B4A20", pH = "#7A3DB8", pH_H2O = "#7A3DB8",
             temp_sq = "#D9272E", precip_sq = "#1F5FBF",
             forest_sq = "#1B8A3A", soil_sq = "#8B4A20",
             pH_sq = "#7A3DB8")
  )
}

.nature_axis_labels <- function() {
  c(temp = "Temp.", precip = "Precip.", forest = "Forest",
    soil = "Soil", pH = "pH", pH_H2O = "pH")
}

.nature_theme <- function(base_size = 7) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = .nature_palette()$ink),
      axis.line = ggplot2::element_line(linewidth = 0.28, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.28, colour = "black"),
      axis.text = ggplot2::element_text(size = base_size - 1),
      axis.title = ggplot2::element_text(size = base_size),
      plot.title = ggplot2::element_text(size = base_size + 1, face = "bold",
                                         hjust = 0),
      legend.title = ggplot2::element_text(size = base_size - 1),
      legend.text = ggplot2::element_text(size = base_size - 2),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = base_size - 1),
      panel.grid.major = ggplot2::element_line(linewidth = 0.15,
                                               colour = "#E5E7EB"),
      panel.grid.minor = ggplot2::element_blank()
    )
}

.nature_label_axes <- function(x) {
  labs <- .nature_axis_labels()
  unname(ifelse(x %in% names(labs), labs[x], x))
}

.nature_analysis_axes <- function(evo, include_quadratic = FALSE) {
  md_axes <- evo$metadata$analysis_axes
  if (!is.null(md_axes)) {
    axes <- md_axes
  } else {
    axes <- .non_intercept_axes(colnames(evo$beta))
  }
  if (!include_quadratic) axes <- axes[!grepl("_sq$|_squared$|\\^2$", axes)]
  axes[axes %in% colnames(evo$beta)]
}

.nature_clade <- function(evo) {
  cl <- evo$clade
  sp <- evo$species
  if (is.null(cl)) {
    z <- rep(LETTERS[1:min(5, length(sp))], length.out = length(sp))
    names(z) <- sp
    return(z)
  }
  cl <- as.character(cl)
  names(cl) <- names(evo$clade)
  keys <- unique(cl)
  map <- setNames(LETTERS[seq_along(keys)], keys)
  out <- unname(map[cl[sp]])
  names(out) <- sp
  out
}

.nature_plot_scales <- function() {
  pal <- .nature_palette()
  list(
    ggplot2::scale_colour_manual(values = c(pal$clade, pal$axis), na.value = pal$grey),
    ggplot2::scale_fill_manual(values = c(pal$clade, pal$axis), na.value = pal$grey)
  )
}

.tree_segments <- function(phylo, species_order = NULL, clade = NULL) {
  if (is.null(phylo)) return(list(edges = data.frame(), tips = data.frame()))
  .require_pkg("ape", "Nature-style phylogeny panels")
  phylo <- ape::reorder.phylo(phylo, "cladewise")
  ntip <- length(phylo$tip.label)
  nn <- phylo$Nnode
  depth <- ape::node.depth.edgelength(phylo)
  if (is.null(depth) || any(!is.finite(depth))) depth <- ape::node.depth(phylo)
  if (is.null(species_order)) species_order <- rev(phylo$tip.label)
  y <- rep(NA_real_, ntip + nn)
  names(y) <- as.character(seq_len(ntip + nn))
  tip_y <- setNames(seq_along(species_order), species_order)
  y[as.character(seq_len(ntip))] <- tip_y[phylo$tip.label]
  for (node in rev(unique(phylo$edge[, 1]))) {
    ch <- phylo$edge[phylo$edge[, 1] == node, 2]
    y[as.character(node)] <- mean(y[as.character(ch)], na.rm = TRUE)
  }
  node_clade <- rep(NA_character_, ntip + nn)
  node_clade[seq_len(ntip)] <- if (is.null(clade)) NA_character_ else clade[phylo$tip.label]
  for (node in rev(unique(phylo$edge[, 1]))) {
    ch <- phylo$edge[phylo$edge[, 1] == node, 2]
    vals <- unique(stats::na.omit(node_clade[ch]))
    node_clade[node] <- if (length(vals) == 1) vals else NA_character_
  }
  edges <- data.frame(
    branch = seq_len(nrow(phylo$edge)),
    parent = phylo$edge[, 1],
    child = phylo$edge[, 2],
    x = depth[phylo$edge[, 1]],
    xend = depth[phylo$edge[, 2]],
    y = y[as.character(phylo$edge[, 2])],
    yend = y[as.character(phylo$edge[, 2])],
    ym = y[as.character(phylo$edge[, 1])],
    clade = node_clade[phylo$edge[, 2]],
    stringsAsFactors = FALSE
  )
  node_ids <- unique(phylo$edge[, 1])
  verts <- data.frame(
    node = node_ids,
    x = depth[node_ids],
    y = vapply(node_ids, function(n) {
      ch <- phylo$edge[phylo$edge[, 1] == n, 2]
      min(y[as.character(ch)], na.rm = TRUE)
    }, numeric(1)),
    yend = vapply(node_ids, function(n) {
      ch <- phylo$edge[phylo$edge[, 1] == n, 2]
      max(y[as.character(ch)], na.rm = TRUE)
    }, numeric(1)),
    clade = node_clade[node_ids],
    stringsAsFactors = FALSE
  )
  nodes <- data.frame(
    node = node_ids,
    x = depth[node_ids],
    y = y[as.character(node_ids)],
    clade = node_clade[node_ids],
    stringsAsFactors = FALSE
  )
  tips <- data.frame(species = phylo$tip.label,
                     x = depth[seq_len(ntip)],
                     y = y[as.character(seq_len(ntip))],
                     clade = if (is.null(clade)) NA_character_ else clade[phylo$tip.label],
                     stringsAsFactors = FALSE)
  list(edges = edges, verts = verts, nodes = nodes, tips = tips,
       max_x = max(depth, na.rm = TRUE))
}

.plot_phylo_tree <- function(evo, tip_metric = NULL, node_metric = NULL,
                             time_axis = TRUE, title = NULL,
                             show_legend = FALSE,
                             fill_name = "value",
                             size_name = "magnitude") {
  pal <- .nature_palette()
  cl <- .nature_clade(evo)
  tree <- .tree_segments(evo$phylo, species_order = rev(evo$species), clade = cl)
  if (nrow(tree$edges) == 0) return(ggplot2::ggplot() + .nature_theme())
  tips <- tree$tips
  if (!is.null(tip_metric)) {
    tips <- merge(tips, tip_metric, by = "species", all.x = TRUE)
  }
  nodes <- tree$nodes
  if (!is.null(node_metric)) {
    nodes <- merge(nodes, node_metric, by = "node", all.x = TRUE)
  }
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = tree$edges,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend, colour = clade),
      linewidth = 0.45, na.rm = TRUE) +
    ggplot2::geom_segment(data = tree$verts,
      ggplot2::aes(x = x, xend = x, y = y, yend = yend, colour = clade),
      linewidth = 0.45, na.rm = TRUE) +
    ggplot2::geom_point(data = tips,
      ggplot2::aes(x = x, y = y, colour = clade),
      shape = 21, size = 1.7, fill = "white", stroke = 0.35) +
    ggplot2::geom_text(data = tips,
      ggplot2::aes(x = x + 0.03 * tree$max_x, y = y, label = species,
                   colour = clade),
      hjust = 0, size = 2.0) +
    ggplot2::scale_colour_manual(values = pal$clade, na.value = pal$grey) +
    ggplot2::guides(colour = "none") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(title = title, x = if (time_axis) "Time from root (Ma)" else NULL,
                  y = NULL) +
    .nature_theme(7) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   legend.position = if (show_legend) "right" else "none",
                   plot.margin = ggplot2::margin(2, 12, 2, 2))
  if (!is.null(tip_metric) && "value" %in% names(tips)) {
    p <- p + ggplot2::geom_point(data = tips,
      ggplot2::aes(x = x, y = y, size = abs(value), fill = value),
      shape = 21, colour = "black", stroke = 0.2, inherit.aes = FALSE) +
      ggplot2::scale_size_continuous(range = c(1, 4), guide = "none") +
      ggplot2::scale_fill_gradient2(low = "#1F5FBF", mid = "white",
                                    high = "#D9272E", midpoint = 0,
                                    na.value = "grey85", name = fill_name)
  }
  if (!is.null(node_metric) && "value" %in% names(nodes)) {
    if (!"n_species" %in% names(nodes)) nodes$n_species <- abs(nodes$value)
    p <- p + ggplot2::geom_point(data = nodes,
      ggplot2::aes(x = x, y = y, size = n_species, fill = value),
      shape = 21, colour = "black", stroke = 0.25, inherit.aes = FALSE,
      na.rm = TRUE) +
      ggplot2::scale_size_continuous(range = c(1.2, 4.5), name = size_name) +
      ggplot2::scale_fill_gradient2(low = "#1F5FBF", mid = "white",
                                    high = "#D9272E", midpoint = 0,
                                    na.value = "grey85", name = fill_name) +
      ggplot2::theme(legend.position = if (show_legend) "right" else "none")
  }
  p
}

.ggplot_grid_grob <- function(plots, positions) {
  grobs <- Map(function(p, i) {
    g <- ggplot2::ggplotGrob(p)
    grid::grobTree(g, vp = grid::viewport(x = positions$x[i], y = positions$y[i],
                                          width = positions$w[i], height = positions$h[i],
                                          just = c("left", "bottom"),
                                          clip = "on"))
  }, plots, seq_along(plots))
  do.call(grid::grobTree, grobs)
}

.p_signif_label <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.001, "***",
    ifelse(p < 0.01, "**",
    ifelse(p < 0.05, "*", ""))))
}

.tree_heatmap_plot <- function(evo, matrix, title = NULL, fill_title = "Beta",
                               limits = NULL) {
  pal <- .nature_palette()
  axes <- colnames(matrix)
  cl <- .nature_clade(evo)
  tree <- .tree_segments(evo$phylo, species_order = rev(rownames(matrix)), clade = cl)
  tree_scale <- 0.55
  tree$edges$x <- tree$edges$x * tree_scale
  tree$edges$xend <- tree$edges$xend * tree_scale
  tree$verts$x <- tree$verts$x * tree_scale
  tree$tips$x <- tree$tips$x * tree_scale
  tree$max_x <- tree$max_x * tree_scale
  beta_long <- .matrix_to_long(matrix, "value", "species", "axis")
  beta_long$y <- tree$tips$y[match(beta_long$species, tree$tips$species)]
  heat_gap <- max(4.0, tree$max_x * 0.120)
  heat_step <- max(0.55, tree$max_x * 0.028)
  beta_long$x <- tree$max_x + heat_gap + match(beta_long$axis, axes) * heat_step
  axis_df <- data.frame(axis = axes,
                        x = tree$max_x + heat_gap + seq_along(axes) * heat_step,
                        y = max(tree$tips$y, na.rm = TRUE) + 0.8)
  lim <- limits %||% max(abs(beta_long$value), na.rm = TRUE)
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = tree$edges,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend, colour = clade),
      linewidth = 0.42, na.rm = TRUE) +
    ggplot2::geom_segment(data = tree$verts,
      ggplot2::aes(x = x, xend = x, y = y, yend = yend, colour = clade),
      linewidth = 0.42, na.rm = TRUE) +
    ggplot2::geom_text(data = tree$tips,
      ggplot2::aes(x = tree$max_x + heat_gap * 0.28, y = y, label = species),
      hjust = 0, size = 1.9) +
    ggplot2::geom_tile(data = beta_long,
      ggplot2::aes(x = x, y = y, fill = value),
      width = heat_step * 0.86, height = 0.80, colour = "white", linewidth = 0.15) +
    ggplot2::geom_text(data = axis_df,
      ggplot2::aes(x = x, y = y, label = .nature_label_axes(axis)),
      angle = 90, size = 1.9, fontface = "bold", vjust = 0.5) +
    ggplot2::scale_colour_manual(values = pal$clade, na.value = pal$grey) +
    ggplot2::scale_fill_gradient2(low = "#1F5FBF", mid = "white",
                                  high = "#D9272E", midpoint = 0,
                                  limits = c(-lim, lim),
                                  name = fill_title, na.value = "grey90") +
    ggplot2::coord_cartesian(xlim = c(0, tree$max_x + heat_gap +
      (length(axes) + 1) * heat_step), clip = "off") +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    .nature_theme(7) +
    ggplot2::guides(colour = "none") +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank(),
                   legend.position = "right",
                   plot.margin = ggplot2::margin(2, 8, 2, 2))
}

.panel_title_grob <- function(label, title, navy) {
  grid::grobTree(
    grid::roundrectGrob(x = grid::unit(0.02, "npc"), y = grid::unit(0.965, "npc"),
                        width = grid::unit(0.045, "npc"),
                        height = grid::unit(0.055, "npc"),
                        r = grid::unit(0.006, "npc"),
                        just = c("left", "top"),
                        gp = grid::gpar(fill = navy, col = navy, lwd = 0.3)),
    grid::textGrob(label, x = grid::unit(0.043, "npc"), y = grid::unit(0.938, "npc"),
                   gp = grid::gpar(col = "white", fontsize = 9, fontface = "bold")),
    grid::textGrob(title, x = grid::unit(0.085, "npc"), y = grid::unit(0.945, "npc"),
                   just = c("left", "center"),
                   gp = grid::gpar(col = navy, fontsize = 8.5, fontface = "bold"))
  )
}

.draw_nature_panel <- function(panel, navy) {
  grid::pushViewport(grid::viewport(x = panel$x, y = panel$y,
                                    width = panel$w, height = panel$h,
                                    just = c("left", "bottom")))
  grid::grid.roundrect(r = grid::unit(0.012, "npc"),
                       gp = grid::gpar(fill = "white", col = navy, lwd = 1.0))
  grid::grid.draw(.panel_title_grob(panel$label, panel$title, navy))
  note_h <- if (!is.null(panel$note)) 0.12 else 0.02
  grid::pushViewport(grid::viewport(x = 0.03, y = 0.05 + note_h,
                                    width = 0.94, height = 0.83 - note_h,
                                    just = c("left", "bottom")))
  if (inherits(panel$plot, "ggplot")) {
    print(panel$plot, newpage = FALSE)
  } else if (inherits(panel$plot, "grob")) {
    grid::grid.draw(panel$plot)
  }
  grid::popViewport()
  if (!is.null(panel$note)) {
    grid::grid.roundrect(x = 0.5, y = 0.065, width = 0.90, height = 0.10,
                         r = grid::unit(0.006, "npc"),
                         gp = grid::gpar(fill = "#F8FAFC", col = navy, lwd = 0.6))
    grid::grid.text(panel$note, x = 0.5, y = 0.065,
                    gp = grid::gpar(col = navy, fontsize = 6.5, fontface = "italic"))
  }
  grid::popViewport()
}

.nature_legend_grob <- function(text = NULL) {
  pal <- .nature_palette()
  text <- text %||%
    "Legend & scale    Beta/Gamma effects: blue = negative, red = positive     Clade colors: A-E     Error bars: 95% credible interval     Dashed line: reference/null expectation     Computed by HmscEcoEvo"
  grid::grobTree(
    grid::roundrectGrob(r = grid::unit(0.01, "npc"),
                        gp = grid::gpar(fill = "white", col = pal$navy, lwd = 0.9)),
    grid::textGrob(text, x = grid::unit(0.02, "npc"), y = grid::unit(0.56, "npc"),
                   just = c("left", "center"),
                   gp = grid::gpar(col = pal$ink, fontsize = 7))
  )
}

.nature_key_grob <- function(lines) {
  pal <- .nature_palette()
  y <- seq(0.82, 0.22, length.out = length(lines))
  children <- list(
    grid::roundrectGrob(r = grid::unit(0.012, "npc"),
                        gp = grid::gpar(fill = "white", col = pal$navy, lwd = 0.9)),
    grid::textGrob("Key takeaways", x = 0.5, y = 0.93,
                   gp = grid::gpar(col = pal$navy, fontsize = 8, fontface = "bold"))
  )
  for (i in seq_along(lines)) {
    wrapped <- paste(strwrap(lines[i], width = 20), collapse = "\n")
    children[[length(children) + 1L]] <- grid::pointsGrob(
      x = 0.14, y = y[i], pch = 21, size = grid::unit(4, "mm"),
      gp = grid::gpar(fill = c(pal$blue, pal$green, pal$orange, pal$purple, pal$red)[(i - 1) %% 5 + 1],
                      col = "white", lwd = 0.5)
    )
    children[[length(children) + 1L]] <- grid::textGrob(
      wrapped, x = 0.24, y = y[i], just = c("left", "center"),
      gp = grid::gpar(col = pal$navy, fontsize = 5.3, fontface = "italic",
                      lineheight = 0.9)
    )
  }
  do.call(grid::grobTree, children)
}

.nature_page <- function(title, subtitle, panels, legend = NULL,
                         width_mm = 300, height_mm = 210) {
  structure(list(title = title, subtitle = subtitle, panels = panels,
                 legend = legend %||% .nature_legend_grob(),
                 width_mm = width_mm, height_mm = height_mm),
            class = "hmsc_ecoevo_nature_page")
}

#' Print a Nature-style HmscEcoEvo figure page
#'
#' @param x A Nature-style page object returned by a diagnostic plotting
#'   function.
#' @param ... Additional arguments, currently ignored.
#' @return Invisibly returns `x`.
#' @export
print.hmsc_ecoevo_nature_page <- function(x, ...) {
  pal <- .nature_palette()
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
  grid::grid.text(x$title, x = 0.5, y = 0.982,
                  gp = grid::gpar(col = "black", fontsize = 19,
                                  fontface = "bold"))
  grid::grid.text(x$subtitle, x = 0.5, y = 0.955,
                  gp = grid::gpar(col = pal$ink, fontsize = 10,
                                  fontface = "italic"))
  for (panel in x$panels) .draw_nature_panel(panel, pal$navy)
  grid::pushViewport(grid::viewport(x = 0.008, y = 0.012,
                                    width = 0.984, height = 0.105,
                                    just = c("left", "bottom")))
  grid::grid.draw(x$legend)
  grid::popViewport()
  invisible(x)
}

.panel <- function(label, title, plot, x, y, w, h, note = NULL) {
  list(label = label, title = title, plot = plot, x = x, y = y, w = w, h = h,
       note = note)
}

.plot_nature_beta_support <- function(metrics, axes) {
  S <- metrics$beta_summary
  S <- S[S$axis %in% axes, , drop = FALSE]
  if (nrow(S) > 0) {
    S <- stats::aggregate(
      S[, intersect(c("beta", "lwr", "upr", "support", "sd"), names(S)), drop = FALSE],
      list(axis = S$axis), mean, na.rm = TRUE
    )
  }
  axis_col <- .nature_palette()$axis
  S$axis_label <- .nature_label_axes(S$axis)
  ggplot2::ggplot(S, ggplot2::aes(beta, axis_label, colour = axis)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lwr, xmax = upr),
                           orientation = "y", width = 0.20,
                           linewidth = 0.35, na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(size = support), na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = axis_col, na.value = "grey50") +
    ggplot2::scale_size_continuous(range = c(1.5, 3.5), limits = c(0, 1)) +
    ggplot2::labs(x = "Posterior mean beta +/- 95% CrI", y = NULL) +
    .nature_theme(7) +
    ggplot2::theme(legend.position = "right")
}

.plot_nature_response_curves <- function(metrics, axes) {
  R <- metrics$response_curves
  R <- R[R$response_axis %in% axes, , drop = FALSE]
  if (nrow(R) == 0) return(ggplot2::ggplot() + .nature_theme())
  keep_species <- unique(R$species)[seq_len(min(5, length(unique(R$species))))]
  R <- R[R$species %in% keep_species, , drop = FALSE]
  R$response_axis <- factor(R$response_axis, levels = axes,
                            labels = .nature_label_axes(axes))
  ggplot2::ggplot(R, ggplot2::aes(x, response, colour = species, group = species)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = response_lwr, ymax = response_upr,
                                      fill = species),
                         alpha = 0.08, colour = NA, show.legend = FALSE,
                         na.rm = TRUE) +
    ggplot2::geom_line(linewidth = 0.45, na.rm = TRUE) +
    ggplot2::facet_wrap(~response_axis, nrow = 1, scales = "free_y") +
    ggplot2::labs(x = "Standardized environmental gradient (SD)",
                  y = "Predicted occurrence") +
    .nature_theme(6.5) +
    ggplot2::theme(legend.position = "right")
}

.plot_nature_pca <- function(evo, axes) {
  B <- evo$beta[, axes, drop = FALSE]
  cl <- .nature_clade(evo)
  pca <- stats::prcomp(B, center = TRUE, scale. = TRUE)
  P <- data.frame(species = rownames(B), clade = cl[rownames(B)], pca$x,
                  check.names = FALSE)
  V <- data.frame(axis = rownames(pca$rotation), pca$rotation,
                  check.names = FALSE)
  hull <- do.call(rbind, lapply(split(P, P$clade), function(d) {
    if (nrow(d) < 3) return(NULL)
    d[grDevices::chull(d$PC1, d$PC2), , drop = FALSE]
  }))
  scale <- max(abs(P[, c("PC1", "PC2")]), na.rm = TRUE) * 0.8
  ggplot2::ggplot(P, ggplot2::aes(PC1, PC2, fill = clade)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey70",
                        linewidth = 0.25) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey70",
                        linewidth = 0.25) +
    ggplot2::geom_polygon(data = hull,
      ggplot2::aes(group = clade), alpha = 0.13, colour = NA,
      inherit.aes = TRUE, na.rm = TRUE) +
    ggplot2::geom_point(size = 2, shape = 21, colour = "white", stroke = 0.25) +
    ggplot2::geom_segment(data = V,
      ggplot2::aes(x = 0, y = 0, xend = PC1 * scale, yend = PC2 * scale),
      inherit.aes = FALSE, arrow = ggplot2::arrow(length = grid::unit(0.09, "in")),
      linewidth = 0.35, colour = "black") +
    ggplot2::geom_text(data = V,
      ggplot2::aes(x = PC1 * scale, y = PC2 * scale,
                   label = .nature_label_axes(axis)),
      inherit.aes = FALSE, size = 2.2) +
    ggplot2::scale_fill_manual(values = .nature_palette()$clade) +
    ggplot2::labs(x = "beta-PC1", y = "beta-PC2") +
    .nature_theme(7) +
    ggplot2::theme(legend.position = "right")
}

.plot_nature_optimum <- function(metrics) {
  O <- metrics$optimum_breadth
  if (nrow(O) == 0) return(ggplot2::ggplot() + .nature_theme())
  cl <- .nature_clade(metrics$source)
  O$clade <- cl[O$species]
  ggplot2::ggplot(O, ggplot2::aes(optimum, breadth_index, colour = clade, fill = clade)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey70") +
    ggplot2::geom_hline(yintercept = stats::median(O$breadth_index, na.rm = TRUE),
                        linetype = 2, colour = "grey70") +
    ggplot2::geom_point(shape = 21, size = 2.4, colour = "white", stroke = 0.25,
                        na.rm = TRUE) +
    ggplot2::geom_text(ggplot2::aes(label = species), size = 2.0,
                       vjust = -0.8, check_overlap = TRUE) +
    ggplot2::scale_colour_manual(values = .nature_palette()$clade) +
    ggplot2::scale_fill_manual(values = .nature_palette()$clade) +
    ggplot2::labs(x = "Niche optimum (SD)", y = "Niche breadth (SD)") +
    .nature_theme(7) +
    ggplot2::theme(legend.position = "right")
}

.plot_nature_vector <- function(metrics, axes = NULL) {
  B <- metrics$beta
  if (!is.null(axes)) B <- B[, axes, drop = FALSE]
  magnitude <- sqrt(rowSums(B^2, na.rm = TRUE))
  axes <- colnames(B)
  theta <- seq(0, 2 * pi, length.out = length(axes) + 1)[-1]
  names(theta) <- axes
  seg <- do.call(rbind, lapply(rownames(B), function(sp) {
    vals <- B[sp, axes]
    data.frame(species = sp,
               x = 0, y = 0,
               xend = sum(vals * cos(theta), na.rm = TRUE),
               yend = sum(vals * sin(theta), na.rm = TRUE),
               magnitude = magnitude[sp],
               stringsAsFactors = FALSE)
  }))
  cl <- .nature_clade(metrics$source)
  seg$clade <- cl[seg$species]
  ggplot2::ggplot(seg) +
    ggplot2::geom_segment(ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                                       colour = clade, linewidth = magnitude),
                          arrow = ggplot2::arrow(length = grid::unit(0.08, "in")),
                          alpha = 0.85) +
    ggplot2::geom_text(ggplot2::aes(x = xend, y = yend, label = species),
                       size = 1.9, check_overlap = TRUE) +
    ggplot2::scale_colour_manual(values = .nature_palette()$clade) +
    ggplot2::scale_linewidth_continuous(range = c(0.25, 1.0), guide = "none") +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "Direction", y = "Direction") +
    .nature_theme(7) +
    ggplot2::theme(legend.position = "right")
}

.branch_probability_table <- function(branch) {
  if (nrow(branch) == 0) return(data.frame(branch = numeric(), P_shift = numeric()))
  z <- branch
  z$P_shift_plot <- z$P_shift
  if (all(is.na(z$P_shift_plot)) && "branch_contrast_outlier_score" %in% names(z)) {
    mx <- max(z$branch_contrast_outlier_score, na.rm = TRUE)
    if (is.finite(mx) && mx > 0) z$P_shift_plot <- pmin(1, z$branch_contrast_outlier_score / mx)
  }
  if (all(is.na(z$P_shift_plot))) z$P_shift_plot <- 0
  out <- stats::aggregate(P_shift_plot ~ branch, z, max, na.rm = TRUE)
  names(out)[2] <- "P_shift"
  out
}

.plot_shift_count_posterior <- function(m) {
  P <- m$shift_count_posterior
  if (is.null(P) || nrow(P) == 0) {
    return(ggplot2::ggplot() + .nature_theme(6.5) +
      ggplot2::labs(x = "Number of shifts", y = "Posterior density"))
  }
  ci <- stats::quantile(P$number_of_niche_shifts, c(0.025, 0.5, 0.975),
                        na.rm = TRUE)
  ggplot2::ggplot(P, ggplot2::aes(number_of_niche_shifts)) +
    ggplot2::geom_density(fill = .nature_palette()$light_blue,
                          colour = .nature_palette()$blue,
                          linewidth = 0.45, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = ci[2], colour = .nature_palette()$blue,
                        linetype = 2, linewidth = 0.35) +
    ggplot2::annotate("text", x = ci[2], y = Inf,
                      label = paste0("median=", round(ci[2], 1)),
                      vjust = 1.5, hjust = -0.05, size = 2.0,
                      colour = .nature_palette()$navy) +
    ggplot2::labs(x = "Number of shifts", y = "Density") +
    .nature_theme(6.2)
}

.plot_evo_shift_rate_tree <- function(evo, m) {
  tree <- .tree_segments(evo$phylo, species_order = rev(evo$species),
                         clade = .nature_clade(evo))
  if (nrow(tree$edges) == 0) return(ggplot2::ggplot() + .nature_theme())
  br <- .branch_probability_table(m$branch_shifts)
  mag <- stats::aggregate(abs(delta_theta) ~ branch, m$branch_shifts,
                          max, na.rm = TRUE)
  names(mag)[2] <- "shift_magnitude"
  br <- merge(br, mag, by = "branch", all.x = TRUE)
  rate <- m$branch_rate_shifts
  if (is.null(rate) || nrow(rate) == 0) {
    rate <- data.frame(branch = tree$edges$branch, rate_ratio = NA_real_,
                       P_rate_shift = NA_real_)
  }
  edges <- merge(tree$edges, br, by = "branch", all.x = TRUE)
  edges <- merge(edges, rate[, intersect(c("branch", "rate_ratio", "P_rate_shift"),
                                         names(rate)), drop = FALSE],
                 by = "branch", all.x = TRUE)
  edges$rate_ratio[is.na(edges$rate_ratio)] <- 1
  edges$P_rate_shift[is.na(edges$P_rate_shift)] <- 0
  node_points <- data.frame(node = edges$child, x = edges$xend, y = edges$yend,
                            shift_magnitude = edges$shift_magnitude,
                            P_shift = edges$P_shift, stringsAsFactors = FALSE)
  root <- setdiff(edges$parent, edges$child)[1]
  root_row <- tree$nodes[tree$nodes$node == root, , drop = FALSE]
  if (nrow(root_row) == 1) {
    node_points <- rbind(data.frame(node = root, x = root_row$x, y = root_row$y,
      shift_magnitude = max(node_points$shift_magnitude, na.rm = TRUE),
      P_shift = mean(node_points$P_shift, na.rm = TRUE),
      stringsAsFactors = FALSE), node_points)
  }
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = edges,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend, colour = clade,
                   linewidth = rate_ratio, alpha = P_rate_shift),
      na.rm = TRUE) +
    ggplot2::geom_segment(data = tree$verts,
      ggplot2::aes(x = x, xend = x, y = y, yend = yend, colour = clade),
      linewidth = 0.35, na.rm = TRUE) +
    ggplot2::geom_point(data = node_points,
      ggplot2::aes(x = x, y = y, size = shift_magnitude, fill = shift_magnitude),
      shape = 21, colour = "black", stroke = 0.22, na.rm = TRUE) +
    ggplot2::geom_text(data = tree$tips,
      ggplot2::aes(x = x + 0.025 * tree$max_x, y = y, label = species),
      hjust = 0, size = 1.8, colour = .nature_palette()$navy) +
    ggplot2::scale_colour_manual(values = .nature_palette()$clade,
                                 na.value = .nature_palette()$grey) +
    ggplot2::guides(colour = "none") +
    ggplot2::scale_linewidth_continuous(range = c(0.2, 1.25),
                                        name = "Rate ratio") +
    ggplot2::scale_alpha_continuous(range = c(0.35, 1), guide = "none") +
    ggplot2::scale_size_continuous(range = c(1.2, 5.2),
                                   name = "Shift magnitude") +
    ggplot2::scale_fill_gradientn(colours = c("#FFF7BC", "#FEC44F", "#F26B21", "#D9272E"),
                                  name = "Shift magnitude",
                                  na.value = "grey85", guide = "none") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = "Time from root (Ma)", y = NULL) +
    .nature_theme(6.4) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   legend.position = "right",
                   legend.text = ggplot2::element_text(size = 4.8),
                   legend.title = ggplot2::element_text(size = 5.2),
                   plot.margin = ggplot2::margin(2, 12, 2, 2))
}

.rate_diagnostics_grob <- function(m) {
  R <- m$rate_metrics
  B <- m$branch_rate_shifts
  clade_rate <- R[R$clade != "global", , drop = FALSE]
  clade_mean <- stats::aggregate(rate_ratio ~ clade, clade_rate, mean, na.rm = TRUE)
  p1 <- ggplot2::ggplot(clade_mean,
    ggplot2::aes(rate_ratio, stats::reorder(clade, rate_ratio), colour = clade)) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "grey60") +
    ggplot2::geom_point(size = 2.0, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = .nature_palette()$clade, na.value = "grey50") +
    ggplot2::labs(x = "Clade/background rate", y = NULL) + .nature_theme(5.8) +
    ggplot2::theme(legend.position = "none")
  p2 <- ggplot2::ggplot(B, ggplot2::aes(branch, P_rate_shift)) +
    ggplot2::geom_col(fill = .nature_palette()$blue, width = 0.75, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 2,
                        colour = .nature_palette()$purple) +
    ggplot2::labs(x = "Branch", y = "P rate shift") + .nature_theme(5.8)
  p3 <- ggplot2::ggplot(clade_rate,
    ggplot2::aes(axis, clade, fill = axis_specific_rate)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.12) +
    ggplot2::scale_x_discrete(labels = .nature_label_axes) +
    ggplot2::scale_fill_gradientn(colours = c("#1F5FBF", "white", "#FEC44F", "#D9272E"),
                                  name = "Rate") +
    ggplot2::labs(x = NULL, y = NULL) + .nature_theme(5.8)
  .ggplot_grid_grob(list(p1, p2, p3), data.frame(
    x = c(0.00, 0.00, 0.52), y = c(0.50, 0.00, 0.00),
    w = c(0.48, 0.48, 0.48), h = c(0.48, 0.46, 0.96)
  ))
}

.diversification_grob <- function(m) {
  V <- m$variance_metrics
  vmean <- if (nrow(V) > 0) stats::aggregate(variance_ratio ~ clade, V, mean, na.rm = TRUE) else data.frame()
  p1 <- ggplot2::ggplot(vmean,
    ggplot2::aes(variance_ratio, stats::reorder(clade, variance_ratio), colour = clade)) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "grey60") +
    ggplot2::geom_point(size = 2.0, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = .nature_palette()$clade, na.value = "grey50") +
    ggplot2::labs(x = "Variance ratio", y = NULL) + .nature_theme(5.8) +
    ggplot2::theme(legend.position = "none")
  p2 <- ggplot2::ggplot(m$dtt, ggplot2::aes(time, disparity, colour = axis))
  if (all(c("bm_lwr", "bm_upr") %in% names(m$dtt))) {
    p2 <- p2 + ggplot2::geom_ribbon(
      ggplot2::aes(x = time, ymin = bm_lwr, ymax = bm_upr),
      fill = .nature_palette()$light_blue, alpha = 0.45,
      colour = NA, inherit.aes = FALSE, na.rm = TRUE)
  }
  p2 <- p2 +
    ggplot2::geom_line(linewidth = 0.45, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = .nature_palette()$axis, na.value = "grey50") +
    ggplot2::labs(x = "Time from root (Ma)", y = "DTT") + .nature_theme(5.8) +
    ggplot2::theme(legend.position = "none")
  eb <- data.frame()
  if (nrow(m$early_burst) > 0) {
    eb <- rbind(eb, data.frame(axis = m$early_burst$axis, metric = "EB alpha",
                               value = m$early_burst$early_burst_parameter))
  }
  if (nrow(m$mdi) > 0) {
    eb <- rbind(eb, data.frame(axis = m$mdi$axis, metric = "MDI",
                               value = m$mdi$MDI))
  }
  p3 <- ggplot2::ggplot(eb, ggplot2::aes(value, axis, colour = metric)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    ggplot2::geom_point(size = 1.7, na.rm = TRUE) +
    ggplot2::facet_wrap(~metric, nrow = 1, scales = "free_x") +
    ggplot2::scale_y_discrete(labels = .nature_label_axes) +
    ggplot2::labs(x = "Estimate", y = NULL) + .nature_theme(5.7) +
    ggplot2::theme(legend.position = "none")
  .ggplot_grid_grob(list(p1, p2, p3), data.frame(
    x = c(0.00, 0.42, 0.00), y = c(0.50, 0.50, 0.00),
    w = c(0.38, 0.58, 0.96), h = c(0.46, 0.46, 0.45)
  ))
}

.plot_regime_tree <- function(evo, regimes) {
  tree <- .tree_segments(evo$phylo, species_order = rev(evo$species),
                         clade = .nature_clade(evo))
  if (is.null(regimes) || nrow(regimes) == 0) {
    regimes <- data.frame(species = evo$species, peak = evo$regime[evo$species])
  }
  tips <- merge(tree$tips, regimes[, intersect(c("species", "peak"), names(regimes)),
                                   drop = FALSE],
                by = "species", all.x = TRUE)
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = tree$edges,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend, colour = clade),
      linewidth = 0.38, na.rm = TRUE) +
    ggplot2::geom_segment(data = tree$verts,
      ggplot2::aes(x = x, xend = x, y = y, yend = yend, colour = clade),
      linewidth = 0.38, na.rm = TRUE) +
    ggplot2::geom_point(data = tips, ggplot2::aes(x = x, y = y, fill = peak),
                        shape = 21, size = 2.0, colour = "black",
                        stroke = 0.22, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = .nature_palette()$clade,
                                 na.value = .nature_palette()$grey) +
    ggplot2::guides(colour = "none") +
    ggplot2::labs(x = "Time from root (Ma)", y = NULL, fill = "Regime") +
    .nature_theme(5.8) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   legend.position = "bottom",
                   legend.text = ggplot2::element_text(size = 4.2),
                   legend.title = ggplot2::element_text(size = 4.8),
                   legend.key.height = grid::unit(2.4, "mm"),
                   legend.key.width = grid::unit(3.2, "mm"))
}

.integration_network_plot <- function(net) {
  axes <- unique(c(net$axis1, net$axis2))
  theta <- seq(0, 2 * pi, length.out = length(axes) + 1)[-1]
  nodes <- data.frame(axis = axes, x = cos(theta), y = sin(theta),
                      stringsAsFactors = FALSE)
  ed <- net[net$axis1 != net$axis2 & is.finite(net$correlation), , drop = FALSE]
  ed <- ed[!duplicated(t(apply(ed[, c("axis1", "axis2")], 1, sort))), , drop = FALSE]
  ed <- merge(ed, nodes, by.x = "axis1", by.y = "axis", all.x = TRUE)
  ed <- merge(ed, nodes, by.x = "axis2", by.y = "axis", all.x = TRUE,
              suffixes = c("1", "2"))
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = ed,
      ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2,
                   linewidth = abs(correlation), colour = correlation),
      alpha = 0.85, na.rm = TRUE) +
    ggplot2::geom_point(data = nodes, ggplot2::aes(x, y),
                        size = 4, shape = 21, fill = "white",
                        colour = .nature_palette()$navy) +
    ggplot2::geom_text(data = nodes, ggplot2::aes(x, y, label = .nature_label_axes(axis)),
                       size = 2.0) +
    ggplot2::scale_colour_gradient2(low = "#1F5FBF", mid = "grey80",
                                    high = "#D9272E", midpoint = 0) +
    ggplot2::scale_linewidth_continuous(range = c(0.15, 1.2), guide = "none") +
    ggplot2::coord_equal() +
    ggplot2::labs(x = NULL, y = NULL) + .nature_theme(5.8) +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank(),
                   legend.position = "none")
}

.convergence_grob <- function(evo, m) {
  p1 <- .plot_regime_tree(evo, m$adaptive_regimes)
  P <- m$phylo_beta_distance
  P$convergent <- "Other"
  if ("convergence_score" %in% names(P)) {
    P$convergent[!is.na(P$convergence_score) & P$convergence_score > 0] <- "Convergent"
  }
  p2 <- ggplot2::ggplot(P, ggplot2::aes(phylo_distance, beta_distance,
                                        colour = convergent)) +
    ggplot2::geom_point(size = 1.1, alpha = 0.75, na.rm = TRUE) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                         colour = "grey30", linewidth = 0.35, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = c(Convergent = .nature_palette()$red,
                                            Other = "grey70")) +
    ggplot2::labs(x = "Phylogenetic distance", y = "Beta distance") +
    .nature_theme(5.8)
  C <- .matrix_to_long(m$beta_cosine_similarity, "cosine", "species1", "species2")
  p3 <- ggplot2::ggplot(C, ggplot2::aes(species1, species2, fill = cosine)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.08) +
    ggplot2::scale_fill_gradient2(low = "#1F5FBF", mid = "white",
                                  high = "#D9272E", midpoint = 0,
                                  name = "Cosine") +
    ggplot2::labs(x = NULL, y = NULL) + .nature_theme(5.2) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,
                                                       size = 3.8),
                   axis.text.y = ggplot2::element_text(size = 3.8))
  p4 <- .integration_network_plot(m$evolutionary_integration)
  .ggplot_grid_grob(list(p1, p2, p3, p4), data.frame(
    x = c(0.00, 0.52, 0.00, 0.52), y = c(0.52, 0.52, 0.00, 0.00),
    w = c(0.48, 0.46, 0.48, 0.46), h = c(0.44, 0.44, 0.44, 0.44)
  ))
}

.plot_gamma_shift_tree <- function(evo, m) {
  tree <- .tree_segments(evo$phylo, species_order = rev(evo$species),
                         clade = .nature_clade(evo))
  if (nrow(tree$edges) == 0) return(ggplot2::ggplot() + .nature_theme())
  br <- m$gamma_branch_shift_probability
  if (is.null(br) || nrow(br) == 0) {
    shift <- m$gamma_shift_probability
    val <- mean(shift$gamma_shift_probability, na.rm = TRUE)
    br <- data.frame(branch = tree$edges$branch,
                     gamma_shift_probability = val,
                     stringsAsFactors = FALSE)
  }
  edges <- merge(tree$edges, br[, intersect(c("branch", "gamma_shift_probability"),
                                           names(br)), drop = FALSE],
                 by = "branch", all.x = TRUE)
  node_points <- data.frame(node = edges$child, x = edges$xend, y = edges$yend,
                            gamma_shift_probability = edges$gamma_shift_probability,
                            stringsAsFactors = FALSE)
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = edges,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend,
                   colour = gamma_shift_probability,
                   linewidth = gamma_shift_probability),
      na.rm = TRUE) +
    ggplot2::geom_segment(data = tree$verts,
      ggplot2::aes(x = x, xend = x, y = y, yend = yend),
      colour = "grey35", linewidth = 0.25, na.rm = TRUE) +
    ggplot2::geom_point(data = node_points,
      ggplot2::aes(x = x, y = y, size = gamma_shift_probability,
                   fill = gamma_shift_probability),
      shape = 21, colour = "black", stroke = 0.22, na.rm = TRUE) +
    ggplot2::geom_text(data = tree$tips,
      ggplot2::aes(x = x + 0.025 * tree$max_x, y = y, label = species),
      hjust = 0, size = 1.8, colour = .nature_palette()$navy) +
    ggplot2::scale_colour_gradientn(colours = c("#2C7BB6", "#ABD9E9", "#FFFFBF",
                                                "#FDAE61", "#D7191C"),
                                    limits = c(0, 1),
                                    name = "Gamma-shift P") +
    ggplot2::scale_fill_gradientn(colours = c("#2C7BB6", "#ABD9E9", "#FFFFBF",
                                              "#FDAE61", "#D7191C"),
                                  limits = c(0, 1),
                                  name = "Gamma-shift P") +
    ggplot2::scale_size_continuous(range = c(1.1, 4.2),
                                   limits = c(0, 1),
                                   name = "Gamma-shift P") +
    ggplot2::scale_linewidth_continuous(range = c(0.25, 1.2), guide = "none") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = "Time from root (Ma)", y = NULL) +
    .nature_theme(6.2) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   legend.position = "right",
                   plot.margin = ggplot2::margin(2, 12, 2, 2))
}

.gamma_diagnostics_grob <- function(m) {
  sf <- m$gamma_sign_flip_index
  gv <- m$gamma_variance_among_clades
  sf_sum <- if (nrow(sf) > 0) {
    stats::aggregate(gamma_sign_flip_index ~ trait, sf, mean, na.rm = TRUE)
  } else {
    data.frame(trait = character(), gamma_sign_flip_index = numeric())
  }
  gv_sum <- if (nrow(gv) > 0) {
    stats::aggregate(gamma_variance_among_clades ~ trait, gv, mean, na.rm = TRUE)
  } else {
    data.frame(trait = character(), gamma_variance_among_clades = numeric())
  }
  p1 <- ggplot2::ggplot(sf_sum,
    ggplot2::aes(gamma_sign_flip_index,
                 stats::reorder(trait, gamma_sign_flip_index))) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    ggplot2::geom_point(colour = .nature_palette()$purple, size = 1.9,
                        na.rm = TRUE) +
    ggplot2::labs(x = "Sign-flip index", y = NULL) + .nature_theme(5.8)
  p2 <- ggplot2::ggplot(gv_sum,
    ggplot2::aes(gamma_variance_among_clades,
                 stats::reorder(trait, gamma_variance_among_clades))) +
    ggplot2::geom_point(colour = .nature_palette()$purple, size = 1.9,
                        na.rm = TRUE) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = gamma_variance_among_clades,
                                       yend = trait),
                          colour = .nature_palette()$purple, linewidth = 0.35,
                          na.rm = TRUE) +
    ggplot2::labs(x = "Gamma variance among clades", y = NULL) +
    .nature_theme(5.8)
  .ggplot_grid_grob(list(p1, p2), data.frame(
    x = c(0.00, 0.50), y = c(0.00, 0.00),
    w = c(0.47, 0.47), h = c(0.96, 0.96)
  ))
}

.turnover_from_group_gamma <- function(group_gamma, group_col) {
  if (is.null(group_gamma) || nrow(group_gamma) == 0 || !group_col %in% names(group_gamma)) {
    return(list(index = data.frame(), matrix = data.frame()))
  }
  d <- group_gamma
  names(d)[names(d) == group_col] <- "clade"
  .gamma_turnover(d)
}

.turnover_summary <- function(idx, label) {
  if (is.null(idx) || nrow(idx) == 0) {
    return(data.frame(group = character(), turnover = numeric(), type = character()))
  }
  g1 <- stats::aggregate(turnover ~ group1, idx, mean, na.rm = TRUE)
  g2 <- stats::aggregate(turnover ~ group2, idx, mean, na.rm = TRUE)
  names(g1)[1] <- names(g2)[1] <- "group"
  out <- stats::aggregate(turnover ~ group, rbind(g1, g2), mean, na.rm = TRUE)
  out$type <- label
  out
}

.turnover_diagnostics_grob <- function(m) {
  clade_idx <- m$trait_function_turnover_index
  if (is.null(clade_idx) || nrow(clade_idx) == 0) {
    clade_idx <- data.frame(group1 = character(), group2 = character(),
                            turnover = numeric())
  }
  regime_turn <- .turnover_from_group_gamma(m$regime_specific_gamma, "regime")
  summary <- rbind(.turnover_summary(clade_idx, "clade"),
                   .turnover_summary(regime_turn$index, "regime"))
  p1 <- ggplot2::ggplot(summary,
    ggplot2::aes(turnover, stats::reorder(group, turnover), colour = type)) +
    ggplot2::geom_point(size = 1.9, na.rm = TRUE) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = turnover * 0.90,
                                        xmax = turnover * 1.10),
                           orientation = "y", width = 0.16,
                           linewidth = 0.3, na.rm = TRUE) +
    ggplot2::labs(x = "Turnover index", y = NULL) + .nature_theme(5.8)
  p2 <- ggplot2::ggplot(clade_idx,
    ggplot2::aes(group1, group2, fill = turnover)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.16) +
    ggplot2::scale_fill_gradient2(low = "#1F5FBF", mid = "white",
                                  high = "#D9272E", midpoint = 1,
                                  name = "TI") +
    ggplot2::labs(x = NULL, y = NULL) + .nature_theme(5.8)
  .ggplot_grid_grob(list(p1, p2), data.frame(
    x = c(0.00, 0.50), y = c(0.00, 0.00),
    w = c(0.47, 0.47), h = c(0.96, 0.96)
  ))
}

.plot_niche_summary_nature <- function(x, ...) {
  evo <- if (inherits(x, "hmsc_niche_metrics")) x$source else as_hmsc_ecoevo(x)
  axes <- .nature_analysis_axes(evo)
  q_axes <- paste0(axes, "_sq")
  q_map <- data.frame(linear = axes[q_axes %in% colnames(evo$beta)],
                      quadratic = q_axes[q_axes %in% colnames(evo$beta)],
                      stringsAsFactors = FALSE)
  calc_axes <- unique(c(axes, q_map$quadratic))
  metrics <- if (inherits(x, "hmsc_niche_metrics")) x else
    calc_niche_metrics(evo, axes = calc_axes, quadratic = q_map, ...)
  B <- evo$beta[, axes, drop = FALSE]
  key <- .nature_key_grob(c(
    "Optima and breadth reveal trade-offs.",
    "Beta direction shows environmental response axes.",
    "Clade hulls expose niche overlap and divergence.",
    "Uncertainty is carried from posterior draws."
  ))
  panels <- list(
    .panel("A", "Phylogeny and species-by-environment beta",
           .tree_heatmap_plot(evo, B, fill_title = "beta"),
           0.008, 0.555, 0.305, 0.370,
           "Red = positive response; blue = negative response."),
    .panel("B", "Beta posterior support and uncertainty",
           .plot_nature_beta_support(metrics, axes),
           0.318, 0.555, 0.290, 0.370,
           "Posterior means with 95% credible intervals and support."),
    .panel("C", "Niche response curves along environmental gradients",
           .plot_nature_response_curves(metrics, axes),
           0.613, 0.555, 0.379, 0.370,
           "Curves show species-specific optima and breadths."),
    .panel("D", "Niche optimum vs. niche breadth",
           .plot_nature_optimum(metrics),
           0.008, 0.145, 0.270, 0.390,
           "Upper right: higher optimum and broader niche."),
    .panel("E", "Multivariate niche space (beta-PCA)",
           .plot_nature_pca(evo, axes),
           0.283, 0.145, 0.270, 0.390,
           "Polygons are clade hulls; arrows are environmental gradients."),
    .panel("F", "Niche vector summary and uncertainty",
           .plot_nature_vector(metrics, axes),
           0.558, 0.145, 0.310, 0.390,
           "Length = response strength; direction = beta association."),
    .panel("", "Key takeaways", key, 0.875, 0.145, 0.117, 0.390)
  )
  .nature_page("Figure 1. Niche summary",
               "Extracting ecological niches from HMSC beta", panels)
}

.plot_phylo_signal_nature <- function(x, ...) {
  m <- if (inherits(x, "hmsc_phylo_signal_metrics")) x else calc_phylo_signal_metrics(x, ...)
  evo <- m$source
  axes <- .nature_analysis_axes(evo)
  G <- m$global_signal
  C <- m$correlogram
  S <- m$scale_dependent_signal
  L <- m$local_signal
  N <- m$node_signal
  R <- m$residual_phylogenetic_map
  CCI <- m$clade_conservatism
  metric_order <- c("Hmsc_rho", "Pagel_lambda", "Blomberg_K",
                    "Abouheif_Cmean", "Moran_I")
  metric_labels <- c(Hmsc_rho = "Hmsc rho",
                     Pagel_lambda = "Pagel's lambda",
                     Blomberg_K = "Blomberg's K",
                     Abouheif_Cmean = "Abouheif Cmean",
                     Moran_I = "Moran's I")
  G <- G[G$axis %in% axes & G$metric %in% metric_order, , drop = FALSE]
  G$axis_label <- factor(.nature_label_axes(G$axis),
                         levels = rev(.nature_label_axes(axes)))
  G$metric_label <- factor(metric_labels[G$metric],
                           levels = metric_labels[metric_order])
  pA <- ggplot2::ggplot(G, ggplot2::aes(estimate, axis_label, colour = axis)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lwr, xmax = upr),
                           orientation = "y", width = 0.16, na.rm = TRUE) +
    ggplot2::geom_point(size = 1.8, na.rm = TRUE) +
    ggplot2::facet_grid(metric_label ~ ., scales = "free_y", space = "free_y",
                        switch = "y") +
    ggplot2::scale_colour_manual(values = .nature_palette()$axis, na.value = "grey50") +
    ggplot2::labs(x = "Estimate and 95% CI", y = NULL) + .nature_theme(6.1) +
    ggplot2::theme(strip.placement = "outside",
                   strip.text.y.left = ggplot2::element_text(angle = 0, size = 5.4),
                   panel.spacing.y = grid::unit(0.08, "lines"),
                   legend.position = "right")
  pB <- ggplot2::ggplot(C, ggplot2::aes(phylo_distance, similarity, colour = axis))
  if (all(c("lwr", "upr") %in% names(C))) {
    pB <- pB + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lwr, ymax = upr, fill = axis),
      alpha = 0.10, colour = NA, show.legend = FALSE, na.rm = TRUE)
  }
  if (!"signif" %in% names(C)) C$signif <- if ("p_value" %in% names(C)) .p_signif_label(C$p_value) else ""
  sig_df <- C[!is.na(C$signif) & C$signif != "", , drop = FALSE]
  pB <- pB +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey65") +
    ggplot2::geom_line(linewidth = 0.55, na.rm = TRUE) +
    ggplot2::geom_point(size = 1.7, na.rm = TRUE) +
    ggplot2::geom_text(data = sig_df, ggplot2::aes(label = signif),
                       vjust = -0.55, size = 2.1, show.legend = FALSE,
                       na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = .nature_palette()$axis, na.value = "grey50") +
    ggplot2::scale_fill_manual(values = .nature_palette()$axis, na.value = "grey50") +
    ggplot2::labs(x = "Phylogenetic distance", y = "Moran's I") + .nature_theme(7)
  scale_df <- data.frame(section = "Scale-dependent rho",
                         label = S$scale, value = S$signal, group = S$axis,
                         stringsAsFactors = FALSE)
  cci_df <- if (nrow(CCI) > 0) {
    data.frame(section = "Clade conservatism index",
               label = CCI$clade, value = CCI$CCI, group = CCI$clade,
               stringsAsFactors = FALSE)
  } else data.frame(section = character(), label = character(),
                    value = numeric(), group = character())
  scale_cci <- rbind(scale_df, cci_df)
  pC <- ggplot2::ggplot(scale_cci, ggplot2::aes(label, value, fill = group)) +
    ggplot2::geom_col(position = "dodge", width = 0.7, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
    ggplot2::facet_wrap(~section, scales = "free_x", nrow = 1) +
    ggplot2::scale_fill_manual(values = c(.nature_palette()$axis, .nature_palette()$clade),
                               na.value = "grey50") +
    ggplot2::labs(x = NULL, y = "Signal") + .nature_theme(6.5) +
    ggplot2::theme(legend.position = "none")
  pD <- .plot_phylo_tree(evo,
    tip_metric = data.frame(species = L$species[L$axis == axes[1]],
                            value = L$local_moran_i[L$axis == axes[1]]),
    title = NULL)
  node_df <- if (nrow(N) > 0) {
    a <- stats::aggregate(node_signal ~ node, N, mean, na.rm = TRUE)
    b <- stats::aggregate(n_species ~ node, N, max, na.rm = TRUE)
    out <- merge(a, b, by = "node", all.x = TRUE)
    names(out)[names(out) == "node_signal"] <- "value"
    out
  } else data.frame(node = integer(), value = numeric(), n_species = numeric())
  pE <- .plot_phylo_tree(evo, node_metric = node_df, show_legend = TRUE,
                         fill_name = "Node signal",
                         size_name = "Descendant tips")
  residual_lim <- max(0.45, max(abs(R$residual_beta), na.rm = TRUE))
  pF <- ggplot2::ggplot(R, ggplot2::aes(axis, species, fill = residual_beta)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.18) +
    ggplot2::scale_fill_gradient2(low = "#1F5FBF", mid = "white",
                                  high = "#D9272E", midpoint = 0,
                                  limits = c(-residual_lim, residual_lim)) +
    ggplot2::scale_x_discrete(labels = .nature_label_axes) +
    ggplot2::labs(x = NULL, y = NULL) + .nature_theme(7)
  key <- .nature_key_grob(c("Global metrics detect overall signal.",
                            "Correlograms locate evolutionary distance.",
                            "Local maps identify tips and lineages.",
                            "Residuals show unexplained structure."))
  panels <- list(
    .panel("A", "Global phylogenetic signal", pA, 0.008, 0.555, 0.330, 0.370),
    .panel("B", "Phylogenetic correlogram", pB, 0.343, 0.555, 0.300, 0.370),
    .panel("C", "Scale-dependent signal and clade conservatism", pC, 0.648, 0.555, 0.344, 0.370),
    .panel("D", "Tip-level local phylogenetic signal", pD, 0.008, 0.145, 0.320, 0.390),
    .panel("E", "Internal node-level signal", pE, 0.333, 0.145, 0.300, 0.390),
    .panel("F", "Phylogenetic residual signal summary", pF, 0.638, 0.145, 0.235, 0.390),
    .panel("", "Key takeaways", key, 0.880, 0.145, 0.112, 0.390)
  )
  .nature_page("Figure 2. Phylogenetic signal localization",
               "Where is phylogenetic structure in HMSC niche responses?",
               panels)
}

.plot_evo_transition_nature <- function(x, ...) {
  m <- if (inherits(x, "hmsc_evo_transition_metrics")) x else calc_evo_transition_metrics(x, ...)
  evo <- m$source
  B <- m$branch_shifts
  branch_prob <- .branch_probability_table(B)
  pA1 <- ggplot2::ggplot(branch_prob,
    ggplot2::aes(stats::reorder(branch, -P_shift), P_shift)) +
    ggplot2::geom_col(fill = .nature_palette()$blue, width = 0.75, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 2, colour = .nature_palette()$purple) +
    ggplot2::labs(x = "Branch index (ranked)", y = "P shift") + .nature_theme(6.2) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
  pA2 <- .plot_shift_count_posterior(m)
  pA <- .ggplot_grid_grob(list(pA1, pA2), data.frame(
    x = c(0.00, 0.64), y = c(0.00, 0.00),
    w = c(0.60, 0.34), h = c(0.96, 0.96)
  ))
  pB <- ggplot2::ggplot(B, ggplot2::aes(abs(delta_theta), branch, colour = axis)) +
    ggplot2::geom_point(ggplot2::aes(size = P_shift), na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = .nature_palette()$axis, na.value = "grey50") +
    ggplot2::labs(x = "Shift magnitude |Delta theta|", y = "Branch") + .nature_theme(7)
  pC <- .plot_evo_shift_rate_tree(evo, m)
  pD <- .rate_diagnostics_grob(m)
  pE <- .diversification_grob(m)
  pF <- .convergence_grob(evo, m)
  panels <- list(
    .panel("A", "Branch-level niche shift probability and total shifts", pA, 0.008, 0.555, 0.300, 0.370),
    .panel("B", "Shift magnitude and axis-specific support", pB, 0.313, 0.555, 0.360, 0.370),
    .panel("C", "Phylogeny: optimal shifts and evolutionary rates", pC, 0.678, 0.555, 0.314, 0.370),
    .panel("D", "Rate diagnostics", pD, 0.008, 0.145, 0.300, 0.390),
    .panel("E", "Diversification diagnostics", pE, 0.313, 0.145, 0.330, 0.390),
    .panel("F", "Convergence and multivariate evolution", pF, 0.648, 0.145, 0.344, 0.390)
  )
  .nature_page("Figure 3. Evolutionary transition and diversification",
               "Shift, rate, expansion, radiation, convergence, and multivariate evolution",
               panels)
}

.plot_trait_mediation_nature <- function(x, ...) {
  m <- if (inherits(x, "hmsc_trait_mediation_metrics")) x else calc_trait_mediation_metrics(x, ...)
  evo <- m$source
  G <- m$gamma
  pA <- ggplot2::ggplot(G, ggplot2::aes(axis, trait, fill = gamma)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.18) +
    ggplot2::scale_fill_gradient2(low = "#1F5FBF", mid = "white",
                                  high = "#D9272E", midpoint = 0) +
    ggplot2::scale_x_discrete(labels = .nature_label_axes) +
    ggplot2::labs(x = "Environmental beta gradients", y = "Traits") + .nature_theme(7)
  pB <- ggplot2::ggplot(m$signal_decomposition,
    ggplot2::aes(component, value, fill = component)) +
    ggplot2::geom_col(width = 0.65, na.rm = TRUE) +
    ggplot2::facet_wrap(~axis, nrow = 1) +
    ggplot2::labs(x = NULL, y = "Phylogenetic signal share") + .nature_theme(6.5) +
    ggplot2::theme(legend.position = "none")
  risk <- merge(m$missing_trait_risk,
                m$trait_explained_r2[, c("axis", "r2")], by = "axis",
                all.x = TRUE)
  pC <- ggplot2::ggplot(risk,
    ggplot2::aes(r2, missing_trait_risk_index, colour = axis, label = axis)) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = 2, colour = "grey60") +
    ggplot2::geom_hline(yintercept = 0.3, linetype = 2, colour = "grey60") +
    ggplot2::geom_point(size = 3, na.rm = TRUE) +
    ggplot2::geom_text(vjust = -0.8, size = 2.2) +
    ggplot2::scale_colour_manual(values = .nature_palette()$axis, na.value = "grey50") +
    ggplot2::labs(x = "Trait-explained variation (R2)",
                  y = "Residual phylogenetic signal") + .nature_theme(7)
  pD <- ggplot2::ggplot(m$trait_phylogenetic_redundancy,
    ggplot2::aes(redundancy, moran_i, label = trait)) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = 2, colour = "grey60") +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 2, colour = "grey60") +
    ggplot2::geom_point(size = 2.6, colour = .nature_palette()$green, na.rm = TRUE) +
    ggplot2::geom_text(vjust = -0.8, size = 2.0, check_overlap = TRUE) +
    ggplot2::labs(x = "Trait explanatory overlap", y = "Trait phylogenetic signal") +
    .nature_theme(7)
  pE <- ggplot2::ggplot(m$trait_omission_sensitivity,
    ggplot2::aes(stats::reorder(omitted_trait, delta_r2), delta_r2, fill = axis)) +
    ggplot2::geom_col(position = "dodge", width = 0.75, na.rm = TRUE) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = .nature_palette()$axis, na.value = "grey50") +
    ggplot2::labs(x = "Trait omitted", y = "Increase in residual signal") +
    .nature_theme(7)
  pF <- .tree_heatmap_plot(evo,
    reshape2_like_matrix(m$residual_beta_tree_heatmap, "species", "axis", "residual_beta"),
    fill_title = "Residual")
  key <- .nature_key_grob(c("Measured traits explain conserved niche structure.",
                            "Hidden residual signal prioritizes new traits.",
                            "Omission tests rank high-impact traits."))
  panels <- list(
    .panel("A", "Gamma: trait-by-beta effect matrix", pA, 0.008, 0.555, 0.330, 0.370),
    .panel("B", "Phylogenetic signal decomposition", pB, 0.343, 0.555, 0.320, 0.370),
    .panel("C", "Missing-trait risk map", pC, 0.668, 0.555, 0.324, 0.370),
    .panel("D", "Trait phylogenetic redundancy", pD, 0.008, 0.145, 0.290, 0.390),
    .panel("E", "Trait omission sensitivity", pE, 0.303, 0.145, 0.290, 0.390),
    .panel("F", "Residual beta structure after full trait set", pF, 0.598, 0.145, 0.270, 0.390),
    .panel("", "Key takeaways", key, 0.875, 0.145, 0.117, 0.390)
  )
  .nature_page("Figure 4. Trait mediation and hidden traits",
               "How much phylogenetic niche structure is explained by measured traits?",
               panels)
}

reshape2_like_matrix <- function(df, row, col, value) {
  rows <- unique(df[[row]])
  cols <- unique(df[[col]])
  M <- matrix(NA_real_, length(rows), length(cols), dimnames = list(rows, cols))
  for (i in seq_len(nrow(df))) M[as.character(df[[row]][i]), as.character(df[[col]][i])] <- df[[value]][i]
  M
}

.plot_gamma_evolution_nature <- function(x, ...) {
  m <- if (inherits(x, "hmsc_gamma_evolution_metrics")) x else calc_gamma_evolution_metrics(x, ...)
  evo <- m$source
  pA <- ggplot2::ggplot(m$clade_specific_gamma,
    ggplot2::aes(axis, trait, fill = gamma)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.15) +
    ggplot2::facet_wrap(~clade, nrow = 1) +
    ggplot2::scale_fill_gradient2(low = "#1F5FBF", mid = "white",
                                  high = "#D9272E", midpoint = 0) +
    ggplot2::scale_x_discrete(labels = .nature_label_axes) +
    ggplot2::labs(x = NULL, y = "Traits") + .nature_theme(6.3) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,
                                                       vjust = 0.5, size = 4.5))
  pB <- ggplot2::ggplot(m$regime_specific_gamma,
    ggplot2::aes(axis, trait, fill = gamma)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.15) +
    ggplot2::facet_wrap(~regime, nrow = 1) +
    ggplot2::scale_fill_gradient2(low = "#1F5FBF", mid = "white",
                                  high = "#D9272E", midpoint = 0) +
    ggplot2::scale_x_discrete(labels = .nature_label_axes) +
    ggplot2::labs(x = NULL, y = "Traits") + .nature_theme(6.3) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,
                                                       vjust = 0.5, size = 4.3))
  pC <- .plot_gamma_shift_tree(evo, m)
  pD <- .gamma_diagnostics_grob(m)
  pE <- ggplot2::ggplot(m$trait_function_network,
    ggplot2::aes(axis, trait, fill = weight)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.16) +
    ggplot2::scale_fill_gradient2(low = "#1F5FBF", mid = "white",
                                  high = "#D9272E", midpoint = 0) +
    ggplot2::scale_x_discrete(labels = .nature_label_axes) +
    ggplot2::labs(x = "Niche axes", y = "Traits") + .nature_theme(7)
  pF <- .turnover_diagnostics_grob(m)
  panels <- list(
    .panel("A", "Clade-specific Gamma", pA, 0.008, 0.555, 0.330, 0.370),
    .panel("B", "Regime-specific Gamma", pB, 0.343, 0.555, 0.320, 0.370),
    .panel("C", "Phylogenetic Gamma-shift probability", pC, 0.668, 0.555, 0.324, 0.370),
    .panel("D", "Sign-flip and Gamma-variance diagnostics", pD, 0.008, 0.145, 0.330, 0.390),
    .panel("E", "Trait-function network", pE, 0.343, 0.145, 0.320, 0.390),
    .panel("F", "Trait-function turnover diagnostics", pF, 0.668, 0.145, 0.324, 0.390)
  )
  .nature_page("Figure 5. Trait-function evolution",
               "Do trait-to-niche relationships evolve across clades and regimes?",
               panels)
}

.plot_population_evolution_nature <- function(x, ...) {
  m <- if (inherits(x, "hmsc_population_evolution_metrics")) x else calc_population_evolution_metrics(x, ...)
  P <- m$population_beta
  pA <- ggplot2::ggplot(P, ggplot2::aes(beta, genetic_value, colour = species)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::geom_point(alpha = 0.65, size = 1.5, na.rm = TRUE) +
    ggplot2::labs(x = "Species-level beta", y = "Population-level beta") +
    .nature_theme(7) + ggplot2::theme(legend.position = "none")
  pB <- ggplot2::ggplot(m$within_species_niche_divergence,
    ggplot2::aes(within_species_niche_divergence, species, fill = species)) +
    ggplot2::geom_col(width = 0.65, na.rm = TRUE) +
    ggplot2::labs(x = "Niche beta divergence", y = NULL) +
    .nature_theme(7) + ggplot2::theme(legend.position = "none")
  V <- m$variance_partition
  component_cols <- intersect(c("species_mean_contribution",
                                "local_adaptation_contribution",
                                "plasticity_contribution",
                                "residual_contribution"), names(V))
  component_labels <- c(species_mean_contribution = "Species mean",
                        local_adaptation_contribution = "Local adaptation",
                        plasticity_contribution = "Plasticity",
                        residual_contribution = "Residual")
  vp_base <- stats::aggregate(V[, component_cols, drop = FALSE],
                              list(species = V$species), mean, na.rm = TRUE)
  all_row <- as.data.frame(as.list(colMeans(vp_base[, component_cols, drop = FALSE],
                                            na.rm = TRUE)))
  all_row$species <- "All spp."
  vp_base <- rbind(vp_base, all_row[, names(vp_base), drop = FALSE])
  vp <- stats::reshape(vp_base,
                       varying = component_cols,
                       v.names = "value", timevar = "component",
                       times = unname(component_labels[component_cols]),
                       direction = "long")
  vp$component <- factor(vp$component,
                         levels = c("Species mean", "Local adaptation",
                                    "Plasticity", "Residual"))
  pC <- ggplot2::ggplot(vp, ggplot2::aes(species, value, fill = component)) +
    ggplot2::geom_col(position = "fill", width = 0.75, na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = c("Species mean" = .nature_palette()$blue,
                                          "Local adaptation" = .nature_palette()$orange,
                                          "Plasticity" = .nature_palette()$green,
                                          "Residual" = .nature_palette()$purple),
                               na.value = "grey75") +
    ggplot2::labs(x = "Species", y = "Proportion") + .nature_theme(7)
  pD <- ggplot2::ggplot(V, ggplot2::aes(LA_PL_ratio, species, colour = species)) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "grey60") +
    ggplot2::geom_point(size = 2.2, na.rm = TRUE) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "LA:PL ratio (log scale)", y = NULL) +
    .nature_theme(7) + ggplot2::theme(legend.position = "none")
  pE <- ggplot2::ggplot(m$reaction_norms,
    ggplot2::aes(environment, response, colour = genotype, group = genotype)) +
    ggplot2::geom_point(size = 1.4, na.rm = TRUE) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                         linewidth = 0.45, na.rm = TRUE) +
    ggplot2::facet_wrap(~species, nrow = 1) +
    ggplot2::labs(x = "Environmental gradient", y = "Trait shift") +
    .nature_theme(6.3)
  pF <- ggplot2::ggplot(m$feedback_strength,
    ggplot2::aes(stats::reorder(path, feedback_strength), feedback_strength, fill = path)) +
    ggplot2::geom_col(width = 0.65, na.rm = TRUE) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Eco-evolutionary path", y = "Feedback strength") +
    .nature_theme(7) +
    ggplot2::theme(legend.position = "none")
  panels <- list(
    .panel("A", "Species-level beta vs. population-level beta", pA, 0.008, 0.555, 0.310, 0.370),
    .panel("B", "Within-species niche divergence", pB, 0.323, 0.555, 0.310, 0.370),
    .panel("C", "Partitioning variance in niche differences", pC, 0.638, 0.555, 0.354, 0.370),
    .panel("D", "Local adaptation to plasticity ratio", pD, 0.008, 0.145, 0.310, 0.390),
    .panel("E", "Genetic niche signal and reaction norms", pE, 0.323, 0.145, 0.310, 0.390),
    .panel("F", "Eco-evolutionary feedback diagnostics", pF, 0.638, 0.145, 0.354, 0.390)
  )
  .nature_page("Figure 6. Population-level and dynamic eco-evolution",
               "Within-species divergence, local adaptation, plasticity, and feedbacks",
               panels)
}

.plot_validation_dashboard_nature <- function(x, ...) {
  m <- if (inherits(x, "hmsc_validation_metrics")) x else calc_validation_metrics(x, ...)
  axes <- .nature_analysis_axes(m$source)
  pA <- ggplot2::ggplot(m$simulation_recovery,
    ggplot2::aes(true_value, estimated_value)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::geom_point(colour = .nature_palette()$blue, size = 0.9, alpha = 0.70,
                        na.rm = TRUE) +
    ggplot2::facet_wrap(~metric, nrow = 1) +
    ggplot2::labs(x = "True", y = "Estimated") + .nature_theme(6.3)
  U <- m$posterior_uncertainty_propagation
  if ("axis" %in% names(U)) U <- U[U$axis %in% c(axes, "all"), , drop = FALSE]
  pB <- ggplot2::ggplot(U,
    ggplot2::aes(value, fill = axis)) +
    ggplot2::geom_histogram(bins = 20, colour = "white", linewidth = 0.15,
                            na.rm = TRUE) +
    ggplot2::facet_wrap(~metric, scales = "free") +
    ggplot2::labs(x = "Posterior uncertainty", y = "Count") + .nature_theme(6.3)
  pC <- ggplot2::ggplot(m$tree_uncertainty_sensitivity,
    ggplot2::aes(metric, value, fill = metric)) +
    ggplot2::geom_violin(width = 0.7, na.rm = TRUE) +
    ggplot2::geom_point(size = 1.4, na.rm = TRUE) +
    ggplot2::labs(x = "Metric", y = "Standardized units") + .nature_theme(7) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
                   legend.position = "none")
  sens <- rbind(.standard_metric_table(m$prior_sensitivity, "prior"),
                .standard_metric_table(m$environment_omission_check, "environment"),
                .standard_metric_table(m$spatial_confounding_check, "spatial"))
  pD <- ggplot2::ggplot(sens, ggplot2::aes(value, metric, colour = source_group)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    ggplot2::geom_point(size = 2, na.rm = TRUE) +
    ggplot2::labs(x = "Percent change / sensitivity", y = NULL) + .nature_theme(7)
  pE <- ggplot2::ggplot(m$block_cross_validation,
    ggplot2::aes(block, value, fill = metric)) +
    ggplot2::geom_col(position = "dodge", width = 0.7, na.rm = TRUE) +
    ggplot2::labs(x = "Spatial fold", y = "Performance") + .nature_theme(7)
  ppc_mcmc <- rbind(.standard_metric_table(m$posterior_predictive_check, "PPC"),
                    .standard_metric_table(m$mcmc_diagnostics, "MCMC"))
  pF <- ggplot2::ggplot(ppc_mcmc,
    ggplot2::aes(stats::reorder(metric, value), value, fill = source_group)) +
    ggplot2::geom_col(position = "dodge", width = 0.7, na.rm = TRUE) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~source_group, scales = "free_y") +
    ggplot2::labs(x = "Diagnostic", y = "Value") + .nature_theme(6.5) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5.5))
  key <- .nature_key_grob(c("Recovery checks calibrate inference.",
                            "Uncertainty is propagated into derived metrics.",
                            "Model sensitivity highlights fragile outputs.",
                            "PPC and MCMC diagnostics guard reliability."))
  panels <- list(
    .panel("A", "Simulation recovery", pA, 0.008, 0.555, 0.330, 0.370),
    .panel("B", "Posterior uncertainty propagation", pB, 0.343, 0.555, 0.310, 0.370),
    .panel("C", "Tree uncertainty sensitivity", pC, 0.658, 0.555, 0.334, 0.370),
    .panel("D", "Model sensitivity checks", pD, 0.008, 0.145, 0.330, 0.390),
    .panel("E", "Block cross-validation", pE, 0.343, 0.145, 0.255, 0.390),
    .panel("F", "Posterior predictive check and MCMC diagnostics", pF, 0.603, 0.145, 0.265, 0.390),
    .panel("", "Key takeaways", key, 0.875, 0.145, 0.117, 0.390)
  )
  .nature_page("Figure 7. Validation and robustness",
               "How reliable are the niche, phylogenetic, evolutionary, and trait-based diagnostics?",
               panels)
}
