#' Save a ggplot object to PNG and PDF
#'
#' @param plot A ggplot object.
#' @param file_base Output path without extension.
#' @param width Width in inches.
#' @param height Height in inches.
#' @param dpi PNG resolution.
#' @return Character vector of written files.
#' @export
hee_save_plot <- function(plot, file_base, width = 7.2, height = 4.8, dpi = 300) {
  .require_pkg("ggplot2", "saving plots")
  dir.create(dirname(file_base), recursive = TRUE, showWarnings = FALSE)
  png_file <- paste0(file_base, ".png")
  pdf_file <- paste0(file_base, ".pdf")
  ggplot2::ggsave(png_file, plot = plot, width = width, height = height,
                  dpi = dpi, bg = "white")
  ggplot2::ggsave(pdf_file, plot = plot, width = width, height = height,
                  bg = "white")
  c(png = normalizePath(png_file, winslash = "/", mustWork = FALSE),
    pdf = normalizePath(pdf_file, winslash = "/", mustWork = FALSE))
}

