#' Render a Word report from an R Markdown file
#'
#' @param rmd Path to an R Markdown report.
#' @param output_file Output DOCX path.
#' @param params Optional R Markdown parameters.
#' @return Output DOCX path.
#' @export
hee_render_word_report <- function(rmd, output_file, params = list()) {
  .require_pkg("rmarkdown", "Word report rendering")
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  rmarkdown::render(
    input = rmd,
    output_format = rmarkdown::word_document(),
    output_file = basename(output_file),
    output_dir = dirname(output_file),
    params = params,
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )
  normalizePath(output_file, winslash = "/", mustWork = FALSE)
}
