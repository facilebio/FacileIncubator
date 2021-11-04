#' Approximate a facile heatmap function
#'
#' This just wraps the [sparrow::mgheatmap2()] function, but can accept a
#' `facile_frame`
#'
#' @export
fheatmap <- function(x, rename_rows = NULL, ...) {
  if (is(x, "facile_frame")) {
    x <- edgeR::calcNormFactors(FacileData::biocbox(x, "DGEList"))
  }
  if (test_character(rename_rows)) {
    stopifnot(
      "need length(2) character vector" = length(rename_rows) == 2,
      "character vector references y$genes columns" = all(rename_rows %in% names(x$genes)))
    rename_rows <- x$genes[, rename_rows]
  }
  if (!is.null(rename_rows)) {
    assert_multi_class(rename_rows, c("data.frame", "tbl"))
    stopifnot(ncol(rename_rows) == 2)
  }
  sparrow::mgheatmap2(x, rename.rows = rename_rows, ...)
}
