#' Approximate a facile heatmap function
#'
#' This just wraps the [sparrow::mgheatmap2()] function, but can accept a
#' `facile_frame`
#'
#' TODO: support batch correction over data retrieved from `x`
#' @export
fheatmap <- function(x, gdb = NULL, rename_rows = NULL, ...) {
  if (is(x, "facile_frame")) {
    if (is(gdb, "GeneSetDb")) {
      fids <- sparrow::featureIds(gdb)
    } else {
      fids <- NULL
    }
    x <- FacileData::biocbox(x, "DGEList", features = fids)
    x <- edgeR::calcNormFactors(x)
  }
  if (test_character(rename_rows)) {
    stopifnot(
      "need length(2) character vector" = {
        length(rename_rows) == 2
      },
      "character vector references y$genes columns" = {
        all(rename_rows %in% names(x$genes))
      })
    rename_rows <- x$genes[, rename_rows]
  }
  if (!is.null(rename_rows)) {
    assert_multi_class(rename_rows, c("data.frame", "tbl"))
    stopifnot(ncol(rename_rows) == 2)
  }
  sparrow::mgheatmap2(x, gdb = gdb, rename.rows = rename_rows, ...)
}
