#' Approximate a facile heatmap function
#'
#' This just wraps the [sparrow::mgheatmap2()] function, but can accept a
#' `facile_frame`
#'
#' TODO: support batch correction over data retrieved from `x`
#' @export
#' @param colors a list of named color vectors. Names should correspond to
#'   columns in row or column annotation dataframes(?)
fheatmap <- function(x, assay_name = NULL, gdb = NULL, rename_rows = NULL, ...,
                     colors = NULL) {
  if (is(x, "facile_frame")) {
    if (is(gdb, "GeneSetDb")) {
      fids <- sparrow::featureIds(gdb)
    } else {
      fids <- NULL
    }
    sample.order <- paste(x$dataset, x$sample_id, sep ="__")
    # TODO: 
    #   1. User should be able to specify assay to use for fheatmap
    #   2. the `class` param (DGEList) should be passed in here, with an
    #      attempt to guess what it is if missing, based on assay_type
    x <- FacileData::biocbox(x, "DGEList", features = fids,
                             assay_name = assay_name)
    x <- edgeR::calcNormFactors(x)
    stopifnot(setequal(sample.order, colnames(x)))
    x <- x[, sample.order]
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
  
  dots <- list(...)
  if (!is.null(dots$top_annotation)) {
    assert_character(dots$top_annotation)
    assert_subset(dots$top_annotation, colnames(x$samples))
    ta <- x$samples[, dots$top_annotation, drop = FALSE]
    tcols <- NULL
    if (is.list(colors)) {
      cnames <- intersect(names(colors), colnames(ta))
      if (length(cnames) > 0) {
        tcols <- colors[cnames]
      }
    }
    taa <- ComplexHeatmap::HeatmapAnnotation(df = ta, col = tcols)
    dots$top_annotation <- taa
  }
  dots$x <- x
  dots$gdb <- gdb
  dots$rename.rows <- rename_rows
  do.call(sparrow::mgheatmap2, dots)
}
