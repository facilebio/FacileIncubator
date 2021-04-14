# See also geneset_shifts.R

#' Compare two (or more) ffsea results.
#'
#' This only works with ffsesa results generated from a ttest dge analysis
#'
#' @export
#' @param ... you can provide more FacileFseaTtestAnalysisResult objects
#' @param genesets The genesets to show. If `NULL` (default), takes the `ntop`
#'   genesets by pvalue. Users can provide a collection,name[,group] tibble
#'   that pulls out specific genesets, and optionally groups them by
#'   user-defined categories.
#' @param features The features used to calucate the effect size of the gsea
#'   result. This can be either `"all"`, or `"leading_edge"` to have this method
#'   pull the genes, or a long data.frame with collection,name,feature_id to
#'   pull these out manually using some externally specified approach.
compare.FacileFseaTtestAnalysisResult <- function(x, y, ..., method = NULL,
                                                  score = c("logFC", "t", "NES"),
                                                  features = NULL) {

}
