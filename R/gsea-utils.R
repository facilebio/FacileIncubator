# NOTE: The functionality developed ere was pushed over to the FacileIncubator
#       for further generalization. Don't forget to add any improvements made
#       here to there.
#
#       The code will remain here for analyses within this project.

#' Extracts the leading_edge genes from one or more results.
#'
#' Only works with FacileFseaTtestAnalysisResult that ran `methods = "fgsea"`.
#' If more than one is passed in, the union of genes in the leading edge will
#' be returned. (or intersection, defined by the `combine_by` parameter)
#'
#' FacileFseaTtestAnalysisResult must be passed in as named arguments if there
#' is more than one.
#'
#' @export
leading_edge <- function(..., .combine_by = c("union", "intersect")) {
  args <- list(...)
  .combine_by <- match.arg(.combine_by)

  # Maybe a list of results was passed in as opposed to individual results
  if (length(args) == 1L &&
      is.list(args[[1L]]) &&
      !is(args[[1]], "FacileTtestFseaAnalysisResult")) {
    args <- args[[1L]]
  }
  keep <- sapply(args, is, "FacileTtestFseaAnalysisResult")
  res.all <- args[keep]
  if (length(res.all) == 0L) {
    stop("No FacileTtestFseaAnalysisResult provided")
  } else if (length(res.all) == 1L) {
    if (is.null(names(res.all))) names(res.all) <- "fgsea"
  } else {
    assert_list(res.all, names = "unique")
  }

  le <- lapply(names(res.all), function(gname) {
    gres <- res.all[[gname]]
    gstats <- tidy(gres, "fgsea")
    out <- gstats %>%
      select(collection, name, N = n, leadingEdge) %>%
      tidyr::unnest(leadingEdge) %>%
      mutate(comp = gname) %>%
      select(comp, everything())
    out
  })

  universe <- bind_rows(le) %>%
    rename(feature_id = leadingEdge) %>%
    group_by(collection, name) %>%
    add_count(feature_id, name = "ncomps") %>%
    ungroup() %>%
    select(-comp)

  if (.combine_by == "intersect") {
    out <- filter(universe, ncomps == length(res.all))
    missed <- universe %>%
      distinct(collection, name) %>%
      anti_join(out, by = c("collection", "name"))
    if (nrow(missed)) {
      warning("lost genesets due to 'intersect' criteria")
    }
  } else {
    out <- universe
  }

  distinct(out, collection, name, feature_id, .keep_all = TRUE)
}

#' Aggregates gene-level stats from genes in leading edge.
#'
#' The logFC and t-statistics of the genes from the leading edge
#' (defined in [extract_leading_edge()]) are used to give an effect size of the
#' geneset shift.
#'
#' This is used in place of the ES/NES mojo.
#'
#' @export
#' @param x The FacileFseaAnalysisResult objects where `fgsea` was
#'   ran.
#' @param leading_edge Precaculated leading edge genes, as in the output from
#'   [leading_edge()]. If `NULL`, the leading_edge will be caluclated for all
#'   genesets tested `x` then scored. The scoring can take a long time if there
#'   are many, so you may want to provide a filtered leading_edge set.
#' @param aggregate the column names from the dge results to summarize, defaults
#'   to `c("t", "logFC")`
#' @param stats A name of a gsea method run in `"x"` to append pvalues from.
#'   If `NULL` (default), none will be added.
#' @return a tibble of geneset scores, aggregated by leading edge genesets
#'   for columns in `aggregate`. If `"fgsea"` is a method run in `"x"`, `ES`,
#'   and `NES` columns will be returned as well.
leading_edge_scores <- function(x, leading_edge = NULL,
                                aggregate = c("t", "logFC"),
                                stats = NULL, ...) {
  assert_class(x, "FacileTtestFseaAnalysisResult")

  if (is.null(leading_edge)) leading_edge <- leading_edge(x)
  assert_multi_class(leading_edge, c("data.frame", "tibble"))
  assert_subset(c("collection", "name", "feature_id"), colnames(leading_edge))

  mgres <- result(x)
  lfc <- sparrow::logFC(mgres)

  assert_character(aggregate, min.len = 1L)
  stopifnot(sapply(aggregate, function(x) is.numeric(lfc[[x]])))

  out <- leading_edge %>%
    group_by(collection, name) %>%
    do({
      xstats <- inner_join(., lfc, by = "feature_id")
      res <- sapply(aggregate, function(a) mean(xstats[[a]]), simplify = FALSE)
      bind_cols(res)
    }) %>%
    ungroup()

  rnames <- sparrow::resultNames(mgres)
  if ("fgsea" %in% rnames) {
    nes <- select(tidy(x, "fgsea"), collection, name, ES, NES)
    out <- left_join(out, nes, by = c("collection", "name"))
  }

  if (test_choice(stats, rnames)) {
    gs <- select(tidy(x, stats), collection, name, pval, starts_with("padj"))
    out <- left_join(out, gs, by = c("collection", "name"))
  }

  out
}

