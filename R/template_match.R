
#' Score genes against a template of expression matching.
#'
#' Follows main idea from [MApckg::PatternHunter()]. Thomas pointed out that
#' this is called "template matching" in the field, which was (I think)
#' originally coined by Pavlidis.
#'
#' Other resources:
#' * Analysis of strain and regional variation in gene expression in mouse brain
#'   https://doi.org/10.1186/gb-2001-2-10-research0042 (pavlidis et al.)
#'
#' * Gene Time Eχpression Warper: a tool for alignment, template matching and
#'   visualization of gene expression time series.
#'   https://doi.org/10.1093/bioinformatics/bti787
#'
#' * Post hoc pattern matching: assigning significance to statistically defined
#'   expression patterns in single channel microarray data.
#'   https://doi.org/10.1186/1471-2105-8-240
#'
#' This is a bad name for this function, but just putting it down here.
#'
#' GOTO: FacileAnalysis
#'
#' @section Usage Examples:
#' Imagine we have several timepoints encoded in a 'group' covariate.
#' you can run an nova to find genes that have some association w/ time,
#' then find ones within that result that match your template (let's say
#' ones that increase over time
#'
#' ```r
#' increasing <- fds %>%
#'   flm_def("group") %>%
#'   fdge() %>%
#'   template_match("ascending)
#' ```
#'
#' @export
#' @param x A FacileAnalysis result (likely FacileAnovaAnalysisResult).
#'   We assume the formula was an intercept/effect model.
#' @param template Defines the template to match against. Can be specified by
#'   name for simple template (like `"ascending"` or `"descending"`), or a
#'   numeric vector that provides the pattern, ie. c(1, 2, 3, 2, 1) for a
#'   "peak". Default `"ascending"`.
#' @return a tibble of features that match the template pattern
template_match <- function(x, template = "ascending",
                           cor_method = c("spearman", "pearson", "kendall"),
                           ...) {
  UseMethod("template_match", x)
}

#' @noRd
#' @export
template_match.FacileAnovaAnalysisResult <- function(
    x, template = "ascending",
    cor_method = c("spearman", "pearson", "kendall"),
    ...) {
  assert_class(x, "FacileAnovaAnalysisResult")
  if (is.character(template)) {
    assert_choice(template, c("ascending", "descending"))
  }
  cor_method <- match.arg(cor_method)
  xres <- tidy(x)

  # we need to define the order of the groups before we look for the template
  # and assume here that the order follows the order of appearance of the
  # group-level covariates in the anova result (which follows factor levels if
  # the thing is a factor)
  group.vals <- local({
    samples. <- samples(x)
    group.var <- param(model(x), "covariate")
    group.vals <- samples.[[group.var]]
    if (is.factor(group.vals)) group.vals <- levels(group.vals)
    intersect(colnames(xres), group.vals)
  })

  M <- cbind(0, as.matrix(xres[, group.vals, drop = FALSE]))
  if (!is.numeric(template)) {
    if (template == "ascending") {
      template <- seq(ncol(M))
    } else if (template == "descending") {
      template <- ncol(M):1
    }
  }
  assert_numeric(template, len = ncol(M))

  cstats <- suppressWarnings({
    apply(M, 1, function(vals) {
      cres <- cor.test(vals, template, method = cor_method)
      c(cres$estimate, cres$stat, cres$p.value)
    })
  })
  cstats <- t(cstats)
  colnames(cstats) <- c("estimate", "statistic", "pval")

  # bring the anova stats with you
  ares.out <- select(xres, feature_id, meta, symbol:padj)
  ares.out <- rename(ares.out, pval_anova = pval, padj_anova = padj)

  out <- bind_cols(
    tibble(feature_id = xres[["feature_id"]]),
    as.data.frame(cstats)) %>%
    mutate(padj = p.adjust(pval, "BH")) %>%
    left_join(ares.out, by = "feature_id")
  # arrange features so that most ascending are at top and most descending
  # at bottom
  ordered <- bind_rows(
    filter(out, estimate >= 0) %>% arrange(desc(estimate), pval_anova),
    filter(out, estimate < 0) %>% arrange(desc(estimate), desc(pval_anova)))
  ordered
}
