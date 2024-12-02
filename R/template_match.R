
#' Score genes against a template of expression matching.
#'
#' Provides scores for genes where increasing values indicate closer
#' concordance to some pattern of interest such as genes that have monotonically
#' increasing, decreasing, or a peak in the middle, etc. This is called
#' "template matching" in the field, which was (I think) originally coined by
#' Pavlidis.
#'
#' This implementation follows the main idea from the `PatternHunter` function
#' in the MApckg package.
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
#' @references \url{https://rdrr.io/github/flajole/MApckg/man/PatternHunter.html}
#' @export
#'
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
  
  if (is.numeric(template)) {
    # The names of the template need to match columnsin xres
    assert_numeric(template, names = "unique")
    missed <- setdiff(names(template), colnames(xres))
    if (length(missed)) {
      stop(
        length(missed), " template groups were not found in the ANOVA result: ", 
        paste(missed, collapse = ","))
    }
    M <- as.matrix(xres[, names(template)])
  } else {
    assert_choice(template, c("ascending", "descending"))
    if (template == "ascending") {
      template <- seq(ncol(M))
    } else if (template == "descending") {
      template <- seq(ncol(M), 1)
    }
    
    # we need to define the order of the groups before we look for the template
    # and assume here that the order follows the order of appearance of the
    # group-level covariates in the anova result (which follows factor levels if
    # the thing is a factor)
    group.vals <- local({
      samples. <- samples(x)
      group.var <- param(model(x), "covariate")
      group.vals <- samples.[[group.var]]
      if (is.factor(group.vals)) group.vals <- levels(group.vals)
      group.vals <- make.names(group.vals)
      group.vals[1] <- paste0("mean.", group.vals[1])
      group.vals[2:length(group.vals)] <- paste0("logFC.", group.vals[2:length(group.vals)])
      intersect(colnames(xres), group.vals)
    })
  
    M <- cbind(0, as.matrix(xres[, group.vals, drop = FALSE]))
  }

  stopifnot(
    is.numeric(template),
    ncol(M) > 2,
    ncol(M) == length(template))
  
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
