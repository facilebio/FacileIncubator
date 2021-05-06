# Originally extracted from the app.model::compare_geneset_shifts function
#
#' Generate geneset shift plots to compare N contrasts
#'
#' This funciton is used primarily to generate geneshift plots to compare
#' the AppSAA model vs 5XFAD, but why not generalize ...
#'
#' We assume that the ffsea results each have a GeneSetDb that has genesets
#' by the same name in them.
#'
#' The row and column order are determined by the order in which ffsea results
#' and genesets are provided in the respective variables.
#'
#' @export
#' @param genesets character vector of geneset names to extract from the
#'   GeneSetDb objects from each results.
#' @param ... arbitrary number of named ffsea results
#' @param stats the name of the gsea result to pull pvalues from
#' @param colors aaaaarrrghhhhhhhhh
geneset_shifts <- function(genesets, ..., gsea.method = NULL,
                           colors = NULL, ymin = 0.005, ymax = 0.025,
                           bw.bg = 0.7, bw.gs = 0.9,
                           ribbon_wrap_n = 18,
                           legend.position = "none",
                           xlims = NULL, with_pvals = TRUE) {
  args <- assert_list(list(...), min.len = 1)
  if (length(args) == 1L) {
    # user may have passed in an already named-list of gsea results
    args <- args[[1L]]
  }
  keep <- sapply(args, is, "FacileTtestFseaAnalysisResult")
  comps <- assert_list(args[keep], min.len = 2, names = "unique")
  if (is.null(gsea.method)) {
    gsea.method <- sparrow::resultNames(FacileAnalysis::result(comps[[1L]]))[1L]
  }
  assert_string(gsea.method)
  for (cname in names(comps)) {
    mg.res <- FacileAnalysis::result(comps[[cname]])
    assert_choice(gsea.method, sparrow::resultNames(mg.res))
  }

  dat.all <- lapply(names(comps), function(cname) {
    gsea <- comps[[cname]]
    mg.res <- FacileAnalysis::result(gsea)
    dge <- param(gsea, "x")

    gs.dat <- lapply(genesets, function(gs) {
      sparrow::geneSet(mg.res, name = gs) %>%
        transmute(group = "geneset", symbol, feature_id, logFC, pval, padj) %>%
        bind_rows(dge %>%
                    tidy() %>%
                    transmute(group = "transcriptome", symbol,
                              feature_id, logFC, pval, padj)) %>%
        as_tibble() %>%
        mutate(dataset = cname, name = gs,
               y = runif(nrow(.), ymin, ymax))
    }) %>% bind_rows()

    stat.dat <- FacileAnalysis::result(mg.res, gsea.method) %>%
      filter(name %in% genesets) %>%
      transmute(dataset = cname, name, logFC = mean.logFC, pval, padj,
                pvalue = sprintf("pvalue: %.02f", pval))

    list(gs = gs.dat, stats = stat.dat)
  })

  dat.gs <- lapply(dat.all, "[[", "gs") %>%
    bind_rows() %>%
    mutate(
      name = factor(name, genesets),
      dataset = factor(dataset, names(comps)),
      group.density = paste0(dataset, ".density")) %>%
    as_tibble()

  dat.stat <- lapply(dat.all, "[[", "stats") %>%
    bind_rows() %>%
    as_tibble() %>%
    mutate(name = factor(name, levels(dat.gs$name)),
           dataset = factor(dataset, levels(dat.gs$dataset)))

  def.cols <- c(
    transcriptome = "lightgrey",
    geneset = "orange")

  if (is.null(colors)) {
    colors <- def.cols
  }
  if (is.na(colors["transcriptome"])) {
    colors["transcriptome"] <- def.cols["transcriptome"]
  }
  if (is.na(colors["geneset"])) {
    colors["geneset"] <- def.cols["geneset"]
  }
  # make sure we have colors for the foreground density
  for (bg.name in paste0(names(comps), ".density")) {
    if (is.na(colors[bg.name])) colors[bg.name] <- def.cols["geneset"]
  }
  # make sure we hae colors for the points
  for (cname in names(comps)) {
    if (is.na(colors[cname])) colors[cname] <- "cornflowerblue"
  }

  if (is.null(xlims)) {
    xlims <- local({
      xl <- range(filter(dat.gs, group == "geneset")$logFC)
      c(xl[1] - 0.5, xl[2] + 0.5)
    })
    if (xlims[1L] > -2) xlims[1L] <- -2
    if (xlims[2L] < 2) xlims[2L] <- 2
  }


  if (!is.null(names(genesets))) {
    label <- setNames(names(genesets), unname(genesets))
    nolabel <- is.na(label)
    if (any(nolabel)) {
      label[nolabel] <- names(label)[nolabel]
    }
    dat.gs$lname <- factor(label[dat.gs$name], unname(label))
    dat.stat$lname <- factor(label[dat.stat$name], unname(label))
  } else {
    dat.gs$lname <- dat.gs$name
    dat.stat$lname <- dat.stat$name
  }

  gg <- ggplot2::ggplot(
    dat.gs,
    ggplot2::aes(x = logFC)) +
    ggplot2::geom_density(
      color = colors["transcriptome"],
      fill = colors["transcriptome"],
      alpha = 0.8,
      bw = bw.bg,
      data = filter(dat.gs, group != "geneset")) +
    ggplot2::geom_vline(
      xintercept = 0, color = "darkgrey", linetype = "dashed",
      size = 0.5) +
    ggplot2::geom_point(
      ggplot2::aes(y = y, color = dataset, text = symbol),
      data = filter(dat.gs, group == "geneset")) +
    ggplot2::geom_density(
      ggplot2::aes(color = group.density), size = 1,
      bw = bw.gs,
      data = filter(dat.gs, group == "geneset")) +
    ggplot2::facet_grid(
      # dataset ~ name,
      dataset ~ lname,
      scales = "free_y",
      switch = "y",
      labeller = ggplot2::labeller(
        lname = ggplot2::label_wrap_gen(ribbon_wrap_n))) +
    # ggplot2::xlim(-8.5, 15.5) +
    scale_y_continuous(position = "right") +
    ggplot2::xlim(xlims) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(
      x = "log2FC",
      y = ""
    )

  if (with_pvals) {
    gg <- gg +
      ggplot2::geom_text(
        mapping = aes(x = -Inf, y = Inf, label = pvalue),
        hjust = -0.1, vjust = 1.2,
        data = dat.stat)
  }
  gg <- gg + ggplot2::theme(legend.position = legend.position)

  list(
    dat.gs = dat.gs,
    dat.stat = dat.stat,
    plot = gg)
}
