# Originally extracted from the app.model::compare_geneset_shifts function

#' Visualizes geneset "activity" across a number of comparisons.ks
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
#' This will eventually go into the FacileAnalysis package
#'
#' @export
#' @param x A single FacileTtestFseaAnalysisResult or a named list of them.
#' @param genesets character vector of geneset names to extract from the
#'   GeneSetDb objects from each results.
#' @param ... arbitrary number of named ffsea results
#' @param stats the name of the gsea result to pull pvalues from
#' @param colors aaaaarrrghhhhhhhhh
#' @param with_stats A named character vector that specifies the gsea statistics
#'   to print in each of the results. The values correspond to the column names
#'   from the result table of the individual GSEA result, and the `names()` are
#'   the labels you want to use for that statistic. How do you know what columns
#'   are available for printing? Look at the column names in the table returned
#'   by `tidy(ffsea.resut, name = gsea.method)`. By default, the nominal
#'   pvalue and FDR are printed. To disable, set to `NULL`.
#' @examples
#' # We'll setup two ffsea (GSEA) results and plot geneset effects from each
#' efds <- FacileData::exampleFacileDataSet()
#' gdb.h <- sparrow::getMSigGeneSetDb("H", "human", id.type = "entrez")
#' dge <- list(
#'   crc = efds %>%
#'     FacileData::filter_samples(indication == "CRC") %>%
#'     FacileAnalysis::flm_def(
#'       covariate = "sample_type", numer = "tumor", denom = "normal",
#'       batch = "sex") %>%
#'     FacileAnalysis::fdge(method = "voom"),
#'  blca = efds %>%
#'     FacileData::filter_samples(indication == "BLCA") %>%
#'     FacileAnalysis::flm_def(
#'       covariate = "sample_type", numer = "tumor", denom = "normal",
#'       batch = "sex") %>%
#'     FacileAnalysis::fdge(method = "voom"))
#' gsea <- lapply(dge, FacileAnalysis::ffsea, gdb.h, "fgsea")
#'
#' gs.1 <- geneset_shifts(
#'   gsea,
#'   c("beta cells" = "HALLMARK_PANCREAS_BETA_CELLS",
#'     "angiogenesis" = "HALLMARK_ANGIOGENESIS"),
#'   gsea.method = "fgsea",
#'   columns = "genesets",
#'   with_stats = c("pvalue" = "pval", "FDR" = "padj", NES = "NES"))
#' gs.1$plot
#'
#' gs.2 <- geneset_shifts(
#'   gsea,
#'   c("beta cells" = "HALLMARK_PANCREAS_BETA_CELLS",
#'     "angiogenesis" = "HALLMARK_ANGIOGENESIS"),
#'   gsea.method = "fgsea",
#'   columns = "results",
#'   with_stats = c("pvalue" = "pval", "FDR" = "padj", NES = "NES"))
#' gs.2$plot
geneset_shifts <- function(x, genesets, gsea.method = NULL,
                           colors = NULL, ymin = 0.005, ymax = 0.025,
                           bw.bg = 0.7, bw.gs = 0.9,
                           ribbon_wrap_n = 18,
                           legend.position = "none",
                           columns = c("genesets", "results"),
                           xlims = NULL,
                           facet_scales = "free",
                           trim = 0.02,
                           with_stats = c("pvalue" = "pval", "FDR" = "padj.by.collection"),
                           ...) {
  if (is(x, "FacileTtestFseaAnalysisResult")) {
    x <- list(result = x)
  } else {
    x <- assert_list(x, min.len = 1, names = "unique")
  }
  kosher <- sapply(x, is, "FacileTtestFseaAnalysisResult")
  if (any(!kosher)) {
    stop("Illegal objects in x, all elements must be a ",
         "FacileTtestFseaAnalysisResult")
  }
  if (is.null(gsea.method)) {
    gsea.method <- sparrow::resultNames(FacileAnalysis::result(x[[1L]]))[1L]
  }
  assert_string(gsea.method)
  columns <- match.arg(columns)

  assert_number(trim, lower = 0, upper = 0.4, null.ok = TRUE)

  for (cname in names(x)) {
    mg.res <- FacileAnalysis::result(x[[cname]])
    assert_choice(gsea.method, sparrow::resultNames(mg.res))
  }

  dat.all <- lapply(names(x), function(cname) {
    gsea <- x[[cname]]
    mg.res <- FacileAnalysis::result(gsea)
    dge <- FacileAnalysis::param(gsea, "x")

    gs.dat <- lapply(genesets, function(gs) {
      sparrow::geneSet(mg.res, name = gs) %>%
        transmute(group = "geneset", symbol, feature_id, logFC, pval, padj) %>%
        bind_rows(FacileAnalysis::tidy(dge) %>%
                    transmute(group = "transcriptome", symbol,
                              feature_id, logFC, pval, padj)) %>%
        as_tibble() %>%
        mutate(dataset = cname, name = gs,
               y = runif(nrow(.), ymin, ymax))
    }) %>% bind_rows()

    stat.dat <- sparrow::result(mg.res, gsea.method) %>%
      filter(name %in% genesets) %>%
      mutate(dataset = cname, .before = 1L)
      # transmute(dataset = cname, name, logFC = mean.logFC, pval, padj,
      #           pvalue = sprintf("pvalue: %.02f", pval))

    # if (!is.null(trim)) {
    #   lfc.bg <- filter(gs.dat, group == "transcriptome")$logFC
    #   lfc.gs <- filter(gs.dat, group != "transcriptome")$logFC
    #   f.posts <- quantile(lfc.bg, c(trim, 1 - trim))
    #   f.posts[1] <- min(f.posts[1], min(lfc.gs))
    #   f.posts[2] <- max(f.posts[2], max(lfc.gs))
    #   gs.dat <- filter(gs.dat, logFC >= f.posts[1] & logFC <= f.posts[2])
    # }

    list(gs = gs.dat, stats = stat.dat)
  })

  dat.gs <- lapply(dat.all, "[[", "gs") %>%
    bind_rows() %>%
    mutate(
      name = factor(name, genesets),
      dataset = factor(dataset, names(x)),
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
  for (bg.name in paste0(names(x), ".density")) {
    if (is.na(colors[bg.name])) colors[bg.name] <- def.cols["geneset"]
  }
  # make sure we hae colors for the points
  for (cname in names(x)) {
    if (is.na(colors[cname])) colors[cname] <- "cornflowerblue"
  }

  # if (is.null(xlims)) {
  #   xlims <- local({
  #     # xl <- range(filter(dat.gs, group == "geneset")$logFC)
  #     xl <- range(dat.gs$logFC)
  #     c(xl[1] - 0.5, xl[2] + 0.5)
  #   })
  #   if (xlims[1L] > -2) xlims[1L] <- -2
  #   if (xlims[2L] < 2) xlims[2L] <- 2
  # }


  if (!is.null(names(genesets))) {
    label <- setNames(names(genesets), unname(genesets))
    nolabel <- is.na(label) | nchar(label) == 0L
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
    ggplot2::scale_y_continuous(position = "right") +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(
      x = "log2FC",
      y = "")

  if (columns == "genesets") {
    gg <- gg +
      ggplot2::facet_grid(
        dataset ~ lname,
        scales = facet_scales,
        switch = "y",
        labeller = ggplot2::labeller(
          lname = ggplot2::label_wrap_gen(ribbon_wrap_n),
          dataset = ggplot2::label_wrap_gen(ribbon_wrap_n)))
  } else {
    gg <- gg +
      ggplot2::facet_grid(
        lname ~ dataset,
        scales = facet_scales,
        switch = "y",
        labeller = ggplot2::labeller(
          lname = ggplot2::label_wrap_gen(ribbon_wrap_n),
          dataset = ggplot2::label_wrap_gen(ribbon_wrap_n)))
  }

  if (is.numeric(xlims)) {
    gg <- gg + ggplot2::xlim(xlims)
  }

  if (!is.null(with_stats)) {
    # values are the  column names from the stats table generated ffsea
    # names are the labels you want to rename the columns to
    labels <- names(with_stats)
    cnames <- unname(with_stats)
    dat.stat$label <- sapply(1:nrow(dat.stat), function(i) {
      paste(sprintf("%s: %0.3f", labels, dat.stat[i, cnames]), collapse = "\n")
    })

    gg <- gg +
      ggplot2::geom_text(
        mapping = ggplot2::aes(
          x = -Inf, y = Inf, label = label),
        hjust = -0.1, vjust = 1.2,
        data = dat.stat)
  }

  gg <- gg + ggplot2::theme(
    legend.position = legend.position,
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank())

  list(
    dat.gs = dat.gs,
    dat.stat = dat.stat,
    plot = gg)
}
