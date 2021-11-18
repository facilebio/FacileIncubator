#' Simplifies ffsea GO analyses
#'
#' We use rrvgo, goSemSim and friends
#'
#' @export
#' @param x A FacileFseaAnalysisResult
#' @param ontology which ontology are were summarizing? BP, MF, or CC
#' @param direction focus on all genes, or only up or down regulated genes
#'   (works for ttest, not ANOVA)
#' @param similarity The method to use to measure similarity,
#'   [see here](https://www.bioconductor.org/packages/release/bioc/vignettes/GOSemSim/inst/doc/GOSemSim.html#semantic-similarity-measurement-based-on-go).
#' @param min_logFC,max_padj parameters used to define signficant genes. If
#'   the TREAT framework was used, defaults to the threshold tested, otherwise
#'   its zero.
go_simplify <- function(x, method = NULL, ontology = c("BP", "MF", "CC"),
                        direction = c("all", "up", "down"),
                        similarity = c("Resnik", "Lin", "Rel", "Jiang", "Wang"),
                        threshold = 0.7, max_padj = 0.10,
                        min_logFC = NULL,
                        score_by = NULL, semdata = NULL, ...,
                        verbose = FALSE) {
  UseMethod("go_simplify", x)
}

#' @noRd
#' @export
go_simplify.default <- function(x, method = NULL,
                                ontology = c("BP", "MF", "CC"),
                                direction = c("all", "up", "down"),
                                similarity = c("Resnik", "Lin", "Rel", "Jiang", "Wang"),
                                threshold = 0.7, max_padj = 0.10,
                                min_logFC = NULL,
                                score_by = NULL,
                                semdata = NULL, ..., verbose = FALSE) {
  stop("default method undefined")
}

#' @noRd
#' @export
go_simplify.FacileFseaAnalysisResult <- function(
  x, method = FacileAnalysis::param(x, "methods")[1L],
  ontology = c("BP", "MF", "CC"),
  direction = c("all", "up", "down"),
  similarity = c("Resnik", "Lin", "Rel", "Jiang", "Wang"),
  threshold = 0.7, max_padj = 0.10, min_size = 1, max_size = Inf,
  score_by = NULL, semdata = NULL, organism = NULL, ..., verbose = FALSE) {
  if (is.null(organism)) organism <- FacileData::organism(x)

  from <- FacileAnalysis::param(x, "x") # the thing ffsea was calculated from

  ontology <- match.arg(ontology)
  direction <- match.arg(direction)

  if (is(from, "FacileAnovaAnalysisResult") && direction != "all") {
    stop("There is no up/down directionality from an ANOVA result")
  }

  x.input <- FacileAnalysis::param(x, "x")
  similarity <- match.arg(similarity)
  assert_number(threshold, lower = 0, upper = 1)
  assert_number(max_padj, lower = 0, upper = 1)

  org.info <- sparrow::species_info(organism)
  orgdb <- sprintf("org.%s.eg.db", org.info$bioc_abbrev)

  mg.method <- .assert_method_ran(x, method, direction = direction, ...)

  # Figure out what score_by metric to use given the gsea method we are using
  # to simplify the GO results from.
  if (method %in% c("ora")) {
    score.opts <- c("pval", "t", "logFC")
  } else if (method %in% c("cameraPR", "camera")) {
    score.opts <- c("logFC", "t", "pval")
  } else if (method %in% c("fgsea")) {
    score.opts <- c("NES", "ES")
  } else {
    stop("Need to inject some knowledge into the GSEA <-> score_by maps")
  }
  if (is.null(score_by)) {
    score_by <- score.opts[1L]
  }
  score_by <- match.arg(score_by, score.opts)

  mstats <- FacileAnalysis::tidy(x, mg.method)
  rstats.all <- transmute(
    mstats,
    collection, name, pval,
    logFC = mean.logFC.trim)

  if (method == "fgsea") {
    rstats.all[["ES"]] <- mstats[["ES"]]
    rstats.all[["NES"]] <- mstats[["NES"]]
  }

  if ("mean.t.trim" %in% colnames(mstats) && !any(is.na(mstats$mean.t.trim))) {
    rstats.all$t <- mstats$mean.t.trim
  }

  # Fishes out the specific branch of the gene ontology tree that we are
  # summarizing here and readjusts the pvalues
  gsets <- .ontology_genesets(x, ontology, ...)
  rstats <- rstats.all %>%
    inner_join(gsets, by = c("collection", "name")) %>%
    mutate(padj = p.adjust(pval, "BH")) %>%
    filter(.data$padj <= .env$max_padj)

  # Let's handle directionality.
  # If we are testing the `ora` method, directionality was handled becasue
  # we are retrieving either `"ora.up"`, `"ora.down"`, or straight-up `"ora"`
  # TODO: This is another place were knowledge of methods to gsea type
  # is useful.
  #
  # For ranking type of GSEA, we'll set the score of the genesets going in the
  # other direction to nothing (0 for score like things, one for pvalues)
  if (method %in% c("camera", "cameraPR", "fgsea")) {
    if (method == "fgsea") {
      up <- rstats[[score_by]] > 0
      lowval <- 0
    } else if (method %in% c("camera", "cameraPR")) {
      lfc.stat <- if (is.numeric(rstats[["t"]])) "t" else "logFC"
      up <- rstats[[lfc.stat]] > 0
      lowval <- if (score_by == "pval") 1 else 0
    }
    nuke <- if (direction == "up") !up else up
    rstats[[score_by]] <- ifelse(nuke, lowval, rstats[[score_by]])
    no.enrich.dir <- all(nuke)
  } else {
    no.enrich.dir <- nrow(rstats) == 0L
  }

  if (no.enrich.dir) {
    # Either there are no genesets going this direction, or significant ones
    # were returned from ora
    return(NULL)
  }

  if (score_by == "pval") {
    rstats[["score"]] <- -log10(rstats[["pval"]])
  } else {
    rstats[["score"]] <- abs(rstats[[score_by]])
  }

  rstats <- arrange(rstats, desc(score))
  scores <- setNames(rstats[["score"]], rstats[["gs_id"]])

  if (verbose) message("Calculating similarity matrix")
  semdata <- .load_semdata(semdata, orgdb, ontology, org.info,
                           verbose = verbose)

  sim.matrix <- rrvgo::calculateSimMatrix(
    names(scores), orgdb = orgdb, semdata = semdata,
    ont = ontology, method = similarity)

  if (verbose) message("Reducing sem similarity")

  sim.reduced <- rrvgo::reduceSimMatrix(sim.matrix, scores,
                                        threshold = threshold, orgdb = orgdb)

  out <- list(
    params = list(
      x = x, method = method, ontology = ontology, direction = direction,
      similarity = similarity, threshold = threshold, max_padj = max_padj,
      score_by = score_by),
    stats = rstats,
    sim_matrix = sim.matrix,
    reduced_terms = sim.reduced)
}

#' Helper function to load semdata object.
#'
#' @noRd
.load_semdata <- function(semdata, orgdb, ontology, org.info = NULL,
                          cache.dir = getOption("semdata.dir", NULL),
                          verbose = FALSE) {
  ontology <- match.arg(ontology, c("BP", "MF", "CC"))
  if (test_string(semdata) &&
      test_file_exists(semdata, "r", extension = "rds")) {
    if (verbose) message("... loading from local fs: ", semdata)
    semdata <- readRDS(semdata)
  }
  if (is(semdata, "GOSemSimDATA")) {
    is.kosher <- isTRUE(ontology == semdata@ont)
    if (!is.null(org.info)) {
      # Check that we got the semdata object, if one was provided
      semdata.species <- subset(semdata@metadata, name == "ORGANISM")$value
      is.kosher <- is.kosher && isTRUE(org.info$species == semdata.species)
    }
    if (is.kosher) {
      if (verbose) message("... semdata is legit")
      return(semdata)
    }
  }

  if (verbose) message("Calling GOSemSim::godata")
  GOSemSim::godata(orgdb, ont = ontology)
}

#' Check that the appropriate ontology category was tested, and returns the
#' GeneSetDb that applies to that.
#'
#' @noRd
.ontology_genesets <- function(x, ontology = c("BP", "MF", "CC"), ...) {
  assert_class(x, "FacileFseaAnalysisResult")
  ontology <- match.arg(ontology)
  mg.res <- FacileAnalysis::result(x)
  gdb.all <- sparrow::geneSetDb(mg.res)

  # We try to match the go ontolgy genesets based on two criteria.
  # 1. Legacy uses of this functionality internally at Denali, when we were
  #    using our multiGSEA package would have had the GO_(BP|MF|CC) set in the
  #    collection of the geneset because this is how it was returned in
  #    multiGSEA::getMSigGeneSetDb with promote_subcategory_to_collection = TRUE
  # 2. Newer functionality uses the sparrow/msigdbr packages to get genesets.
  #    This uses version >= 7.4.1 of MSigDB, which appends the ontology
  #    tree type to the name of the geneset itself (GOBP_<NAME>)
  gsets <- sparrow::geneSets(gdb.all)
  omatched <-
    grepl(paste0("GO_", ontology), gsets[["collection"]]) |
    grepl(paste0("^GO", ontology, "_"), gsets[["name"]])
  gsets <- gsets[omatched,]
  if (nrow(gsets) == 0) {
    stop("The geneset collection was not tested: ", ontology)
  }
  dplyr::select(gsets, collection, name, gs_id)
}


#' @noRd
.assert_method_ran <- function(x, method, direction = c("all", "up", "down"), ...) {
  assert_class(x, "FacileFseaAnalysisResult")
  direction <- match.arg(direction)
  if (method == "ora" && direction != "all") {
    method <- paste0(method, ".", direction)
  }
  mg.res <- FacileAnalysis::result(x)
  match.arg(method, sparrow::resultNames(mg.res))
}
