#' Classify samples based on expression of X vs Y chromosome genes.
#'
#' The massiR package does this for microarrays
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4080740/
#'
#' @export
#' @examples
#' y <- exampledgelist
#' classify_sex(y, species = "human")
classify_sex <- function(x, sex_genes = NULL, ...) {
  UseMethod("classify_sex", x)
}

#' @noRd
#' @export
#' @param x a log normalized matrix of expression
#' @param female logical vector indicating which rows corrsepond to
#'   female-specific genes
classify_sex.matrix <- function(x, sex_genes = NULL, ...) {
  assert_data_frame(sex_genes)
  assert_character(sex_genes[["feature_id"]])
  sex_genes[["sex"]] <- assert_set_equal(
    tolower(sex_genes[["sex"]]),
    c("male", "female")
  )
  
  xsex <- x[rownames(x) %in% sex_genes$feature_id, , drop = FALSE]
  sg <- dplyr::filter(sex_genes, .data$feature_id %in% rownames(xsex))
  
  female <- sg$sex == "female"
  nf <- sum(female)
  if (nf == 0) {
    stop("No female marker genes found in dataset")
  }
  if (nf == nrow(xsex)) {
    stop("No male marker genes found in dataset")
  }
  
  xsex <- xsex[sg$feature_id, , drop = FALSE]
  fscore <- colMeans(xsex[female, , drop = FALSE])
  mscore <- colMeans(xsex[!female, , drop = FALSE])
  
  fdiff <- fscore - mscore
  pred <- ifelse(fdiff > 1, "female", "ambiguous")
  pred <- ifelse(fdiff < -1, "male", pred)
  pred <- dplyr::tibble(
    sample_id = colnames(xsex),
    fdiff_score = fdiff,
    predicted_sex = pred
  )
  
  datadf <- lapply(1:ncol(x), function(i) {
    dplyr::mutate(sg, sample_id = colnames(x)[i])
  }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(value = as.vector(xsex))
  
  list(
    data = datadf,
    prediction = pred
  )
}

#' @noRd
#' @export
classify_sex.DGEList <- function(
    x,
    sex_genes = NULL,
    species = NULL,
    plot = TRUE,
    ...
) {
  assert_class(x, "DGEList")
  if (is.null(sex_genes)) {
    sex_genes <- sex_genes(species)
  }
  x <- edgeR::calcNormFactors(x)
  X <- edgeR::cpm(x, log = TRUE, prior.count = 3)
  out <- classify_sex(X, sex_genes = sex_genes)
  if (plot) {
    xdat <- out$data |>
      dplyr::inner_join(out$prediction, by = "sample_id") |>
      dplyr::mutate(
        label = ifelse(predicted_sex == "ambiguous", sample_id, "")
      )
    gg <- ggplot2::ggplot(xdat) +
      ggplot2::aes(x = predicted_sex, y = value) +
      ggplot2::geom_boxplot(outlier.size = 0) +
      ggplot2::geom_jitter(ggplot2::aes(color = predicted_sex), width = 0.2) +
      ggplot2::facet_wrap(~symbol) +
      ggrepel::geom_label_repel(ggplot2::aes(label = label))
    out$plot <- gg
    print(gg)
  }
  out
}

#' @noRd
#' @export
classify_sex.SummarizedExperiment <- function(
    x,
    sex_genes = NULL,
    species = NULL,
    plot = TRUE,
    ...
) {
  assert_class(x, "SummarizedExperiment")
  anames <- SummarizedExperiment::assayNames(x)
  aname <- if ("counts" %in% anames) "counts" else anames[1L]
  y <- edgeR::DGEList(
    counts = SummarizedExperiment::assay(x, aname),
    samples = as.data.frame(SummarizedExperiment::colData(x)),
    genes = as.data.frame(SummarizedExperiment::rowData(x)),
  )
  classify_sex(y, species = species, plot = plot, ...)
}

#' Genes associated with sex-specific expression.
#'
#' The genes used in massiR seeded this:
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4080740/
#'
#' @export
#' @return a tibble of female/male genes with human/mouse ensembl id's
sex_genes <- function(species = "human") {
  spinfo <- sparrow::species_info(species)
  assert_choice(spinfo$alias, c("mouse", "human"))
  id_col <- paste0(spinfo$alias, "_id")
  
  out <- dplyr::tribble(
    ~sex,     ~symbol,   ~human_id,            ~mouse_id,
    "female", "XIST",    "ENSG00000229807",    "ENSMUSG00000086503",
    "male",   "EIF1AY",  "ENSG00000198692",    NA_character_,
    "male",   "NLGN4Y",  "ENSG00000165246",    NA_character_,
    "male",   "DDX3Y",   "ENSG00000067048",    "ENSMUSG00000069045",
    "male",   "UTY",     "ENSG00000183878",    "ENSMUSG00000068457",
    "male",   "KDM5D",   "ENSG00000012817",    "ENSMUSG00000056673",
    "male",   "RPS4Y1",  "ENSG00000129824",    NA_character_,
  ) # fmt: skip
  
  out |>
    dplyr::select(sex, symbol, feature_id = .data[[id_col]]) |>
    dplyr::filter(!is.na(feature_id))
}
