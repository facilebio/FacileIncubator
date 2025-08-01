#' Species identifier patterns
#' @export
species_ensembl_regex <- function() {
  dplyr::tribble(
    ~species,   ~gene,         ~transcript,
    "human",    "ENSG\\d+",    "ENST\\d+",
    "mouse",    "ENSMUSG\\d+", "ENSMUST\\d+",
    "rat",      "ENSRNOG\\d+", "ENSRNOT\\d+",
    "cyno",     "ENSMFAG\\d+", "ENSMFAT\\d+"
  ) # fmt: skip
}

#' Infer spcies from ensembl identifers patterns
#'
#' @export
#' @examples
#' ids <- list(
#'   human = c("ENSG00000229807", "ENSG00000198692", "ENSG00000165246"),
#'   mouse = c("ENSMUSG00000086503", "ENSMUSG00000069045")
#' )
#'
#' lapply(ids, infer_species_from_ens, collapse = FALSE)
#' lapply(ids, infer_species_from_ens, collapse = TRUE)
#'
#' infer_species_from_ens(unlist(ids), collapse = FALSE)
#' infer_species_from_ens(unlist(ids), collapse = TRUE)
infer_species_from_ens <- function(x, collapse = TRUE, ...) {
  assert_character(x, pattern = "^ENS[A-Z]+\\d+$")
  assert_flag(collapse)
  regexes <- species_ensembl_regex()
  out <- dplyr::tibble(feature_id = x, species = "unknown")
  for (i in 1:nrow(regexes)) {
    info <- regexes[i, ]
    matched <- grepl(info$gene, x) | grepl(info$transcript, x)
    out <- out |>
      dplyr::mutate(
        species = dplyr::case_when(
          species == "unknown" & matched ~ info$species,
          species != "unknown" & matched ~ "ambiguous",
          .default = species
        )
      )
  }
  if (collapse) {
    nambig <- sum(out$species == "amibuous")
    nunk <- sum(out$species == "unknown")
    n <- nrow(out)
    if (nambig > 0L || nunk > 0L) {
      msg <- sprintf(
        "%d / %d ambiguous; %d / %d unknown identifiers",
        nambig,
        n,
        nunk,
        n
      )
      warning(msg)
    }
    props <- table(out$species) / n
    majority <- props > 0.75
    out <- if (!any(majority)) "ambiguous" else names(props)[majority]
  }
  
  out
}