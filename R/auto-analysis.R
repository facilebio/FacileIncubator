#' Executes the predefined analysis defined in an experiment setup table
#'
#' @export
#' @param definitions path to the comparison definition file
#' @param board the pinboard to serialize the results to
#' @examples
#' efds <- FacileData::exampleFacileDataSet()
#' esamples <- FacileData::with_sample_covariates(FacileData::samples(efds))
#' adefs <- read.csv(
#'   system.file("extdata", "example-auto-analysis.csv",
#'               package = "FacileIncubator"))
#' results <- run_analyses(adefs, samples = esamples)
run_analyses <- function(definitions, board = pins::board_temp(), key = NULL,
                         gdb = NULL, ffsea_methods = NULL,
                         semdata = NULL,
                         BPPARAM = BiocParallel::bpparam(), ...) {
  assert_class(BPPARAM, "BiocParallelParam")
  if (test_string(definitions)) {
    assert_file_exists(definitions, "r", extension = "csv")
    definitions <- read.csv(definitions)
  }
  assert_multi_class(definitions, c("data.frame", "tbl"))
  if (!is.null(board)) {
    assert_class(board, "pins_board")
    assert_string(key, null.ok = TRUE)
  }

  # every comparison needs a unique name
  assert_character(definitions[["name"]], min.chars = 5, unique = TRUE)

  # ensure interaction comparisons name previously named direct ones
  direct <- NULL
  direct.comps <- subset(definitions, type == "direct")
  interaction <- NULL
  i.comps <- subset(definitions, type == "interaction")
  if (nrow(i.comps)) {
    assert_subset(i.comps$numer, direct.comps$name)
    assert_subset(i.comps$denom, direct.comps$name)
  }

  # Run direct comparisons
  if (nrow(direct.comps)) {
    direct <- lapply(seq_len(nrow(direct.comps)), function(i) {
      params <- as.list(direct.comps[i,])
      cres <- run_comparison(params, board = board, gdb = gdb,
                             ffsea_methods = ffsea_methods, semdata = semdata,
                             BPPARAM = BPPARAM, ...)
      if (!is.null(board)) {
        pin_write_auto_analysis(board, cres, key, params$name)
      }
      cres
    })
    names(direct) <- direct.comps[["name"]]
  }

  # Run interactions
  if (nrow(i.comps)) {
    interaction <- lapply(seq_len(nrow(i.comps)), function(i) {
      params <- as.list(i.comps[i,])
      cres <- run_comparison(params, board = board, gdb = gdb,
                             ffsea_methods = ffsea_methods, semdata = semdata,
                             direct_dge = direct, BPPARAM = BPPARAM, ...)
      if (!is.null(board)) {
        pin_write_auto_analysis(board, cres, key, params$name)
      }
      cres
    })
    names(interaction) <- i.comps[["name"]]
  }

  list(
    direct = direct,
    interaction = interaction,
    board = board)
}

#' Performs dge, gsea, and GO simplification in one shot for a given analysis.
#'
#' @param params a list of params used for the analysis. This is typically a
#'   single row of a comparison table
#' @param samples the facile_frame that includes the sample info
#' @param board the pinboard to serialize the results to. If you don't wan to
#'   use a fancy remote one, you can always use a local directory.
run_comparison <- function(params, samples = NULL, board = pins::board_temp(),
                           name = NULL,
                           gdb = NULL, ffsea_methods = NULL,
                           go_ontologies = c("BP", "MF", "CC"),
                           go_simplify_method = "ora",
                           go_simplify_sem_sim_method = "Resnik",
                           go_simplify_max_padj = 0.10,
                           go_simplify_sem_threshold =0.6,
                           go_semdata = NULL, BPPARAM = BiocParallel::bpparam(),
                           direct_dge = NULL,
                           ...) {
  if (FALSE) {
    params <- analysis_definitions()[1,]
  }
  params <- lapply(as.list(params), function(x) if (is.na(x)) NULL else x)
  assert_list(params, names = "unique")
  assert_class(BPPARAM, "BiocParallelParam")

  browser()

  if (params$type == "direct") {
    if (is.null(samples)) {
      samples <- study_samples(params$sample_key)
    }
    assert_class(samples, "facile_frame")
    assert_string(params$covariate)
    assert_string(params$numer, null.ok = TRUE)
    assert_string(params$denom, null.ok = TRUE)
    assert_string(params$batch, null.ok = TRUE)
    if (!is.null(params$batch)) {
      params$batch <- trimws(strsplit(params$batch, ",")[[1L]])
    }
    if (test_string(params$subset_covariate)) {
      vals <- trimws(strsplit(params$subset_values, ",")[[1L]])
      keep <- samples[[params$subset_covariate]] %in% vals
      samples <- samples[keep,,drop = FALSE]
    }
    flm <- FacileAnalysis::flm_def(samples, covariate = params$covariate,
                                   numer = params$numer, denom = params$denom,
                                   batch = params$batch)
    dge <- FacileAnalysis::fdge(flm, ...)
  } else {
    cmp.x <- direct_dge[[params$numer]][["dge"]]
    cmp.y <- direct_dge[[params$denom]][["dge"]]
    dge <- FacileAnalysis::compare(cmp.x, cmp.y)
  }

  if (test_character(ffsea_methods)) {
    gsea <- FacileAnalysis::ffsea(dge, gdb, methods = ffsea_methods,
                                  BPPARAM = BPPARAM, ...)
    samples <- FacileData::samples(gsea)
  } else {
    gsea <- NULL
  }

  simplify.go <- !is.null(gsea) &&
    (go_simplify_method %in% sparrow::resultNames(gsea$result)) &&
    test_character(go_ontologies, min.len = 1L) &&
    test_subset(go_ontologies, c("BP", "MF", "CC"))

  if (simplify.go) {
    fds. <- FacileData::fds(samples)
    go.all <- BiocParallel::bplapply(go_ontologies, function(ont) {
      # todo: insert check to see if there are any dge genes up or down, and
      # filter `dirs` accordingly
      dirs <- c("up", "down")
      sapply(dirs, function(direction) {
        if (go_simplify_method == "ora") {
          resname <- paste0("ora.", direction)
        } else {
          resname <- go_simplify_method
        }
        if (resname %in% sparrow::resultNames(FacileAnalysis::result(gsea))) {
          dout <- tryCatch({
            FacileIncubator::go_simplify(
              gsea, method = go_simplify_method, direction = direction,
              ontology = ont,
              similarity = go_simplify_sem_sim_method,
              threshold = go_simplify_sem_threshold,
              max_padj = go_simplify_max_padj,
              semdata = go_semdata[[ont]],
              organism = FacileData::organism(fds.),
              verbose = TRUE)
          }, error = function(e) geterrmessage())
        } else {
          dout <- NULL
        }
        dout
      }, simplify = FALSE)
    }, BPPARAM = BPPARAM)
    names(go.all) <- go_ontologies
  } else {
    go.all <- NULL
  }

  if (!test_string(name)) {
    name <- params$name
  }

  out <- list(
    name = name,
    dge = dge,
    gsea = gsea,
    go = go.all)
  class(out) <- c("FacileAutoAnalysis", class(out))
  out
}

#' The element names found in an auto-analysis result.
#'
#' These are used to iterate over the things to save/load
#' @noRd
.auto_analysis_elements <- function() {
  c("dge", "gsea", "go")
}

#' Analysis types to serialize with FacileAnalysis::fsave/fload
#' @noRd
.auto_analysis_facile <- function() {
  c("dge", "gsea")
}

#' Save the analysis result to a pinboard
#'
#' @export
#' @param board the pinboard
#' @param x the FacileAutoanalysis result
#' @param name the name of the objec to pin. Defaults to `x$name`.
#' @param key a unique string that prefixes all of the objects serialized by
#'   this project
pin_write_auto_analysis <- function(board, x, name = NULL, key = NULL,
                                    elements = NULL, ...) {
  assert_class(board, "pins_board")
  assert_class(x, "FacileAutoAnalysis")
  assert_string(key, null.ok = TRUE)
  if (is.null(name)) name <- x[["name"]]
  assert_string(name)

  if (is.null(elements)) elements <- .auto_analysis_elements()
  assert_subset(elements, .auto_analysis_elements())

  if ("gsea" %in% elements) {
    # the dge result is embedded within the gsea result, so we can skip that
    elements <- setdiff(elements, "dge")
  }

  use.fsave <- .auto_analysis_facile()

  saved <- sapply(elements, function(fname) {
    res <- FacileAnalysis::unfds(x[[fname]])
    pin.name <- paste(name, fname, sep = "___")
    if (!is.null(key)) {
      pin.name <- paste(key, pin.name, sep = "___")
    }
    out <- try({
      pins::pin_write(
        board,
        res,
        name = pin.name,
        type = "qs",
        title = sprintf("'%s' result from '%s' analysis in '%s' project",
                        fname, name, key))
    }, silent = TRUE)
    if (is(out, "try-error")) NULL else out
  })

  saved
}

#' Restore the pinned analysis
#' @export
#' @return a list of analysis results, with `$name` `$dge`, `$gsea`, `$go`
#'   elements.
pin_read_auto_analysis <- function(board, name, key = NULL, definitions = NULL,
                                   elements = NULL, fds. = NULL, ...) {
  assert_class(board, "pins_board")
  assert_string(name)
  assert_string(key, null.ok = TRUE)

  all.pins <- pins::pin_list(board)
  prefix <- paste0(name, "___")
  if (!is.null(key)) {
    prefix <- paste(key, prefix, sep = "___")
  }
  these.pins <- all.pins[grepl(prefix, all.pins)]

  info.all <- dplyr::tibble(
    pin_name = these.pins,
    analysis_name = name,
    result_type = sub(prefix, "", pin_name),
    use_fload = result_type %in% .auto_analysis_facile())

  if (nrow(info.all) == 0L) {
    stop(sprintf("No results found for key:name `%s`:`%s`", key, name))
  }

  if (!is.null(elements)) {
    info.all <- subset(info.all, result_type %in% elements)
  }
  if (nrow(info.all) == 0L) {
    stop("All result types are filtered out")
  }

  if (is.null(fds.) && is.data.frame(definitions)) {
    params <- definitions[definitions$name == name,]
    if (nrow(params) == 1L) {
      fsamples <- study_samples(params$sample_key)
      fds. <- FacileData::fds(fsamples)
    } else {
      warning("Definition not found, no fds is set: ", name)
      fsamples <- NULL
      fds. <- NULL
    }
  }

  out <- sapply(info.all$result_type, function(atype) {
    info <- subset(info.all, result_type == atype)
    out <- pins::pin_read(board, info$pin_name)
    if (info$use_fload) {
      out <- FacileAnalysis::refds(out, fds.)
    }
    out
  }, simplify = FALSE)

  if ("gsea" %in% names(out)) {
    out[["dge"]] <- param(out[["gsea"]], "x")
  }

  out
}
