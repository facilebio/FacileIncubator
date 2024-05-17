#' Approximate a facile heatmap function
#'
#' This just wraps the [sparrow::mgheatmap2()] function, but can accept a
#' `facile_frame`
#'
#' TODO: support batch correction over data retrieved from `x`
#' @export
#' @param colors a list of named color vectors. Names should correspond to
#'   columns in row or column annotation dataframes(?)
fheatmap <- function(x, assay_name = NULL, gdb = NULL, rename_rows = NULL, ...,
                     colors = NULL) {
  if (is(x, "facile_frame")) {
    if (is(gdb, "GeneSetDb")) {
      fids <- sparrow::featureIds(gdb)
    } else {
      fids <- NULL
    }
    sample.order <- paste(x$dataset, x$sample_id, sep ="__")
    # TODO: 
    #   1. User should be able to specify assay to use for fheatmap
    #   2. the `class` param (DGEList) should be passed in here, with an
    #      attempt to guess what it is if missing, based on assay_type
    x <- FacileData::biocbox(x, "DGEList", features = fids,
                             assay_name = assay_name)
    x <- edgeR::calcNormFactors(x)
    stopifnot(setequal(sample.order, colnames(x)))
    x <- x[, sample.order]
  }
  if (test_character(rename_rows)) {
    stopifnot(
      "need length(2) character vector" = {
        length(rename_rows) == 2
      },
      "character vector references y$genes columns" = {
        all(rename_rows %in% names(x$genes))
      })
    rename_rows <- x$genes[, rename_rows]
  }
  if (!is.null(rename_rows)) {
    assert_multi_class(rename_rows, c("data.frame", "tbl"))
    stopifnot(ncol(rename_rows) == 2)
  }
  
  dots <- list(...)
  if (!is.null(dots$top_annotation)) {
    assert_character(dots$top_annotation)
    assert_subset(dots$top_annotation, colnames(x$samples))
    ta <- x$samples[, dots$top_annotation, drop = FALSE]
    tcols <- NULL
    if (is.list(colors)) {
      cnames <- intersect(names(colors), colnames(ta))
      if (length(cnames) > 0) {
        tcols <- colors[cnames]
      }
    }
    taa <- ComplexHeatmap::HeatmapAnnotation(df = ta, col = tcols)
    dots$top_annotation <- taa
  }
  
  dots$x <- x
  dots$gdb <- gdb
  dots$rename.rows <- rename_rows
  
  do.call(fheatmap2, dots)
}

fheatmap2 <- function(x, gdb = NULL, col = NULL,
                      aggregate.by = c("none", "ewm", "ewz", "zscore"),
                      split = TRUE, scores = NULL, gs.order = NULL,
                      name = NULL, rm.collection.prefix = TRUE,
                      rm.dups = FALSE, recenter = FALSE, rescale = FALSE,
                      center = FALSE, scale = FALSE,
                      uncenter = FALSE, unscale = FALSE, rename.rows = NULL,
                      zlim = NULL, transpose = FALSE, ...) {
  X <- sparrow:::as_matrix(x, ...)
  stopifnot(
    ncol(X) > 1L,
    !any(is.na(X)))
  if (is.null(scores)) {
    aggregate.by <- match.arg(aggregate.by)
  } else {
    stopifnot(
      is.character(aggregate.by),
      length(aggregate.by) == 1L,
      aggregate.by %in% scores$method)
  }
  
  # if (!is.null(gdb) && aggregate.by != "none") {
  #   if (!is(gdb, "GeneSetDb")) {
  #     gdb <- GeneSetDb(gdb)
  #   }
  # }
  
  drop1.split <- missing(split)
  stopifnot(is.logical(split) && length(split) == 1L)
  if (!is.null(scores)) stopifnot(is.data.frame(scores))
  if (!missing(zlim) && !is.null(zlim)) {
    stopifnot(
      is.numeric(zlim),
      length(zlim) == 2L,
      zlim[1L] < zlim[2])
  }
  
  X <- scale_rows(X, center = center, scale = scale)
  center. <- if (missing(uncenter)) attr(X, "scaled:center") else uncenter
  scale. <- if (missing(unscale)) attr(X, "scaled:scale") else unscale
  
  gdbc.df <- NULL
  if (!is.null(gdb)) {
    if (!is.data.frame(gdb)) {
      gdbc <- suppressWarnings(conform(gdb, X, ...))
      gdbc.df <- as.data.frame(gdbc) # keep only genes that matched in gdb.df
    } else {
      # maintain order user wanted.
      gdbc.df <- dplyr::filter(gdb, .data$feature_id %in% rownames(X))
    }
    
    # Order genesets in requested (if any) order
    if (!is.null(gs.order)) {
      assert_character(gs.order, min.len = 1)
      gs.order <- unique(c(gs.order, gdbc.df[["name"]]))
      gs.order <- intersect(gs.order, gdbc.df[["name"]])
      assert_set_equal(gs.order, gdbc.df[["name"]])
      name. <- factor(gdbc.df[["name"]], gs.order)
      gdbc.df <- gdbc.df[order(name.),,drop = FALSE]
    }
    
    # Set this up so we can order the data.frame in the way requested by user
    gdbc.df$key <- sparrow::encode_gskey(gdbc.df)
  }
  
  if (aggregate.by == "none") {
    if (!is.null(gdbc.df)) {
      ridx <- if (rm.dups) unique(gdbc.df$feature_id) else gdbc.df$feature_id
      # We may have a sparse matrix at this point, turning it to dense for now,
      # but need to fix.
      X <- X[ridx,,drop=FALSE]
      if (is.numeric(recenter)) recenter <- recenter[ridx]
      if (is.numeric(center)) center <- center[ridx]
      split <- if (split) gdbc.df$key else NULL
    }
  } else {
    stop("Haven't vetted geneset scoring in this version")
    if (is.null(scores)) {
      X <- scoreSingleSamples(gdb, X, methods = aggregate.by, as.matrix = TRUE,
                              center = FALSE, scale = FALSE,
                              uncenter = center., unscale = scale., ...)
    } else {
      xs <- data.table::setDT(scores[scores[['method']] == aggregate.by,,drop=FALSE])
      xs[, key := sparrow::encode_gskey(xs)]
      xw <- data.table::dcast(xs, key ~ sample_id, value.var = "score")
      xw <- unique(xw, by = "key")
      X <- as.matrix(xw[, -1, with = FALSE])
      rownames(X) <- xw[[1]]
    }
    # If we want to split, it (only?) makes sense to split by collection
    split <- if (split) sparrow::split_gskey(rownames(X))$collection else NULL
  }
  
  if (!isFALSE(recenter) || !isFALSE(rescale)) {
    X <- scale_rows(X, center = recenter, scale = rescale)
    isna <- which(is.na(X), arr.ind = TRUE)
    if (nrow(isna) > 0L) {
      na.rows <- unique(isna[, "row"])
      if (length(na.rows) == nrow(X)) {
        stop("All rows removed after `scale`")
      }
      warning(length(na.rows), " features NA'd during `scale`, ",
              "these are removed", immediate. = TRUE)
      X <- X[-na.rows,,drop = FALSE]
      split <- split[-na.rows]
    }
  }

  dots <- list(...)
  # side_annos <- intersect(c("right_annotation", "left_annotation"), names(dots))
  side_annos <- intersect(c("right_annotation"), names(dots))
  split <- NULL
  for (aname in side_annos) {
    ranno <- dots[[aname]]
    assert_class(ranno, "data.frame")
    rmissing <- setdiff(rownames(X), rownames(ranno))
    if (length(rmissing)) {
      stop("Missing row level annotations for: ", paste(missing, collapse = ","))
    }
    ranno <- ranno[rownames(ranno) %in% rownames(X),,drop = FALSE]
    X <- X[rownames(ranno),,drop=FALSE]
    rcols <- NULL
    if (is.list(colors)) {
      cnames <- intersect(names(colors), colnames(ranno))
      rcols <- colors[cnames]
    }
    if (!is.null(gdbc.df)) {
      split <- factor(gdbc.df$name, unique(gdbc.df$name))
    }
    ra <- ComplexHeatmap::HeatmapAnnotation(
      df = ranno,
      col = rcols,
      which = "row")
    dots[[aname]] <- ra
  }
  
  
  # What kind of colorscale are we going to use?
  # If this is 0-centered ish, we use a red-white-blue scheme, otherwise
  # we use viridis.
  if (is.null(col)) {
    # Is 0 close to the center of the score distribution?
    qtile.X <- quantile(X, c(0.25, 0.75))
    zero.center <- (qtile.X[1L] < 0 && qtile.X[2L] > 0) || any(recenter)
    if (zero.center) {
      if (missing(zlim)) {
        fpost <- quantile(abs(X), 0.975)
        zlim <- c(-fpost, fpost)
      } else if (is.null(zlim)) {
        zlim <- c(min(X), max(X))
      } else {
        stopifnot(zlim[1L] < 0, zlim[2L] > 0)
      }
      col <- circlize::colorRamp2(
        c(zlim[1L], 0, zlim[2L]),
        # c('#1F294E', '#F7F7F7', '#6E0F11')
        c("navy", "white", "firebrick")
      )
    } else {
      if (missing(zlim)) {
        fpost <- quantile(X, c(0.025, 0.975))
      } else if (is.null(zlim)) {
        fpost <- c(min(X), max(X))
      } else {
        stopifnot(all(zlim >= 0), all(zlim <= 1))
        fpost <- quantile(X, zlim)
      }
      # Higher granularity for viridis colorRamp
      breaks <- quantile(X, seq(0, 1, by = 0.05))
      if (fpost[1L] > breaks[2L] || fpost[2L] < breaks[20L]) {
        stop("Illegal values for zlim")
      }
      breaks[1] <- fpost[1]
      breaks[21] <- fpost[2]
      col <- circlize::colorRamp2(breaks, viridis::viridis(21))
    }
  }
  stopifnot(is.function(col))
  
  # if (drop1.split && !is.null(split) && length(unique(split)) == 1L) {
  #   split <- NULL
  # }
  
  # if (rm.collection.prefix) {
  #   if (aggregate.by != 'none') {
  #     rownames(X) <- split_gskey(rownames(X))$name
  #   } else {
  #     if (!is.null(split)) {
  #       # The order of the splits should be preserved up until this point.
  #       # Since this is our final "look" at the split character vector, let's
  #       # set this as a factor with the levels set in the order of their first
  #       # appearance.
  #       split <- split_gskey(split)$name
  #       split <- factor(split, unique(split))
  #     }
  #   }
  # }
  
  ## Catch Heatmap arguments in `...` and build a list do do.call() them down
  ## into the function call.
  dot.args <- list(...)
  hm.args.default <- as.list(formals(Heatmap))
  
  if (is.null(name)) {
    name <- if (aggregate.by == 'none') 'value' else 'score'
  }
  hm.args <- dot.args[intersect(names(dot.args), names(hm.args.default))]
  hm.args[['matrix']] <- X
  hm.args[['col']] <- col
  hm.args[['row_split']] <- split
  hm.args[['name']] <- name
  if (is.null(hm.args[["cluster_row_slices"]]) && !is.null(gs.order)) {
    hm.args[["cluster_row_slices"]] <- FALSE
  }
  
  row.labels <- rownames(X)
  if (!is.null(rename.rows)) {
    has.meta <- is(x, "DGEList") ||
      is(x, "EList") ||
      is(x, "SummarizedExperiment") ||
      is(x, "eSet")
    is.string <- is.character(rename.rows) && length(rename.rows) == 1L
    if (aggregate.by == "none") {
      if (has.meta && is.string) {
        metadf <- fdata(x, as.df = TRUE)
        metadf <- data.frame(rn = rownames(x), to = metadf[[rename.rows]],
                             stringsAsFactors = FALSE)
        if (!is.null(metadf$to)) {
          row.labels <- rownames(renameRows(X, xref = metadf, ...))
        } else {
          warning("rename.rows column not found in metadata for x")
        }
      } else {
        row.labels <- rownames(renameRows(X, rename.rows, ...))
      }
    } else {
      if (!(is.data.frame(rename.rows) && ncol(rename.rows) == 2)) {
        warning("rename.rows parameter must be a 2 column data.frame when ",
                "aggregate.by != 'none'", immediate. = TRUE)
      } else {
        if (rm.collection.prefix && any(grepl(";", rename.rows[[1]]))) {
          rr <- rename.rows
          rr[[1L]] <- sub("^.*;;?", "", rename.rows[[1L]])
          rename.rows <- rbind(rename.rows, rr)
        }
        row.labels <- rownames(renameRows(X, rename.rows, ...))
      }
    }
  }
  hm.args[["row_labels"]] <- row.labels
  
  H <- do.call(ComplexHeatmap::Heatmap, hm.args)
  H
}
