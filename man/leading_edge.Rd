% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gsea-utils.R
\name{leading_edge}
\alias{leading_edge}
\title{Extracts the leading_edge genes from one or more results.}
\usage{
leading_edge(..., .combine_by = c("union", "intersect"))
}
\arguments{
\item{...}{A single (or list of) \code{FacileTtestAnalysisResult} objects that
have been run with \code{"fgsea"} that we will extract the leading edge genes
from per geneset.}

\item{.combine_by}{When multiple \code{FacileTtestAnalysisResult} objects are
provided in \code{...}, are the leadeing edge genes the \code{"union"} of the genes
from each geneset across results, or the \code{"intersect"}-ion?}
}
\description{
Only works with FacileFseaTtestAnalysisResult that ran \code{methods = "fgsea"}.
If more than one is passed in, the union of genes in the leading edge will
be returned. (or intersection, defined by the \code{combine_by} parameter)
}
\details{
FacileFseaTtestAnalysisResult must be passed in as named arguments if there
is more than one.
}
