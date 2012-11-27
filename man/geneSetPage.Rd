\name{geneSetPage}
\alias{geneSetPage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  A function to create an HTML page for each gene set, as well as the
  HTML pages for each significant gene set.
}
\description{
  This is intended to be an internal function to \code{runRomer}, and is
  not intended to be called by end users. However, the \dots argument to
  \code{runRomer} allows one to pass arguments to lower level functions,
  so the arguments are described here.
}
\usage{
geneSetPage(rslts, genesets, eset, fit, file, cutoff = 0.05, dir = ".",
subdir = ".", columns = NULL, colnames = NULL, col = NULL, caption =
NULL, annot = NULL, nam = NULL, fitind = NULL, bline = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rslts}{
    The results from running \code{\link{romer}} on one gene set.
}
  \item{genesets}{
    Character. A vector of gene symbols for one gene set.
}
  \item{eset}{
   An \code{\link{ExpressionSet}} containing normalized, summarized
  gene expression data.
}
  \item{fit}{
 An \code{MArrayLM} object, containing the fitted data.
}
  \item{file}{
    Filename for the resulting HTML page.
}
  \item{cutoff}{
    Numeric. The cutoff for significance for a given gene set. Defaults
    to 0.05.
}
  \item{dir}{
    The directory to write the results. Defaults to the working directory.
}
  \item{subdir}{
    The subdirectory to write the individual gene set results. Defaults
    to the working directory. 
}
  \item{columns}{
    Numeric. The columns of the \code{ExpressionSet} to use for the
    individual gene set output pages. See \code{dataAndHeatmapPage} for
    more information.
}
  \item{colnames}{
    Character. Alternative column names for the resulting heatmap. See \code{dataAndHeatmapPage} for
    more information.
}
  \item{col}{
    A vector of colors for the heatmap. Defaults to \code{\link{bluered}}.
}
  \item{caption}{
    Caption to put at the top of the HTML page.
}
  \item{annot}{
    Character. The name of the annotation package for the arrray
    used. If \code{NULL}, the row names of the heatmap will be the
    probeset IDs unless the \sQuote{nam} argument below is used.
}
  \item{nam}{
    Character. The row names for the resulting heatmap. If \code{NULL}
    and no \sQuote{annot} argument is used, the probeset IDs will be
    used for the row names.
}
  \item{fitind}{
    Numeric. The columns of the \code{MArrayLM} object to use for
    the individual HTML tables.
  }
  \item{bline}{Defaults to \code{NULL}. Otherwise, a numeric vector
    indicating which columns of the data are the baseline samples. The
    data used for the heatmap will be centered by subtracting the mean
    of these columns from all data.
    }
  \item{\dots}{
    Allows arguments to be passed to lower-level functions. See \code{dataAndHeatmapPage} and
    \code{gsHeatmap} for available arguments.
}
}
\details{
  This function creates a \sQuote{midlevel} HTML table that contains
  each gene set that was significant, with a link to an HTML table that
  shows data for each gene in that gene set (with annotation), as well
  as a heatmap showing the expression levels. Normally this is not run
  by end users, but is called as part of the \code{runRomer}
  function.
}
\value{
  Nothing is returned. Called only for the side effect of creating HTML tables.
}
\author{
  James W. MacDonald <jmacdon@u.washington.edu>
}
