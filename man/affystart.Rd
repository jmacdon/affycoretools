\name{affystart}
\alias{affystart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Pre-processing for Affymetrix Data }
\description{
  This function is designed to automatically read in all cel files in a
  directory, make all pre-processing QC plots and compute expression measures.
}
\usage{
affystart(..., filenames=NULL, groups = NULL, groupnames = NULL, plot =
TRUE, pca = TRUE,  squarepca = FALSE, plottype="pdf", express = c("rma",
"mas5", "gcrma"), addname = NULL, output = c("txt", "xls"), annotate =
FALSE, ann.vec =  c("SYMBOL","GENENAME","ENTREZID","UNIGENE","REFSEQ"),
phenoData = new("phenoData"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{ Requires that all variables be named }
  \item{filenames}{ If not all cel files in a directory will be used,
    pass a vector of filenames.}
  \item{groups}{ An integer vector indicating the group assignments
    for the PCA plot. }
  \item{groupnames}{ A character vector with group names for PCA
    legend.}
  \item{plottype}{ What type of plot to save. Can be "pdf","postscript",
    "png","jpeg", or "bmp". Defaults to "pdf". Note that "png" and
    "jpeg" may not be available on a given computer. See the
    help page for \code{capabilities} and  \code{png} for more information.}
  \item{plot}{ Should density and degradation plots be made? Defaults to
    \code{TRUE}.}
  \item{pca}{ Should a PCA plot be made? Defaults to \code{TRUE}. }
  \item{squarepca}{Should the y-axis of the PCA plot be made comparable
    to the x-axis? This may aid in interpretation of the PCA
    plot. Defaults to \code{FALSE}}
  \item{express}{ One of either rma, mas5, gcrma. Defaults to
    rma. Partial matching OK. }
  \item{addname}{ Used to append something to the name of the pca plot
    and the expression values output file (e.g., if function is run
    twice using different methods to compute expression values)}
  \item{output}{Output expression values in either a text ("txt") file
    or an Excel spreadsheet ("xls"). The latter requires the
    \code{xlsReadWrite} package.}
  \item{annotate}{Boolean. Add annotation data to the output file?}
  \item{ann.vec}{A character vector of annotation data to add to the
    output file.}
  \item{phenoData}{A \code{phenoData} object.}
}

\value{
  Returns an \code{exprSet}.
}
\seealso{\code{plotHist}, \code{plotDeg}, \code{plotPCA}
}

\author{James W. MacDonald <jmacdon@med.umich.edu> }

\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{manip }% __ONLY ONE__ keyword per line
