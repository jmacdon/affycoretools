\name{plotDeg}
\alias{plotDeg}
\alias{plotHist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions to Plot Density and RNA Degradation Plots }
\description{
 These functions make density and RNA degradation plots with automatic
 placement of legends.
}
\usage{
plotDeg(dat, filenames = NULL)
plotHist(dat, filenames = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{An \code{AffyBatch} object, or in the case of
    \code{plotHist}, a matrix (e.g., from a call to
    \code{read.probematrix}. Note that \code{plotDeg} requires an
    \code{AffyBatch} object to work correctly. }
  \item{filenames}{Filenames that will be used in the legend of the
    resulting plot. If \code{NULL} (the default), these names will be
    extracted from the sampleNames slot of the \code{AffyBatch} object.}
}
\value{
These functions are called only for the side effect of making the
plots. Nothing else is returned.
}
\author{James W. MacDonald <jmacdon@med.umich.edu>}

\examples{
library("affydata")
data(Dilution)
plotDeg(Dilution)
plotHist(Dilution)
}
\keyword{ hplot }% at least one, from doc/KEYWORDS
