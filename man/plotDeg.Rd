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
plotDeg(dat, filenames)
plotHist(dat, filenames)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{An \code{AffyBatch} object, or in the case of
    \code{plotHist}, a matrix (e.g., from a call to
    \code{read.probematrix}. Note that \code{plotDeg} requires an
    \code{AffyBatch} object to work correctly. }
  \item{filenames}{A vector of .cel filenames. Usually built using \code{list.celfiles()} }
}
\value{
These functions are called only for the side effect of making the
plots. Nothing else is returned.
}
\author{James W. MacDonald }

\examples{
library(affy)
data(affybatch.example)
plotDeg(affybatch.example, sampleNames(affybatch.example))
plotHist(affybatch.example, sampleNames(affybatch.example))
}
\keyword{ hplot }% at least one, from doc/KEYWORDS
