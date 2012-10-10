\name{writeFit}
\alias{writeFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to output annotated fit data from limma}
\description{
  This function is designed to take an \code{ExpressionSet} an
  annotation package and an \code{lmFit} object, and output an annotated
  text file containing t-statistics, p-values, and fold change data for
  all contrasts.}
\usage{
writeFit(fit, annotation, eset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{A \code{lmFit} object, created by the limma package.}
  \item{annotation}{An annotation package, specific for the chip used in
    the analysis.}
  \item{eset}{An \code{ExpressionSet} object containing expression
    values.}
}
\details{
  This function is designed to output annotation data as well as
  statistics (p-values, fold change, t-statistics) for all probes on a
  chip.
}
\value{
  A \code{data.frame} is returned.
}
\author{James W. MacDonald <jmacdon@u.washington.edu>}
\seealso{\code{\link[limma:write.fit]{write.fit}}
}
\keyword{manip}

