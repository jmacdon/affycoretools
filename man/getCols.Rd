\name{getCols}
\alias{getCols}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Correct Ordering of Contrasts}
\description{
 A function to determine the correct ordering of contrasts for
 \code{vennSelect}. This is an internal function and should not be
 called by the end user
}
\usage{
getOrd(design, contrasts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design}{A design matrix, usually produced by a call to
    \code{model.matrix}.}
  \item{contrasts}{A contrasts matrix, either produced by hand or as a
    result of a call to \code{\link[limma]{makeContrasts}}}
}

\value{
 This function returns the correct order for the contrasts.
}

\author{James W. MacDonald}
\keyword{ internal }% __ONLY ONE__ keyword per line
