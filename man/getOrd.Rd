\name{getOrd}
\alias{getOrd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Correct Ordering of Contrasts}
\description{
 A function to determine the correct ordering of contrasts for
 \code{vennSelect}. This is an internal function and should not be
 called by the end user
}
\usage{
getOrd(x, design)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A \code{TestResults} object from a call to \code{decideTests}.}
  \item{design}{A design matrix, usually produced by a call to
    \code{model.matrix}.}
}

\value{
 This function returns the correct order for the contrasts.
}

\author{James W. MacDonald <jmacdon@med.umich.edu>}
\keyword{ internal }% __ONLY ONE__ keyword per line
