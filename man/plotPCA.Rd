\name{plotPCA}
\alias{plotPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A Function to Make a PCA Plot from an exprSet }
\description{
 This function makes a PCA plot from an exprSet
}
\usage{
plotPCA(eset, groups, groupnames, addtext = NULL, x.coord = NULL,
y.coord = NULL, screeplot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{ An \code{exprSet} }
  \item{groups}{A numeric vector delineating group membership for samples. }
  \item{groupnames}{A character vector describing the different
    groups. }
  \item{addtext}{A character vector of additional text to be placed just
    above the plotting symbol for each sample. This is helpful if there
    are a lot of samples for identifying e.g., outliers.}
  \item{x.coord}{Pass an x-coordinate if automatic legend placement fails}
  \item{y.coord}{Pass a y-coordinate if automatic legend placement
    fails}
  \item{screeplot}{Boolean: Plot a \code{\link[stats]{screeplot}}
    instead of a PCA plot? Defaults to \code{FALSE}}
  
}

\value{
 This function returns nothing. It is called only for the side effect of
 producing a PCA plot or screeplot.
}
\author{ James W. MacDonald}
\keyword{hplot }% at least one, from doc/KEYWORDS
