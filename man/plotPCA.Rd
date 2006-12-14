\name{plotPCA}
\alias{plotPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A Function to Make a PCA Plot from an ExpressionSet }
\description{
 This function makes a PCA plot from an ExpressionSet, exprSet or matrix
}
\usage{
plotPCA(eset, groups = NULL, groupnames = NULL, addtext = NULL, x.coord = NULL,
y.coord = NULL, screeplot = FALSE, squarepca = FALSE, pch = NULL, col =
NULL, pcs = c(1,2), legend = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{ An \code{ExpressionSet}, \code{exprSet} or matrix }
  \item{groups}{A numeric vector delineating group membership for
    samples. Default is \code{NULL}, in which case default plotting symbols and
    colors will be used.}
  \item{groupnames}{A character vector describing the different
    groups. Default is \code{NULL}, in which case the sample names will
    be used.}
  \item{addtext}{A character vector of additional text to be placed just
    above the plotting symbol for each sample. This is helpful if there
    are a lot of samples for identifying e.g., outliers.}
  \item{x.coord}{Pass an x-coordinate if automatic legend placement fails}
  \item{y.coord}{Pass a y-coordinate if automatic legend placement
    fails}
  \item{screeplot}{Boolean: Plot a \code{\link[stats]{screeplot}}
    instead of a PCA plot? Defaults to \code{FALSE}}
  \item{squarepca}{Should the y-axis of the PCA plot be made comparable
    to the x-axis? This may aid in interpretation of the PCA
    plot. Defaults to \code{FALSE}}
  \item{pch}{A numeric vector indicating what plotting symbols to
    use. Default is \code{NULL}, in which case default plotting symbols
    will be used. Note that this argument will override the 'groups'
    argument.}
  \item{col}{A numeric or character vector indicating what color(s) to
    use for the plotting symbols. Default is \code{NULL} in which case
    default colors will be used. Note that this argument will override
    the 'groups' argument.}
  \item{pcs}{A character vector of length two indicating which principal
    components to plot. Defaults to the first two principal components.}
  \item{legend}{Boolean. Should a legend be added to the plot? Defaults
    to \code{TRUE}}
  \item{...}{Further arguments to be passed to \code{plot}. See the help page
    for \code{plot} for further information.}
}

\value{
 This function returns nothing. It is called only for the side effect of
 producing a PCA plot or screeplot.
}
\author{ James W. MacDonald <jmacdon@med.umich.edu>}
\keyword{hplot }% at least one, from doc/KEYWORDS
