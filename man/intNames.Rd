\name{intNames}
\alias{intNames}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create Names for Venn Diagram Intersections}
\description{
  This function is designed to create the names for all the
  intersections of a Venn diagram based on the names of all the single
  comparisons. This is an internal function and is not intended to be
  called by the end user.
}
\usage{
intNames(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A \code{TestResults} object from a call to \code{decideTests}.}
}
\value{
 A vector of names.
}
\author{James W. MacDonald}
\keyword{ internal}
