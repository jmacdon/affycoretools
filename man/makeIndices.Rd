\name{makeIndices}
\alias{makeIndices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create Indices for Venn Diagrams}
\description{
  This function is used to create indices for making Venn diagrams or
  \code{vennCounts} objects. This is an internal method and is not
  intended to be called by the end user. 
}
\usage{
makeIndices(x, method = "same")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A \code{TestResults} object from a call to \code{decideTests}.}
  \item{method}{One of "same", "up", "down", "both", indicating
    direction of differential expression. See details of
    \code{vennSelect} for more information.}
}
\value{
 A list containing \code{TRUE} and \code{FALSE} values, to be used by \code{vennSelect}.
}
\author{ James W. MacDonald <jmacdon@u.washington.edu>}
\keyword{internal }% __ONLY ONE__ keyword per line
