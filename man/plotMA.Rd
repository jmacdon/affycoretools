\name{plotma}
\alias{plotma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A Function to Make an MA Plot }
\description{
 This function makes an MA plot from either a matrix or two vectors. 
}
\usage{
plotma(x, y = NULL, main = "", log = TRUE, out = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A matrix of two columns or a vector of numbers. If a matrix
    with more than two columns is passed, only the first two will be used.}
  \item{y}{A vector of numbers. If \code{NULL}, \code{x} should be a
    matrix with two columns.}
  \item{main}{A title for the plot. }
  \item{log}{Boolean describing if the data are log2 transformed.}
  \item{out}{Boolean describing if the M and A values should be returned.}
}

\value{
 If \code{out} is \code{TRUE}, a matrix of M and A values is returned.
}
\author{ James W. MacDonald}
\keyword{hplot }% at least one, from doc/KEYWORDS
