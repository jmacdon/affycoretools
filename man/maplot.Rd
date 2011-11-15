\name{maplot}
\alias{maplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A Function to make MA plots from all arrays. }
\description{
 This function creates an MA plot for all arrays in either an
 ExpressionSet or a matrix. A 'baseline' array is created using the
 median expression for each gene, and each array is then compared to the
 baseline array.
}
\usage{
maplot(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An ExpressionSet or matrix containing log-transformed
    array data.}
}

\value{
 No output. Used only for the side effect of creating MA plots.
}
\author{ James W. MacDonald <jmacdon@med.umich.edu>}
\keyword{hplot }% at least one, from doc/KEYWORDS
