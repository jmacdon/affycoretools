\name{convert.back}
\alias{convert.back}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A Function to convert M and A values back to R and G values }
\description{
 This function takes M and A values and converts back to log2 R and G
 values.
}
\usage{
convert.back(M, A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{A vector of log2 R/G values.}
  \item{A}{A vector of log2 R + G values.}
}

\value{
 A matrix of log2 R and G values.
}
\author{ James W. MacDonald}
\keyword{hplot }% at least one, from doc/KEYWORDS
