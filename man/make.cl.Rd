\name{make.cl}
\alias{make.cl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A Function to Make a Classlabel Vector for Plotting }
\description{
 This function takes a vector of filenames and makes a classlabel vector
 for plotting functions of \code{affystart}. This is an internal
 function and is not intended to be called by the end user.
}
\usage{
make.cl(filenames)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filenames}{A vector of .cel filenames }
}

\value{
  \item{cl }{A vector of numbers to be used for plotting colors and line
  types}
 
}
\author{James W. MacDonald <jmacdon@med.umich.edu> }

\keyword{internal }% at least one, from doc/KEYWORDS
