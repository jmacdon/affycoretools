\name{getans2}
\alias{getans2}
\alias{getfldfil}
\alias{getpfil}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Interactive Functions for limma2anaffy }
\description{
  These are internal functions used to filter the number of genes for
  \code{limma2annaffy} as well as to automagically get required
  annotation packages for \code{annaffy}.
}
\usage{
getans2(msg, allowed = c("a", "c"))
getfldfil()
getpfil()
}
\arguments{
  \item{msg}{A character string to pass to the end user.}
  \item{allowed}{A character vector of acceptable responses.}
}

\details{
  These functions were added primarily for the error checking they
  afford. Now the end user can filter the number of genes based on
  p-value and fold change, and incorrect values entered will not result
  in unexpected behaviour.
}

\author{James W. MacDonald <jmacdon@med.umich.edu> }

\keyword{ internal }
