\name{getUniqueLL}
\alias{getUniqueLL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Subset a Vector of Probesets}
\description{
  This function will take a vector of Affy IDs and return a vector of
  Entrez IDs that have replicated IDs removed. The resulting vector will
  still have the corresponding Affy IDs appended as names, which is
  important for some functions.
}
\usage{
getUniqueLL(probes, annot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probes}{ A vector of probe IDs}
  \item{annot}{The annotation package for the chip used}
}
\details{
 Subsetting a set of Affy IDs to unique Entrez Gene IDs is a common
 thing to do prior to doing a hypergeometric test. Functions such as
 \code{\link[Category]{hyperGtest}} can use un-named vectors of Entrez
 IDs (e.g., unique(getLL(probeIDs, annot))), but there is some
 functionality that requires the Entrez Gene IDs to be in a named
 vector, with the names being the associated Probeset IDs.

 As an example, \code{hyperGoutput} will only work correctly if the
 input Entrez ID vector is named with the associated Probeset IDs.
}
\value{
 A named vector of unique Entrez IDs
}

\author{James W. MacDonald <jmacdon@med.umich.edu>}

\keyword{manip}
