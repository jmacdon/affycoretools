\name{makeFCList}
\alias{makeFCList}
\title{Make a Fold Change List}
\description{
  This function is used to create a \code{list} containing named lists of fold
  change values suitable to be passed to \code{probes2table}, and
  ultimately \code{\link[annaffy]{aafTableAnn}} in the annaffy
  package. This is an internal function and is not intended to be called
  by end users.
}
\usage{
makeFCList(x, indices, comps, dat)
}
\arguments{
  \item{x}{An indicator matrix showing up or down regulation based on
    fold change, usually from a call to \code{foldFilt}. See details for
    more information.}
  \item{indices}{A \code{list} created by a call to \code{makeIndices}.}
  \item{comps}{A list containing all the comparisons to be
    made. Each list item should be a vector of length two. This should
    be identical to the 'comps' argument used in the call to
    \code{foldFilt}}
  \item{dat}{A \code{matrix} containing the expression values, usually
    extracted using the \code{exprs} accessor.}
}
\details{
  The purpose of this function is to create a list of named lists that
  will be suitable to be passed to \code{probes2table}, which will
  result in a column of fold change values with a reasonable name when
  output as an HTML or text table using annaffy.
}
\author{James W. MacDonald <jmacdon@med.umich.edu>}
\keyword{internal}
