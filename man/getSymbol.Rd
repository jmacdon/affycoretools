\name{getSymbol}
\alias{getSymbol}
\title{Check to see if the Symbol is Available at a Biomart}
\description{
  This function is designed to check to see if the 'symbol' annotation
  is available at a given Biomart. This is an internal function and is
  not intended to be called by the end user.
}
\usage{getSymbol(mart)}
\arguments{
  \item{mart}{A 'mart' connection, typically created by a
  call to \code{\link[biomaRt]{useMart}}.}
}
\details{
  This function is designed to query a Biomart to see if the 'symbol'
  annotation is available for a particular species. This is an internal
  function and is not intended to be called by the end user.
}
\value{
  Either the symbol attribute (typically something like hgnc\_symbol), or
  'nosymbolavailable'.
}
\author{James W. MacDonald <jmacdon@med.umich.edu>}
\keyword{internal}
