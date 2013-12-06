\name{getOrder}
\alias{getOrder}
\title{Order Probesets Based on Fold Change Values}
\usage{
  getOrder(fcList, order.by)
}
\arguments{
  \item{fcList}{A \code{list} of named lists, generally
  produced by a call to \code{makeFCList}}

  \item{order.by}{One of 'sum', 'max', 'median', or 'mean'.
  This orders the output for those tables that have
  multiple fold change values based on the summary
  statistic chosen. Defaults to 'sum'.}
}
\description{
  This function is used to create a \code{list} containing
  the correct order to output a set of probesets based on
  the fold change values for those probesets. This is an
  internal function and is not intended to be called by end
  users.
}
\details{
  The purpose of this function is to come up with some
  ordering for probesets that will be ouput by a call to
  \code{vennSelectFC}. The ordering for probesets unique to
  a given comparison is simple, since there is only one
  fold change to be considered. However, when there are
  more than one comparison, we have to decide how the
  probesets should be ordered. The choices available
  include those that seemed reasonable to me when I wrote
  this function. Presumably there are other reasonable
  choices that may need to be included.
}
\author{
  James W. MacDonald <jmacdon@u.washington.edu>
}
\keyword{internal}

