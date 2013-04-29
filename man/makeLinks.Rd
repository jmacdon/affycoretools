\name{entrezLinks}
\alias{makeLinks}
\alias{affyLinks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Add links to data for output with ReportingTools}
\description{
  This function is intended to be used for adding links to data that
  will be published using the ReportingTools package. The function
  currently adds links to Affymetrix probeset IDs and Entrez Gene IDs.
}
\usage{
  entrezLinks(df, ...)
  affyLinks(df, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{A data.frame, usually created by a call to
    \code{topTable}. This function is also used within
    \code{vennSelect2}, which creates such a data.frame.}
  \item{\dots}{Requires that all variables be named.}
}

\details{
  These functions are not designed to be called directly. Instead, they
  are used as part of the ReportingTools pipeline, to create links to
  either the Entrez Gene or Affymetrix website for data in an HTML
  table. See the basic ReportingTools vignette, specificially the part
  about .modifyDF
}

\value{
  Returns a \code{data.frame} with links added.
}

\author{James W. MacDonald <jmacdon@u.washington.edu>}
\keyword{hplot}% at least one, from doc/KEYWORDS
\keyword{manip}% __ONLY ONE__ keyword per line
