\name{annBM}
\alias{annBM}
\alias{linksBM}
\title{Select Available Annotation from a Biomart }
\description{
  These functions are designed to do two things that are useful for an end
  user. If called with no arguments, they will output a character vector
  of annotation sources that are typically available from a Biomart
  database. If called with a 'mart' connection (typically created by a
  call to \code{\link[biomaRt]{useMart}}), they will return a character
  vector of annotation sources that exist for that particular Biomart
  and species. If called with a 'mart' connection and a character vector
  of annotation sources, they will return a list that is intended to be
  used by other functions for creating HTML pages. This last function
  doesn't have any real utility for the end user.
}
\usage{annBM(mart, annot, species)
linksBM(mart, annot, affyid = FALSE, ann.source = NULL)
}
\arguments{
  \item{mart}{A 'mart' connection, typically created by a
    call to \code{\link[biomaRt]{useMart}}.}
  \item{annot}{A character vector of annotation sources. This is not
    typically useful for an end user to specify.}
  \item{affyid}{Boolean. Are the IDs being annotated Affymetrix IDs?}
  \item{ann.source}{Character. The filter name to use if we are
    annotating Affymetrix IDs.}
  \item{species}{A species name, of the form e.g., 'hsapiens'}
}
\details{
  The purpose of these functions is to either give an example of typical
  annotation sources that may be available at a particular Biomart, or
  to output those sources that are known to exist at a Biomart.

  \code{linksBM} is intended to list those annotation sources that may
  be turned into hyperlinks whereas \code{annBM} is intended to list
  those annotation sources that will not be linked.

  These functions have only a few of the possible annotation sources, and
  currently there is no simple way to extend these sources. Additions to
  the list are possible, however. Please contact me if there is
  something in particular that should be included in either list.
}
\value{
  Normally called by an end user to output a character vector of
  annotation sources.
}
\author{
  James W. MacDonald <jmacdon@med.umich.edu>
}
\examples{
annBM()
}
\keyword{manip}
