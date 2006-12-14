\name{houtSel}
\alias{houtSel}
\alias{houtAll}
\alias{houtSplit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Internal functions for hyperGoutput}
\description{
 These functions are internal functions for \code{hyperGoutput} and are
 not intended to be called directly by the end user.
}
\usage{
houtSel(tab, fit, comp, statistics, html, text)
houtAll(tab, eset, fit, comp, statistics, html, text)
houtSplit(tab, eset, fit, comp, statistics, html, text)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{ A \code{data.frame}, resulting from a call to
    \code{\link[GOstats]{probeSetSummary}}}
  \ietm{eset}{An \code{link[Biobase]{ExpressionSet}} object}
  \item{fit}{An \code{\link[limma]{MArrayLM}} object, resulting from a
    call to \code{\link[limma]{eBayes}}}
  \item{comp}{Numeric vector of length one, designating which comparison
    to select from the \code{\link[limma]{MArrayLM}} object}
  \item{statistics}{Any or all of 'tstat', 'pval', 'FC', used to control
    the output of statistics from the \code{\link[limma]{MArrayLM}} object} 
  \item{html}{Boolean. Output HTML tables? Defaults to \code{TRUE}}
  \item{text}{Boolean. Output text tables? Defaults to \code{TRUE}}
}
\details{
  These functions control the output and ordering of the probesets for
  the \code{hyperGoutput} function.}
}
\value{
 Nothing is returned. These functions are called solely for their side
 effect, the output of HTML and/or text tables.
}
\author{ James W. MacDonald <jmacdon@med.umich.edu>}

\keyword{internal }% __ONLY ONE__ keyword per line
