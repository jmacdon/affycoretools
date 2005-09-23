\name{probes2table}
\alias{probes2table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Convert Affy Probe ids to Annotated HTML Table }
\description{
 A function to convert a vector of Affy ids to an annotated HTML table.
}
\usage{
probes2table(eset, probids, lib, FC = NULL, tstat = NULL, pval = NULL,
 anncols = aaf.handler()[c(1:3, 7:8, 10:13)], html = TRUE, text = FALSE,
 express = TRUE, save = FALSE, filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{ An \code{exprSet} containing Affy expression
  values.}
  \item{probids}{ A vector of probe ids. }
   \item{lib}{ An annotation package for the Affy chips used. }
  \item{FC}{ Output fold changes in table? Defaults to \code{FALSE}. }
  \item{tstat}{ Output t-statistics in table? Defaults to \code{FALSE}. }
  \item{pval}{ Output (adjusted) p-values in table?
    Defaults to \code{FALSE}. }
  \item{anncols}{ A vector of things to annotate,
    produced by a call to aaf.handler().}
  \item{html}{ Output data in HTML tables? Defaults to \code{TRUE}. }
  \item{text}{ Output data in text tables? Defaults to \code{TRUE}. }
  \item{express}{ Output expression values in table?
    Defaults to \code{TRUE}. }
  \item{save}{ Should tables be saved as R objects for
    further processing? Defaults to \code{FALSE}. }
  \item{filename}{ Filename of the resulting HTML table.}
}

\value{
 If \code{save} is \code{TRUE}, a \code{data.frame} is saved containing
 the data.
}
\author{ James W. MacDonald <jmacdon@med.umich.edu> }

\seealso{ \code{\link[limma]{topTable}}, \code{\link[annaffy]{aafTableAnn}} }

\keyword{ manip }% at least one, from doc/KEYWORDS
