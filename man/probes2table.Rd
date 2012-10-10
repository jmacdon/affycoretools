\name{probes2table}
\alias{probes2table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Convert Affy Probe ids to Annotated HTML Table }
\description{
 A function to convert a vector of Affy ids to an annotated HTML table.
}
\usage{
probes2table(eset, probids, lib, otherdata = NULL,
 anncols = aaf.handler()[c(1:3, 6:7, 9:12)], html = TRUE, text = FALSE,
 express = TRUE, save = FALSE, filename, title = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{ An \code{ExpressionSet} containing Affy expression
  values.}
  \item{probids}{ A vector of probe ids. }
   \item{lib}{ An annotation package for the Affy chips used. }
   \item{otherdata}{A *named* list of additional information to include
     in the resulting table. Examples would be t-statistics, p-values,
     fold change, etc. Each list item should be a vector the same length
     as the probids vector. The name associated with each list item will
     be used as the column name in the resulting table.}
  \item{anncols}{ A vector of things to annotate,
    produced by a call to aaf.handler().}
  \item{html}{ Output data in HTML tables? Defaults to \code{TRUE}. }
  \item{text}{ Output data in text tables? Defaults to \code{TRUE}. }
  \item{express}{ Output expression values in table?
    Defaults to \code{TRUE}. }
  \item{save}{ Should tables be saved as R objects for
    further processing? Defaults to \code{FALSE}. }
  \item{filename}{ Filename of the resulting HTML table.}
  \item{title}{Title for HTML table. If \code{NULL}, the filename will
    be used.}
}

\value{
 If \code{save} is \code{TRUE}, a \code{data.frame} is saved containing
 the data.
}
\author{ James W. MacDonald <jmacdon@u.washington.edu> }

\seealso{ \code{\link[limma:toptable]{topTable}},
  \code{\link[annaffy]{aafTableAnn}} } 

\keyword{ manip }% at least one, from doc/KEYWORDS
