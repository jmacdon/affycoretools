\name{hyperG2annaffy}
\alias{hyperG2annaffy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{HTML tables from GOIDs}
\description{
  Output HTML tables containing the 'enriched' genes for each GO term
  resulting from a call to \code{hyperGtable}.
}
\usage{
hyperG2annaffy(probids, lib, eset, fit = NULL, subset = NULL, comp = 1, type = "MF",
pvalue = 0.05, min.count = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probids}{ A vector of Affymetrix probe IDs}
  \item{lib}{ An annotation package (e.g., \code{hgu95av2})}
  \item{eset}{An \code{ExpressionSet}}
  \item{fit}{ An \code{\link[limma]{lmFit}} object. Only necessary if
    statistics are desired in the resulting table. Defaults to \code{NULL}.}
  \item{subset}{A numeric vector used to select GO terms to output (see
    description for more information). Defaults to \code{NULL}}
  \item{comp}{Which contrast/parameter estimate should be used to
    extract the relevant statistics? Only used if \code{fit} is not
    \code{NULL}. See description for more information.}
  \item{type}{One of "MF", "CC", "BP", indicating molecular function,
    cellular component, or biological process, respectively.}
  \item{pvalue}{The significance level used to choose GO terms}
  \item{min.count}{ The minimum number of a given GO term that must be
    on the chip in order to choose that GO term. This protects against
    very low p-values that result from the situation where there are
    very few genes with a given GO term on the chip, but one or two are
    found in the set of significant genes.}
}
\details{
 This function is used to create HTML tables based on the output of
 \code{hyperGtable}. The basic idea is as follows; as part of an
 analysis, say \code{hyperGtable} was used to create a table of
 'enriched' GO terms. Unfortunately, the table only lists GO terms and
 the number of probesets that are annotated to those GO terms, and the
 client may be interested in knowing what probesets are enriched for
 each (or some) GO term.

 The default behaviour is to output an HTML table for each GO term,
 containing the probesets that are annotated at that term (and that are
 in the set of significant genes). If only some of the GO terms are of
 interest, one may use the \code{subset} argument to select only
 particular rows. In addition, if the relevant t-statistics, p-values
 and fold changes are of interest, one can also use the \code{fit}
 argument to point to an \code{\link[limma]{lmFit}} object that contains
 these data, as well as the \code{comp} argument to indicate which
 parameter or contrast to use. Note that the \code{comp} argument
 defaults to 1, so the first parameter or contrast will be extracted by default.
}
\value{
 This function is used only for the side effect of creating HTML tables.
}
\author{ James W. MacDonald <jmacdon@u.washington.edu>}
\keyword{ manip }% __ONLY ONE__ keyword per line
\keyword{htest}
