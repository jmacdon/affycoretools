\name{foldFilt}
\alias{foldFilt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Output Fold Change Data }
\description{
  This function is designed to take an \code{exprSet} and some
  comparisons and output either HTML tables, text files, or both.
}
\usage{
foldFilt(object, fold = 1, groups, comps, compnames, save = FALSE, text
= TRUE, html = TRUE, filterfun = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An \code{exprSet} object}
  \item{fold}{The log fold change cutoff to use. Note that this is log
    base two.}
  \item{groups}{A vector of group identifiers. Probably easiest to use a
    numeric vector}
  \item{comps}{A list containing all the comparisons to be
    made. Each list item should be a vector of length two. See details
    for more information.}
  \item{compnames}{A character vector of the names for each of the
    comparisons to be made. This will be the name of the resulting HTML
    or text file.}
  \item{save}{Boolean. If \code{TRUE}, a list will be returned. The
    first item in the list will be a vector showing the number of
    'significant' genes for each comparison. The second item will be a
    matrix of -1's, 0's and 1's indicating a significant difference, and
    the direction of the difference. The first item is useful for
    creating Sweave - based reports and the second is useful for making
    Vennn diagrams using the \code{vennDiagram} from the limma package.}
  \item{html}{Boolean - if \code{TRUE}, output HTML tables}
  \item{text}{Boolean - if \code{TRUE}, output text tables}
  \item{filterfun}{A filtering function to be by
    \code{\link[genefilter]{genefilter}} to filter the data using
    additional criteria. See details for more information}
}

\value{
  Returns a list; see above for the elements of the list. This function
  is mainly called for the side effect of outputting HTML or text files containing
  annotated 'significant' gene lists.
}
\details{
  This function is useful for outputting annotated gene lists for
  multiple fold change comparisons. The genes will be ordered by the
  absolute fold change. Note that this function is essentially a wrapper
  to call \code{annaffy}, so is only useful for Affymetrix GeneChips for
  which there is an annotation package.

  Without attaching a data file to this package, it is not possible to
  give a working example. Instead, here is a 'for instance'.

  Say you have an \code{exprSet} containing four Affy HG-U133Plus2
  chips. There is no replication, and you simply want to output genes
  with a two-fold or greater difference between the first chip and each
  of the last three (the first chip is the control, and the other three
  are experimentals). The \code{exprSet} is called eset.

  Additionally, say we don't want any genes called significant if both
  of the samples have very low expression. We can set up a filter using
  the \pkg{genefilter} package.

  f1 <- kOverA(1,6)
  filt <- filterfun(f1)

  foldFilt(eset, groups=1:4, comps=list(c(2, 1), c(3, 1), c(4, 1)),
  compnames=c("Expt1-Cont","Expt2-Cont","Expt3-Cont"), filterfun = filt)

  This will output three HTML tables called 'Expt1-Cont.html', etc.,
  each containing sorted genes that have two-fold or greater differences
  between the two samples.
  }
\author{James W. MacDonald <jmacdon@med.umich.edu> }

\keyword{manip }% __ONLY ONE__ keyword per line
