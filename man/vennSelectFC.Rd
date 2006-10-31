\name{vennSelectFC}
\alias{vennSelectFC}
\title{Select and Output Gene Lists Based on Venn Diagrams}
\description{
  This function is designed to output text and/or HTML tables based on
  the results of a call to \code{foldFilt}. The general idea being that
  one might want to create a Venn diagram showing probesets that are
  unique to particular comparisons, or consistent between comparisons,
  and then might want to output the probesets that are contained in each
  cell of the Venn diagram.
}
\usage{
vennSelectFC(eset, x, comps, order.by = "sum", method =
"same", text = TRUE, html = TRUE, ...)
}
\arguments{
  \item{eset}{A \code{\link[Biobase]{ExpressionSet}} object.}
  \item{x}{An indicator matrix showing up or down regulation based on
    fold change, usually from a call to \code{foldFilt}. See details for
    more information.}
  \item{comps}{A list containing all the comparisons to be
    made. Each list item should be a vector of length two. This should
    be identical to the 'comps' argument used in the call to
    \code{foldFilt} See details for more information.}
  \item{order.by}{One of 'sum', 'max', 'median', or 'mean'. This orders
    the output for those tables that have multiple fold change values
    based on the summary statistic chosen. Defaults to 'sum'.}
   \item{method}{ One of "same", "both", "up", "down", "sameup", or
     "samedown". See details for more information.}
   \item{text}{Boolean. Output text tables? Defaults to \code{TRUE}}
   \item{html}{Boolean. Output HTML tables? Defaults to \code{TRUE}}
   \item{...}{Used to pass other arguments to \code{probes2table}, in
    particular, to change the argument to \code{anncols} 
    which controls the columns of hyperlinks to online databases (e.g.,
    Entrez Gene, etc.). See \code{\link[annaffy]{aaf.handler}} for more
    information.}
}
\details{
  The purpose of this function is to output the probesets listed in a
  Venn diagram that has been produced by a call to \code{foldFilt}. A
  small example would be as follows:

  Assume an \code{ExpressionSet} exists that contains expression values
  for three Affymetrix chips, say a control, and two experimentals. One
  might want to know what probesets are different between each of the
  experimentals and the control, and those that are different between
  both of the experimentals and the control. We first make the
  comparisons, based on a fold change of 2 (or a difference of 1 on the
  log scale).

  comps <- list(c(1,2), c(1,3))

  This list indicates what comparisons we want. In this case 1vs2 and
  1vs3.

  out <- foldFilt(eset, fold = 1, groups = 1:3, comps = comps,
  compnames=c("Control vs experimental1", "Control vs experimental2"),
  save = TRUE)

  By setting save = TRUE, we are saving a list, the first item being a
  vector of the number of probesets in each comparison, the second item
  being an indicator matrix showing up or down regulation based on a
  two-fold difference. We could make a Venn diagram using this matrix
  with \code{vennCounts2} and \code{\link[limma]{vennDiagram}}. If we
  then wanted to output the probesets in each cell of that Venn diagram,
  we could use \code{vennSelectFC} as follows:

  vennSelectFC(eset, out[[2]], comps)

  One thing to note here is that the names of the resulting tables as
  well as the columns containing the fold change values will be
  extracted from the column names of the indicator matrix. This matrix
  will get its column names from the 'compnames' argument to
  \code{foldFilt}, so it is best to use reasonable names here. Also note
  that any character used in the 'compnames' argument that is not a
  valid character for a file name will be stripped out.
}
\value{
  Called only for the side effect of outputting HTML and/or text tables.
}
\author{James W. MacDonald <jmacdon@med.umich.edu>}
\keyword{manip}
