\name{vennSelect}
\alias{vennSelect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Select and Output Genelists Based on Venn Diagrams }
\description{
  This function has two purposes. The main purpose is to select the
  lists of genes from each section of a Venn diagram and output in
  HTML and text tables. The second purpose is to produce Venn diagrams
  where the criteria for the intersections are that the genes pass a
  certain significance level as well as being differentially expressed
  in the same direction.
}
\usage{
vennSelect(eset, design, x, contrasts,  method = "same", foldFilt = NULL,
fit = NULL, indices.only = FALSE, save = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{ A \code{\link[Biobase]{exprSet}} object.}
  \item{design}{ A \code{design} matrix, usually from a call to
    \code{model.matrix}. See details for more information.}
  \item{x}{A \code{\link[limma]{TestResults}} object, usually from a
    call to \code{\link[limma]{designTests}}.}
  \item{contrasts}{A contrasts matrix, produced either by hand, or by a
    call to \code{\link[limma]{makeContrasts}}}
  \item{method}{ One of "same", "both", "up", "down", "sameup", or
    "samedown". See details for more information.}
  \item{foldFilt}{A fold change to use for filtering. Default is \code{NULL},
    meaning no filtering will be done.}
   \item{fit}{An \code{\link[limma]{MArrayLM}} object, from a call to
     \code{\link[limma]{eBayes}}. The default is \code{NULL}. This is
     only required if \code{foldFilt} is not \code{NULL}, as the fold
     change data are extracted from this object.}
  \item{indices.only}{Boolean - Should just the indices indicating genes
    that fulfill the criteria be output?}
  \item{save}{Boolean - If \code{TRUE}, output a count of genes that
    fulfill the criteria. Useful for e.g., Sweave-type reports.}
  \item{...}{Used to pass other arguments to \code{probes2table}, in
    particular, to change the argument to \code{anncols} 
    which controls the columns of hyperlinks to online databases (e.g.,
    Entrez Gene, etc.). See \code{\link[annaffy]{aaf.handler}} for more information.}
}
\details{
  As mentioned in the description, this function has multiple
  purposes. The main purpose is to output HTML and text tables with
  lists of genes that fulfill the criteria of a call to
  \code{\link[limma]{decideTests}} as well as the direction of
  differential expression.

  Some important things to note: First, the names of the HTML and text
  tables are extracted from the \code{colnames} of the
  \code{TestResults} object, which come from the contrasts matrix,
  it is important to use something descriptive. Second, the method is
  analogous to the \code{include} argument from
  \code{\link[limma]{vennCounts}} or
  \code{\link[limma]{vennDiagram}}. Choosing "both" will select genes
  that are differentially expressed in one or more comparisons,
  regardless of direction. Choosing "up" or "down" will select genes
  that are only differentially expressed in one direction. Choosing
  "same" will select genes that are differentially expressed in the same
  direction. Choosing "sameup" or "samedown" will select genes that are
  differentially expressed in the same direction as well as 'up' or
  'down'.

  Note that this is different than sequentially choosing "up"
  and then "down". For instance, a gene that is upregulated in one
  comparison and downregulated in another comparison will be listed in
  the intersection of those two comparisons if "both" is chosen, it will
  be listed in only one comparison for both the "up" and "down" methods,
  and it will be listed in the union (e.g., not selected) if "same" is
  chosen.

  Calling the function normally will result in the output of HTML and
  text tables:

  vennSelect(eset, fit, design, x)

  Calling the function with indices.only set to \code{TRUE} will result
  in a list of \code{TRUE} and \code{FALSE} values. This is normally
  only used by \code{vennCounts2} to output a
  \code{\link[limma]{VennCounts}} object, which can then be used to plot
  a Venn Diagram.

  vennSelect(x, indices.only = TRUE)

  Calling the function with save set to \code{TRUE} will output both
  HTML and text tables as well as a vector of counts for each
  comparison. This is useful when using the function programmatically
  (e.g., when making reports using Sweave).

  out <- vennSelect(eset, fit, design, x, save = TRUE)
}
\value{
 Normally called only for the side effect of producing HTML and text
 tables. However, setting indices.only to \code{TRUE} will output a list of \code{TRUE}
 and \code{FALSE} values (normally used by \code{vennCounts2}), and
 setting save to \code{TRUE} will output a vector of counts.
 }
\author{James W. MacDonald <jmacdon@med.umich.edu>}

\keyword{ manip}% at least one, from doc/KEYWORDS
