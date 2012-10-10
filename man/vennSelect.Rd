\name{vennSelect}
\alias{vennSelect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Select and Output Genelists Based on Venn Diagrams }
\description{
  This function is designed to output text and/or HTML tables based on
  the results of a call to \code{\link[limma]{decideTests}}.
}
\usage{
vennSelect(eset, design, x, contrast, fit, method = "same", adj.meth =
"BH", stat = "fstat", otherstats = c("pval","FC"), order.by = "pval",
foldFilt = NULL, save = FALSE, titleadd = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{ A \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} object.}
  \item{design}{ A \code{design} matrix, usually from a call to
    \code{model.matrix}. See details for more information.}
  \item{x}{A \code{\link[limma]{TestResults}} object, usually from a
    call to \code{\link[limma]{decideTests}}.}
  \item{contrast}{A contrasts matrix, produced either by hand, or by a
    call to \code{\link[limma]{makeContrasts}}}
   \item{fit}{An \code{\link[limma:marraylm]{MArrayLM}} object, from a call to
     \code{\link[limma:ebayes]{eBayes}}.}
  \item{method}{ One of "same", "both", "up", "down", "sameup", or
    "samedown". See details for more information.}
  \item{adj.meth}{Method to use for adjusting p-values. Default is 'BH',
    which corresponds to 'fdr'. Ideally one would set this value to be
    the same as was used for \code{\link[limma]{decideTests}}.}
  \item{stat}{The statistic to report in the resulting HTML
    tables. Choices are 'fstat', 'tstat', and \code{NULL}. Ideally, the
    statistic chosen would correspond to the method used in
    \code{\link[limma]{decideTests}}. In other words, if one used
    methods such as 'separate' or 'hierarchical', which are based on a
    t-statistic, one should choose 'tstat', however, if one used
    'nestedF', the logical choice would be 'fstat'.}
  \item{otherstats}{Other statistics to be included in the HTML
    tables. Choices include 'pval' and 'FC'.}
  \item{order.by}{Which statistic should be used to order the probesets?
    Choices include 'fstat', 'tstat', 'pval', and 'FC'. Note that if
    'FC' is chosen and there are more than one set of fold changes, the
    first will be used.}
  \item{foldFilt}{A fold change to use for filtering. Default is \code{NULL},
    meaning no filtering will be done.}
  \item{save}{Boolean - If \code{TRUE}, output a count of genes that
    fulfill the criteria. Useful for e.g., Sweave-type reports.}
  \item{...}{Used to pass other arguments to \code{probes2table}, in
    particular, to change the argument to \code{anncols} 
    which controls the columns of hyperlinks to online databases (e.g.,
    Entrez Gene, etc.). See \code{\link[annaffy]{aaf.handler}} for more
    information.}
  \item{titleadd}{Additional text to add to the title of the HTML
    tables. Default is NULL, in which case the title of the table will
    be the same as the filename.}
}
\details{
  The purpose of this function is to output HTML and text tables with
  lists of genes that fulfill the criteria of a call to
  \code{\link[limma]{decideTests}} as well as the direction of
  differential expression.

  Some important things to note: First, the names of the HTML and text
  tables are extracted from the \code{colnames} of the
  \code{TestResults} object, which come from the contrasts matrix, so
  it is important to use something descriptive. Second, the method
  argument is analogous to the \code{include} argument from
  \code{\link[limma:venn]{vennCounts}} or
  \code{\link[limma:venn]{vennDiagram}}. Choosing "both" will select genes
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

  Calling the function with save set to \code{TRUE} will output both
  HTML and text tables as well as a vector of counts for each
  comparison. This is useful when using the function programmatically
  (e.g., when making reports using Sweave).

  out <- vennSelect(eset, fit, design, x, save = TRUE)

  An alternative would be to use \code{vennCounts2} and
  \code{\link[limma:venn]{vennDiagram}} to output a Venn diagram, which is
  probably more reasonable since the tables being output are supposed to
  be based on a Venn diagram.
}
\value{
 Normally called only for the side effect of producing HTML and text
 tables. However, setting save to \code{TRUE} will output a vector of
 counts that can be used for making Sweave-style reports.
 }
\author{James W. MacDonald <jmacdon@u.washington.edu>}

\keyword{ manip}% at least one, from doc/KEYWORDS
