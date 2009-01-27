\name{vennSelectBM}
\alias{vennSelectBM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Select and Output Genelists Based on Venn Diagrams using biomaRt }
\description{
  This function is designed to output HTML and text tables based on
  the results of a call to \code{\link[limma]{decideTests}}. This
  function is very similar to \code{vennSelect}, except it uses the
  \code{biomaRt} package to annotate genes, and the \code{annotate}
  package to create the HTML table.
}
\usage{
vennSelectBM(eset, design, x, contrast, fit, method = "same", adj.meth = "BH",
stat = "fstat", otherstats = c("pval", "FC"), order.by = "pval",
foldFilt = NULL, save = FALSE, species, links = linksBM()[1:3],
otherdata = annBM()[1:3], ann.source = "entrezgene",
html = TRUE, text = TRUE, affyid = FALSE, ...)
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
   \item{fit}{An \code{\link[limma]{MArrayLM}} object, from a call to
     \code{\link[limma]{eBayes}}.}
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
  \item{species}{The species name. This must be in a particular format
    for biomaRt. An example for human is "hsapiens", or for mouse
    "mmusculus".}
  \item{links}{A character vector of things to annotate with hyperlinks
    to online databases. See \code{linksBM} for possible values.}
  \item{otherdata}{A character vector of things to annotate with text
    only (i.e., no hyperlinks). See \code{annBM} for possible values.}
  \item{ann.source}{The annotation source of the IDs that will be used
    to annotate the genes. The default value is "entrezgene". See
    details for other possibilities.}
  \item{html}{Boolean. Output HTML tables? Defaults to \code{TRUE}}
  \item{text}{Boolean. Output text tables? Defaulst to \code{TRUE}}
  \item{affyid}{Boolean. Are the IDs used to annotate these data
    Affymetrix IDs?}
  \item{...}{Used to pass other variables to e.g., \code{htmlpage}.}
}
\details{
  The purpose of this function is to output HTML tables with
  lists of genes that fulfill the criteria of a call to
  \code{\link[limma]{decideTests}} as well as the direction of
  differential expression.

  The IDs that will be used to annotate the genes depend on the source
  of the data. If, for example, one is using an Affymetrix chip that
  doesn't have a BioC annotation package, then the IDs will be
  Affymetrix IDs. To find out the correct name to use for the ann.source
  argument, one can create a connection to a Biomart database using
  \code{\link[biomaRt]{useMart}} and then deduce the correct argument by
  the output from \code{listFilters(mart)}. It will usually be
  something starting with 'affy', and contain the name of the chip.

  If one is using one of the re-mapped CDFs from MBNI at University of
  Michigan, then the IDs to use depend on the mapping used to create the
  CDF. At this time, only three types of CDFs can be used; EntrezGene,
  UniGene, and RefSeq. One can determine the correct ann.source argument
  by creating a connection to a Biomart database, and then calling
  \code{linksBM(mart, linksBM())[[3]]}. 

  Some important things to note: First, the names of the HTML 
  tables are extracted from the \code{colnames} of the
  \code{TestResults} object, which come from the contrasts matrix, so
  it is important to use something descriptive. Second, the method
  argument is analogous to the \code{include} argument from
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

  Calling the function normally will result in the output of HTML
  tables. 

  Calling the function with save set to \code{TRUE} will output
  HTML tables as well as a vector of counts for each
  comparison. This is useful when using the function programmatically
  (e.g., when making reports using Sweave).

  out <- vennSelectBM(eset, fit, design, x, <other arguments>,  save = TRUE)

  An alternative would be to use \code{vennCounts2} and
  \code{\link[limma]{vennDiagram}} to output a Venn diagram, which is
  probably more reasonable since the tables being output are supposed to
  be based on a Venn diagram.
}
\value{
 Normally called only for the side effect of producing HTML
 tables. However, setting save to \code{TRUE} will output a vector of
 counts that can be used for making Sweave-style reports.
 }
\author{James W. MacDonald <jmacdon@med.umich.edu>}

\keyword{ manip}% at least one, from doc/KEYWORDS
