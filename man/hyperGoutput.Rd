\name{hyperGoutput}
\alias{hyperGoutput}
\title{ Output Tables Based on Hypergeometric Test}
\description{
  This function will output various tables containing probesets that are
  annotated to a particular GO, KEGG, or PFAM term. The tables are
  based on the results from a call to hyperGtest.
}
\usage{
hyperGoutput(hyptObj, eset, pvalue, categorySize, sigProbesets, fit = NULL,
subset = NULL,comp = 1, output = c("significant", "all", "split"),
statistics = c("tstat", "pval", "FC"), html = TRUE, text = TRUE, ...)
}

\arguments{
  \item{hyptObj}{A \code{HyperGResult} object, usually produced by a call to
    \code{\link[Category]{hyperGTest}}} 
  \item{eset}{An \code{ExpressionSet} object}
  \item{pvalue}{The p-value cutoff used for selecting significant GO
    terms. If not specified, it will be extracted from the
    \code{HyperGResult} object}
  \item{categorySize}{Number of terms in the universe required for a
    term to be significant. See details for more information}
  \item{sigProbesets}{Vector of probeset IDs that were significant in
    the original analysis.} 
  \item{fit}{An \code{\link[limma:marraylm]{MArrayLM}} object, produced from a
    call to \code{\link[limma:ebayes]{eBayes}}}
  \item{subset}{Numeric vector used to select particular tables to
    output. The default is to output tables for all terms. See details
    for more information}
  \item{comp}{Numeric vector of length one, used to indicate which
    comparison in the \code{\link[limma:marraylm]{MArrayLM}} object to use for
    extracting relevant statistics. See details for more information}
  \item{output}{One of 'selected', 'all', or 'split'. See details for
    more information}
  \item{statistics}{Which statistics to output in the resulting
    tables. Choices include 'tstat', 'pval', or 'FC', corresponding to
    t-statistics, p-values, and fold change, respectively}
  \item{html}{Boolean. Output HTML tables? Defaults to \code{TRUE}}
  \item{text}{Boolean. Output text tables? Defaults to \code{TRUE}}
  \item{\dots}{Allows end user to pass further arguments. The most
    notable would be an \code{anncols} argument, passed to
    \code{probes2table} to control the hyperlinked annotation
    columns. See \code{\link[annaffy]{aaf.handler}} for more information}
}
\details{
 This function is designed to be used to output the results from a
 hypergeometric test for over-represented terms. This function would be
 used at the end of an analysis such as:

 1.) Compute expression values
 2.) Fit a model using \code{limma}
 3.) Output significant probesets using \code{limma2annaffy}
 4.) Perform hypergeometric test using
 \code{\link[Category:HyperGResult-accessors]{hyperGTest}}

 At step 4, one can output a list of the over-represented terms using
 \code{\link[Category:HyperGResult-accessors]{htmlReport}}. One might then be interested in
 knowing which probesets contributed to the significance of a particular
 term, which is what this function is designed to do.

 One argument that can be passed to \code{\link[Category:HyperGResult-accessors]{htmlReport}}
 (and also to \code{hyperGoutput}) is \code{categorySize}, which gives a
 lower bound for the number of probesets with a particular term in the
 universe. In other words, assume that a particular GO term is annotated
 to three probesets on a given chip. If, after doing a t-test to detect
 differentially expressed probesets, one of those probesets were found
 to be significantly differentially expressed and was then used to do a
 hypergeometric test, that GO term would be significant, with a small
 p-value. However, this is probably not very strong evidence that the GO
 term is actually over-represented, since there were only three to begin
 with. By setting \code{categorySize} to a sensible value (such as 10),
 this situation can be avoided.

 This function will output HTML and/or text tables containing annotation
 information about each probeset as well as the expression values. In
 addition, if limma were used to fit the model, the relevant statistics
 (t-statistic, p-value, fold change) can also be output in the table by
 passing the \code{\link[limma:marraylm]{MArrayLM}} object that resulted
 from a 
 call to \code{\link[limma:ebayes]{eBayes}}. The \code{statistics}
 argument can 
 be used to control which statistics are output.

 By default \code{hyperGoutput} will output tables for all significant
 terms, which may end up being quite a few tables. Usually only a few
 terms are of interest, so there is a \code{subset} argument that can be
 used to select only those terms. This argument follows directly from
 the order of the table output by \code{\link[Category:HyperGResult-accessors]{htmlReport}} or
 \code{\link[GOstats:GOHyperGResult-class]{summary}}. For instance, if the first, third and
 fifth terms in the HTML table output by
 \code{\link[Category:HyperGResult-accessors]{htmlReport}} were of interest, one would use
 subset=c(1,3,5).

 One critical step prior to the hypergeometric test is to subset the
 probesets to unique Entrez Gene IDs. It should be noted however, that
 the functions used by \code{hypergOutput} will output all the probesets
 annotated to a particular term. The \code{output} argument is used to
 control this behavior. If output = "significant" (the default), then only
 those probesets that correspond to the original subsetting will be
 output. If output = "all", then all probesets will be output (grouped
 by Entrez ID), with the
 'significant' probeset first. If output = "split", then all the probesets
 will be output, with all the 'significant' probesets first, followed by
 the other probesets, grouped by Entrez ID.

 Note that the 'significant' probesets come from one of two
 sources. First, one can pass a character vector of probeset IDs
 corresponding to those that were significant in the original analysis
 (recommended). Second, if the \code{geneIds} slot of the \code{GOHyperGParams}
 object containes a named vector of Entrez Gene IDs, then the names from
 that vector will be used. This can be accomplished by using either
 \code{\link[genefilter]{findLargest}} or \code{getUniqueLL}.

 Since the \code{geneIds} are by definition a unique set of Entrez Gene
 IDs, any duplicate probeset IDs will have been removed, so the first
 method is to be preferred for accuracy.
 
}
\value{
  This function returns no value, and is called solely for the side
  effect of outputting HTML and/or text tables.}

\author{ James W. MacDonald <jmacdon@med.umich.edu> }

\seealso{ \code{\link[Category]{hyperGTest}},
  \code{\link[Category:HyperGResult-accessors]{htmlReport}},
  \code{\link[GOstats]{probeSetSummary}}} 
\keyword{manip }

