\name{entrezLinks}
\alias{entrezLinks}
\alias{affyLinks}
\alias{goLinks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Add links to data when using ReportingTools}
\description{
  These functions are intended to add links to the Affymetrix, Entrez
  Gene, and AmiGO databases when creating HTML tables using ReportingTools.
}
\usage{
entrezLinks(df, ...)
affyLinks(df, ...)
goLinks(df, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{A data.frame, usually created using the \code{select}
    function of the AnnotationDbi package. For Entrez ID data, the
    column name must be ENTREZID. For Affy data, the column name must be
    PROBEID, and for GO data the column name must be Term. Any other
    names will fail.}
  \item{...}{Allows one to pass arbitrary arguments to lower level
    functions. Currently unsupported.}
}
\details{
  These functions are not actually intended to be called
  directly. Instead, they are used as targets for the .modifyDF argument
  of the \code{publish} function of ReportingTools. See the example
  below for more detail.
}
\value{
  A data.frame is returned, with links included.
}
\author{James W. MacDonald <jmacdon@u.washington.edu>}
\examples{
\dontrun{
## say we have an ExpressionSet from HuGene 1.0 ST array
## and fit a model using limma
fit <- lmFit(eset, design)
fit2 <- eBayes(fit)
## now annotate the genes using select()
fit2$genes <- select(hugene10sttranscriptcluster.db, featureNames(eset),
c("ENTREZID","SYMBOL"))
## and create an HTML page with links to Affy and Entrez
out <- topTable(fit2, coef=2)
htab <- HTMLReport("The title","a_short_name")
publish(out, htab, .modifyDF = list(affyLinks, entrezLinks))
finish(htab)}
}
\keyword{ manip}% at least one, from doc/KEYWORDS
