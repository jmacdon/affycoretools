\name{vennCounts2}
\alias{vennCounts2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute Counts for Venn Diagram}
\description{
  This function is designed to compute counts for a Venn diagram. It is
  slightly different from \code{\link[limma]{vennCounts}} in the
  additional ability to compute counts for genes that are differentially
  expressed in the same direction.
}
\usage{
vennCounts2(x, method = "same", fit = NULL, foldFilt = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A \code{\link[limma]{TestResults}} object, produced by a
    call to \code{\link[limma]{decideTests}} or \code{foldFilt}.}
  \item{method}{ One of "same", "both", "up", "down". See details for
    more information.}
  \item{fit}{An \code{\link[limma]{MArrayLM}} object, produced by a call
    to \code{\link[limma]{lmFit}} and \code{\link[limma]{eBayes}}. Only
    necessary if 'foldFilt' = \code{TRUE}.}
  \item{foldFilt}{A fold change to filter samples. This is primarily
    here for consistency with the corresponding argument in
    \code{vennSelect}.}
}
\details{
  The function \code{\link[limma]{vennCounts}} will return identical
  results except for the "same" method. This will only select those
  genes that both pass the criteria of \code{\link[limma]{decideTests}}
  as well as being differentially expressed in the same direction. Note
  that this is different from the "both" method, which simply requires
  that a given gene be differentially expressed in e.g., two different
  comparisons without any requirement that the direction be the same.
}
\value{
  A \code{\link[limma]{VennCounts}} object.
}
\author{ James W. MacDonald <jmacdon@med.umich.edu>}
\examples{
tstat <- matrix(rt(300,df=10),100,3)
tstat[1:33,] <- tstat[1:33,]+2
clas <- classifyTestsF(tstat,df=10,p.value=0.05)
a <- vennCounts2(clas)
print(a)
vennDiagram(a)
}
\keyword{ manip }% at least one, from doc/KEYWORDS
\keyword{ hplot }% __ONLY ONE__ keyword per line
