\name{probes2tableBM}
\alias{probes2tableBM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Convert Affy Probe ids to Annotated HTML Table using biomaRt }
\description{
 A function to convert a vector of Affy ids to an annotated HTML or text
 table. This function is very similar to \code{probes2table}, except it
   uses the \code{biomaRt} package to annotate genes, and the
   \code{annotate} package to create the HTML table.
}
\usage{
probes2tableBM(eset, probids, species, filename, otherdata = NULL,
 links = linksBM()[1:3], otherann = annBM()[1:3], ann.source = "entrezgene",
 express = TRUE, html = TRUE, text = TRUE, affyid = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{ An \code{ExpressionSet} containing Affy expression
  values.}
  \item{probids}{ A vector of probe ids. }
  \item{species}{The species name. This must be in a particular format
    for biomaRt. An example for human is "hsapiens" or for mouse is
    "mmusculus".}
  \item{filename}{File name of the resulting HTML table.}
  \item{otherdata}{A *named* list of additional information to include
    in the resulting table. Examples would be t-statistics, p-values,
    fold change, etc. Each list item should be a vector the same length
    as the probids vector. The name associated with each list item will
    be used as the column name in the resulting table.}
  \item{links}{A character vector of things to annotate with hyperlinks
    to online databases. See \code{linksBM} for possible values.}
  \item{otherann}{A character vector of things to annotate with text
    only (i.e., no hyperlinks). See \code{annBM} for possible values.}
  \item{ann.source}{The annotation source of the IDs that will be used
    to annotate the genes. The default value is "entrezgene". See
    details for other possibilities.}
  \item{express}{ Output expression values in table?
    Defaults to \code{TRUE}. }
  \item{html}{Boolean. Output HTML table? Defaults to \code{TRUE}}
  \item{text}{Boolean. Output text table? Defaults to \code{TRUE}}
  \item{affyid}{Boolean. Are the IDs used to annotate these data
    Affymetrix IDs?}
  }
\details{
  This function is designed to output HTML tables based on a set of
  IDs. This function currently only supports Affymetrix data. It
  is designed for Affymetrix chips that don't have an
  annotation package, which includes data that have been analyzed using
  the 're-mapped' CDFs supplied to BioC by MBNI at University of
  Michigan.

  The IDs that will be used to annotate the genes depend on the source
  of the data. If, for example, one is using an Affymetrix chip that
  doesn't have a BioC annotation package, then the IDs will be
  Affymetrix IDs. To find out the correct name to use for the ann.source
  argument, one can create a connection to a Biomart database using
  \code{\link[biomaRt]{useMart}} and then get a list of available Affy
  arrays using \code{\link[biomaRt]{listFilters}}.

  If one is using one of the re-mapped CDFs from MBNI at University of
  Michigan, then the IDs to use depend on the mapping used to create the
  CDF. At this time, only three types of CDFs can be used; EntrezGene,
  UniGene, and RefSeq. One can determine the correct ann.source argument
  by creating a connection to a Biomart database, and then calling
  \code{linksBM(mart, linksBM())[[3]]}.
}

\value{
 This function is only used for the side effect of outputting an HTML table.
}
\author{ James W. MacDonald <jmacdon@u.washington.edu> }

\seealso{ \code{\link[limma:toptable]{topTable}}}

\keyword{ manip }% at least one, from doc/KEYWORDS
