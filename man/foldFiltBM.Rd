\name{foldFiltBM}
\alias{foldFiltBM}
\title{Output Fold Change Data using biomaRt}
\usage{
  foldFiltBM(object, fold = 1, groups, comps, compnames,
    species, links = linksBM()[1:3],
    otherann = annBM()[1:3], filterfun = NULL,
    ann.source = "entrezgene", affyid = FALSE, html = TRUE,
    text = TRUE, save = FALSE)
}
\arguments{
  \item{object}{An \code{ExpressionSet} object}

  \item{fold}{The log fold change cutoff to use. Note that
  this is log base two.}

  \item{groups}{A vector of group identifiers. Probably
  easiest to use a numeric vector}

  \item{comps}{A list containing all the comparisons to be
  made. Each list item should be a vector of length two.
  See details for more information.}

  \item{compnames}{A character vector of the names for each
  of the comparisons to be made. This will be the name of
  the resulting HTML or text file.}

  \item{species}{The species name. This must be in a
  particular format for biomaRt. An example for human is
  "hsapiens" or for mouse is "mmusculus".}

  \item{links}{A character vector of things to annotate
  with hyperlinks to online databases. See \code{linksBM}
  for possible values.}

  \item{otherann}{A character vector of things to annotate
  with text only (i.e., no hyperlinks). See \code{annBM}
  for possible values.}

  \item{filterfun}{A filtering function created by
  \code{\link[genefilter]{genefilter}} to filter the data
  using additional criteria. See details for more
  information}

  \item{ann.source}{The annotation source of the IDs that
  will be used to annotate the genes. The default value is
  "entrezgene". See details for other possibilities.}

  \item{affyid}{Boolean. Are the IDs used to annotate these
  data Affymetrix IDs?}

  \item{html}{Boolean. Output HTML tables? Defaults to
  \code{TRUE}}

  \item{text}{Boolean. Output text tables? Defaults to
  \code{TRUE}}

  \item{save}{Boolean. If \code{TRUE}, a list will be
  returned. The first item in the list will be a vector
  showing the number of 'significant' genes for each
  comparison. The second item will be a matrix of -1's, 0's
  and 1's indicating a significant difference, and the
  direction of the difference. The first item is useful for
  creating Sweave - based reports and the second is useful
  for making Vennn diagrams using \code{vennDiagram} from
  the limma package.}
}
\value{
  Returns a list; see above for the elements of the list.
  This function is mainly called for the side effect of
  outputting HTML or text files containing annotated
  'significant' gene lists.
}
\description{
  This function is designed to take an \code{ExpressionSet}
  and some comparisons and output HTML tables. It is very
  similar to \code{foldFilt} except it uses the
  \code{biomaRt} package to annotate genes and the annotate
  package to create the HTML table(s).
}
\details{
  This function is useful for outputting annotated gene
  lists for multiple fold change comparisons. The genes
  will be ordered by the absolute fold change.

  This function currently only supports Affymetrix data. It
  is designed for Affymetrix chips that don't have an
  annotation package, which includes data that have been
  analyzed using the 're-mapped' CDFs supplied to BioC by
  MBNI at University of Michigan.

  The IDs that will be used to annotate the genes depend on
  the source of the data. If, for example, one is using an
  Affymetrix chip that doesn't have a BioC annotation
  package, then the IDs will be Affymetrix IDs. To find out
  the correct name to use for the ann.source argument, one
  can create a connection to a Biomart database using
  \code{\link[biomaRt]{useMart}} and then get a list of
  available Affy arrays using \code{getAffyArrays}.

  If one is using one of the re-mapped CDFs from MBNI at
  University of Michigan, then the IDs to use depend on the
  mapping used to create the CDF. At this time, only three
  types of CDFs can be used; EntrezGene, UniGene, and
  RefSeq. One can determine the correct ann.source argument
  by creating a connection to a Biomart database, and then
  calling \code{linksBM(mart, linksBM())[[3]]}.

  One can also protect against selecting probesets that
  have very small expression values for all samples (which
  likely have a large fold change due to noise, rather than
  signal) by using the filterfun argument. An example would
  be:

  f <- kOverA(1, 6)

  filt <- filterfun(f)

  Then add filterfun = filt as an argument to the call to
  \code{foldFilt}.
}
\author{
  James W. MacDonald <jmacdon@u.washington.edu>
}
\keyword{manip}

