\name{limma2annaffy}
\alias{limma2annaffy}
\alias{limma2annaffy.na}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to Create HTML Tables from limma Objects }
\description{
  This function is designed to take an \code{ExpressionSet} and an \code{lmFit}, \code{model.matrix},
  and contrast object from limma and convert into HTML tables using
  annaffy. The alternate function \code{limma2annaffy.na} is designed to
  be run without user intervention.
}
\usage{
limma2annaffy(eset, fit, design, contrast, lib, adjust = "fdr",
 anncols = aaf.handler()[c(1:3, 7:8, 10:13)], number = 30, pfilt = NULL,
 fldfilt= NULL,tstat = TRUE, pval = TRUE, FC = TRUE,
 expression = TRUE, html = TRUE, text = FALSE, save = FALSE, addname =
NULL, interactive = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{ An \code{ExpressionSet} containing affymetrix
    expression values. }
  \item{fit}{An \code{lmFit} object. }
  \item{design}{ A \code{model.matrix} object. }
  \item{contrast}{A contrasts matrix from limma. }
  \item{lib}{ An annotation package for the Affy chips used. }
  \item{adjust}{ Multiplicity adjustment. Choices are
    "fdr","holm","hommel","bonferroni", or "none". Partial matching allowed. }
  \item{anncols}{A vector of things to annotate,
    produced by a call to aaf.handler().}
  \item{number}{ Number of genes to output to table. See
    details for more information. }
  \item{pfilt}{ A p-value to filter output. See details
    for more information. }
  \item{fldfilt}{ A fold change to filter output. See
    details for more information. }
  \item{tstat}{Boolean: Output t-statistics in table? Defaults to \code{FALSE}. }
  \item{pval}{Boolean: Output (adjusted) p-values in table?
    Defaults to \code{FALSE}. }
  \item{FC}{Boolean: Output fold changes in table? Defaults to \code{FALSE}. }
  \item{expression}{ Boolean: Output expression values in
    table? Defaults to \code{TRUE}. }
  \item{html}{Boolean: Output data in HTML tables? Defaults to \code{TRUE}. }
  \item{text}{Boolean: Output data in text tables? Defaults to \code{TRUE}. }
  \item{save}{ Boolean: Save tables as R objects for further
    processing? Defaults to \code{FALSE}. }
  
  \item{addname}{ A character vector to add to the end of
  the automatically generated output file names. Useful for multiple
  calls to eliminate over-writing of existing HTML or text tables.}
  \item{interactive}{Boolean: Is this an interactive call, or run as
    part of a script (e.g., in an \code{Sweave} document)? Defaults to \code{TRUE}}
}
\details{
  This function is designed to automatically output HTML or text tables,
  with filenames taken from the column names of the contrast matrix. The
  number of genes output can be controlled several different
  ways. First, if pfilt and fldfilt are both \code{NULL}, the top genes
  will be output based on the \code{number} variable. Otherwise, the
  genes are filtered based on p-value, fold change, or both. If the
  genes are filtered this way, the number of genes to be output will be
  listed and the filter(s) can then be adjusted if necessary.

  This function currently only supports Affymetrix data.
}
\value{
 If \code{save} is \code{TRUE}, a list of tables from
 \code{\link[limma]{topTable}} will be output.
}

\author{ James W. MacDonald <jmacdon@med.umich.edu> }

\seealso{ \code{\link[limma]{topTable}}, \code{\link[annaffy]{aafTableAnn}}}

\keyword{ manip }% at least one, from doc/KEYWORDS
