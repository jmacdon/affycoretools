% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/geneSetTables.R
\name{runRomer}
\alias{runRomer}
\alias{runRomer,ExpressionSet-method}
\alias{runRomer,DGEList-method}
\alias{runRomer,EList-method}
\title{A function to run the romer function on a set of contrasts.}
\usage{
runRomer(object, ...)

\S4method{runRomer}{ExpressionSet}(
  object,
  fit,
  setloc,
  annot = NULL,
  design = NULL,
  contrast = NULL,
  wts = NULL,
  save = TRUE,
  baseline.hmap = TRUE,
  affy = TRUE,
  ...
)

\S4method{runRomer}{DGEList}(
  object,
  fit,
  setloc,
  design = NULL,
  contrast = NULL,
  save = TRUE,
  baseline.hmap = TRUE,
  ...
)

\S4method{runRomer}{EList}(
  object,
  fit,
  setloc,
  design = NULL,
  contrast = NULL,
  save = TRUE,
  baseline.hmap = TRUE,
  ...
)
}
\arguments{
\item{object}{An ExpressionSet, DGEList, or EList object}

\item{...}{Used to pass arguments to lower-level functions. See
\code{outputRomer} \code{geneSetPage}, \code{dataAndHeatmapPage} and
\code{gsHeatmap} for available arguments.}

\item{fit}{A fitted model from either limma (e.g., MArrayLM) or edgeR (e.g., DGEGLM)}

\item{setloc}{A character vector giving the path for gene set RData files
(see description for more information), or a named list (or list of lists),
where the top-level names consist of gene set grouping names (like KeGG or GO),
the next level names consist of gene set names (like NAKAMURA_CANCER_MICROENVIRONMENT_UP),
and the list items themselves are gene symbols, matching the expected capitalization for
the species being used (e.g., for human, they are ALL CAPS. For most other species only the
First Letter Is Capitalized).}

\item{annot}{Character. The name of the array annotation package. If NULL, the
annotation data will be extracted from the fData slot (for ExpressionSets) or the
genes list (for DGEList or EList objects).}

\item{design}{A design matrix describing the model fit to the data. Ideally this should be
a cell-means model (e.g., no intercept term), as the design and contrast matrices are used
to infer which data to include in the output heatmaps. There is no guarantee that this will work
correctly with a treatment-contrasts parameterization (e.g., a model with an intercept).}

\item{contrast}{A contrast matrix describing the contrasts that were
computed from the data. This contrast should have colnames, which will be
used to create parts of the resulting directory structure.}

\item{wts}{Optional weights vector - if array weights were used to fit the
model, they should be supplied here as well.}

\item{save}{Boolean. If true, after running the \code{\link{romer}} step,
the results will be saved in a file 'romer.Rdata', which can be used as
input for \code{outputRomer} to create HTML tables. Since
\code{\link{romer}} can take a long time to run, it is advantageous to keep
the default.}

\item{baseline.hmap}{Boolean. If \code{TRUE}, then the resulting heatmaps
will be centered by subtracting the mean of the baseline sample. As an
example, in a contrast of treatment A - treatment B, the mean of the
treatment B samples will be subtracted. The heatmap colors then represent
the fold change between the A and B samples.}

\item{affy}{Boolean; are these Affymetrix arrays? If \code{TRUE}, the output tables
will contain links to the netaffx site.}
}
\value{
If save is TRUE, return a list that can be re-processed using \code{outputRomer}.
this is useful in cases where you might need to re-run multiple times.

Nothing is returned. This function is called only for the
side-effects of creating output HTML files in the working and
sub-directories.
}
\description{
This function automates both running \code{\link{romer}} on a set of
contrasts as well as the creation of output HTML tables that can be used to
explore the results.  The basic idea here is that one might have used limma
to fit a model and compute some contrasts, and then want to do a GSEA using
\code{\link{romer}}.
}
\details{
The \code{\link{romer}} expects as input a list or lists of gene symbols
that represent individual gene sets. One example is the various gene sets
from the Broad Institute that are available at
http://bioinf.wehi.edu.au/software/MSigDB/, which are distributed as RData
files. The default assumption for this function is that the end user will
have downloaded these files, and the setloc argument simply tells
\code{runRomer} where to find them.

Alternatively, user-based gene sets could be created (these should consist
of lists of character vectors of gene symbols - see one of the Broad gene
sets for an example).

This function will run \code{\link{romer}} using all the gene sets in the
referenced directory, on all the contrasts supplied, and then output the
results in a (default) 'genesets' subdirectory. There will be an HTML file
in the working directory with a (default) filename of 'indexRomer.html' that
will point to individual HTML files in the genesets subdirectory, which will
point to individual files in subdirectories within the genesets subdirectory
(named after the colnames of the contrast matrix).
}
\section{Methods (by class)}{
\itemize{
\item \code{ExpressionSet}: Perform gene set analysis using microarray data.

\item \code{DGEList}: Perform gene set analysis using RNA-Seq data processed using edgeR.

\item \code{EList}: Perform gene set analysis using RNA-Seq data processed using voom.
}}

\author{
James W. MacDonald <jmacdon@u.washington.edu>

James W. MacDonald <jmacdon@u.washington.edu>
}
