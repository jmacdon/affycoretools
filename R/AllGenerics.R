##' @rdname makeVenn
##' @export
setGeneric("makeVenn", function(object, ...) standardGeneric("makeVenn"))

##' Method to annotate ExpressionSets automatically
##'
##' This function fills the featureData slot of the ExpressionSet automatically, which is
##' then available to downstream methods to provide annotated output. Annotating
##' results is tedious, and can be surprisingly difficult to get right. By annotating
##' the data automatically, we remove the tedium and add an extra layer of security
##' since the resulting ExpressionSet will be tested for validity automatically (e.g., annotation
##' data match up correctly with the expression data). Current choices for the annotation
##' data are a ChipDb object (e.g., hugene10sttranscriptcluster.db) or an AffyGenePDInfo
##' object (e.g., pd.hugene.1.0.st.v1). In the latter case, we use the parsed Affymetrix
##' annotation csv file to get data. This is only intended for those situations where the
##' ChipDb package is not available, and in particular is only available for those packages
##' that contain the parsed annotation csv files (generally, Gene ST arrays, Exon ST arrays and
##' Clariom/HTA/MTA/RTA arrays).
##' @param object An ExpressionSet to which we want to add annotation.
##' @param x Either a ChipDb package (e.g., hugene10sttranscriptcluster.db),
##' or a pdInfoPackage object (e.g., pd.hugene.1.0.st.v1).
##' @param ... Allow users to pass in arbitrary arguments. Particularly useful for
##' passing in columns, multivals, and type arguments for methods.
##' @examples
##' \dontrun{
##' dat <- read.celfiles(filenames = list.celfiles())
##' eset <- rma(dat)
##' ## annotate using ChipDb
##' eset <- annotateEset(eset, hgu10sttranscriptcluster.db)
##' ## or AffyGenePDInfo
##' eset <- annotateEset(eset, pd.hugene.1.0.st.v1)
##' }
##' @export

setGeneric("annotateEset", function(object, x, ...) standardGeneric("annotateEset"))

## functions to convert either data.frame or DataFrames into AnnotatedDataFrames
setGeneric("dfToAnnoDF", function(x) standardGeneric("dfToAnnoDF"))

setGeneric(".getMainProbes", function(object, x, ...) standardGeneric(".getMainProbes"))


##' A function to run the romer function on a set of contrasts.
##' 
##' This function automates both running \code{\link{romer}} on a set of
##' contrasts as well as the creation of output HTML tables that can be used to
##' explore the results.  The basic idea here is that one might have used limma
##' to fit a model and compute some contrasts, and then want to do a GSEA using
##' \code{\link{romer}}.
##' 
##' The \code{\link{romer}} expects as input a list or lists of gene symbols
##' that represent individual gene sets. One example is the various gene sets
##' from the Broad Institute that are available at
##' http://bioinf.wehi.edu.au/software/MSigDB/, which are distributed as RData
##' files. The default assumption for this function is that the end user will
##' have downloaded these files, and the setloc argument simply tells
##' \code{runRomer} where to find them.
##' 
##' Alternatively, user-based gene sets could be created (these should consist
##' of lists of character vectors of gene symbols - see one of the Broad gene
##' sets for an example).
##' 
##' This function will run \code{\link{romer}} using all the gene sets in the
##' referenced directory, on all the contrasts supplied, and then output the
##' results in a (default) 'genesets' subdirectory. There will be an HTML file
##' in the working directory with a (default) filename of 'indexRomer.html' that
##' will point to individual HTML files in the genesets subdirectory, which will
##' point to individual files in subdirectories within the genesets subdirectory
##' (named after the colnames of the contrast matrix).
##' @param object An ExpressionSet, DGEList, or EList object
##' @return If save is TRUE, return a list that can be re-processed using \code{outputRomer}.
##' this is useful in cases where you might need to re-run multiple times.
##' @author James W. MacDonald <jmacdon@u.washington.edu>
##' @export
setGeneric("runRomer", function(object, ...) standardGeneric("runRomer"))


##' A function to create HTML output from the results of running romer on a set
##' of contrasts.
##' 
##' This function is actually intended to be a sub-function of \code{runRomer},
##' but can hypothetically run by itself if the \code{\link{romer}} step has
##' already been done.
##' 
##' This function is intended to be an internal function for \code{runRomer}.
##' However, it is possible that \code{runRomer} errored out after saving the
##' results from running \code{\link{romer}} on a set of contrasts, and all that
##' remains is to create the output HTML.
##' 
##' Please note that the first two arguments to this function have certain
##' expectations. The rsltlst should be the output from running
##' \code{\link{romer}}. If using the saved output from \code{runRomer}, one
##' should first \code{load} the 'romer.Rdata' file, which will introduce a list
##' object with the name 'romerlst' into the working directory, so the first
##' argument should be rsltlst = romerlst.
##' 
##' Second, see the code for runRomer, specifically the line that creates the
##' 'sets' object, which will show how to create the correct genesetlst object.
##' 
##' @param object An ExpressionSet, DGEList or EList object containing normalized, summarized
##' gene expression data.
##' @param fit An MArrayLM or DGEGLM object, containing the fitted data.
##' @author James W. MacDonald <jmacdon@@u.washington.edu>
##' @export
setGeneric("outputRomer", function(object, fit, ...) standardGeneric("outputRomer"))

