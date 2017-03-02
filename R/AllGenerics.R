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
##' data match up correctly with the expression data). Current choices for the annoation
##' data are a ChipDb object (e.g., hugene10sttranscriptcluster.db) or an AffyGenePDInfo
##' object (e.g., pd.hugene.1.0.st.v1). In the latter case, we use the parsed Affymetrix
##' annotation csv file to get data. This is only intended for those situations where the
##' ChipDb package is not available.
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
