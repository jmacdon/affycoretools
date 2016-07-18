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


##' @param columns For ChipDb method; what annotation data to add. Use the \code{\link[AnnotationDbi]{columns}}
##' function to see what choices you have. By default we get the ENTREZID, SYMBOL and GENENAME.
##' @param multivals For ChipDb method; this is passed to \code{\link[AnnotationDbi]{mapIds}} to control how 1:many
##' mappings are handled. The default is 'first', which takes just the first result. Other valid
##' values are 'list' and 'CharacterList', which return all mapped results.
##' @return An ExpressionSet that has annotation data added to the featureData slot.
##' @author Jim MacDonald
##' @describeIn annotateEset Annotate an ExpressionSet using a ChipDb package for annotation data.
##' @export
setMethod("annotateEset", c("ExpressionSet","ChipDb"),
          function(object, x, columns = c("PROBEID","ENTREZID","SYMBOL","GENENAME"), multivals = "first"){
    ## Doing mapIds(chipDb, featureNames(object), "PROBEID","PROBEID") fails for many packages, and wastes compute cycles
    ## so we just add that by hand.
    addcol <- FALSE
    if(any(columns == "PROBEID")) {
        addcol <- TRUE
        columns <- columns[-grep("PROBEID", columns)]
        collst <- list(PROBEID = featureNames(object))
    }
    multivals <- switch(multivals,
                        first = "first",
                        list = "CharacterList",
                        CharacterList = "CharacterList",
                        stop("The multivals argument should be 'first', 'list' or ;CharacterList'", call. = FALSE))
    annolst <- lapply(columns, function(y) mapIds(x, featureNames(object), y, "PROBEID", multiVals = multivals))
    if(addcol) annolst <- c(collst, annolst)
    anno <- switch(multivals,
                   first = as.data.frame(annolst),
                   DataFrame(annolst))
    names(anno) <- if(addcol) c("PROBEID", columns) else columns
    if(!isTRUE(all.equal(row.names(anno), featureNames(object))))
        stop(paste("There appears to be a mismatch between the ExpressionSet and",
                   "the annotation data.\nPlease ensure that the summarization level",
                   "for the ExpressionSet and the annotation package are the same.\n"),
             call. = FALSE)
    andf <- dfToAnnoDF(anno)
    featureData(object) <- andf
    stopifnot(validObject(object))
    object
})


##' @param type For pdInfoPackages; either 'core' or 'probeset', corresponding to the 'target' argument
##' used in the call to \code{\link[oligo]{rma}}.
##' @describeIn annotateEset Annotate an ExpressionSet using an AffyGenePDInfo package.
##' @export
setMethod("annotateEset", c("ExpressionSet","AffyGenePDInfo"),
          function(object, x, type = "core", ...){
    .dataFromNetaffx(object, x, type, ...)
})


##' @describeIn annotateEset Annotate an ExpressionSet using an AffyHTAPDInfo package.
##' @export
setMethod("annotateEset", c("ExpressionSet","AffyHTAPDInfo"),
          function(object, x, type = "core", ...){
    .dataFromNetaffx(object, x, type, ...)
})


##' @describeIn annotateEset Annotate an ExpressionSet using an AffyExonPDInfo package.
##' @export
setMethod("annotateEset", c("ExpressionSet","AffyExonPDInfo"),
          function(object, x, type = "core", ...){
    .dataFromNetaffx(object, x, type,...)
})

.dataFromNetaffx <- function(object, x, type = "core", ...){
    typeToGet <- switch(type,
                        core = "netaffxTranscript.rda",
                        probeset = "netaffxProbeset.rda",
                        stop("Type must be either 'core' or 'probeset'", call. = FALSE))
    load(system.file(paste0("/extdata/", typeToGet), package = annotation(x)))
    annot <- pData(get(sub("\\.rda", "", typeToGet)))
    anno <- as.data.frame(do.call(rbind, lapply(strsplit(annot$geneassignment, " // "), "[", 1:3)))
    names(anno) <- c("ID", "SYMBOL","GENENAME")
    row.names(anno) <- switch(type,
                              probeset = as.character(annot$probesetid),
                              core = as.character(annot$transcriptclusterid))
    anno$PROBEID <- row.names(anno)
    anno <- anno[,c(4,1:3)]
    anno <- anno[featureNames(object),]
    if(!isTRUE(all.equal(row.names(anno), featureNames(object))))
        stop(paste("There appears to be a mismatch between the ExpressionSet and",
                   "the annotation data.\nPlease ensure that the summarization level",
                   "for the ExpressionSet and the 'type' argument are the same.\n",
                   "See ?annotateEset for more information on the type argument.\n"),
             call. = FALSE)
    andf <- dfToAnnoDF(anno)
    featureData(object) <- andf
    stopifnot(validObject(object))
    object
}

##' @describeIn annotateEset Method to capture character input.
##' @export
setMethod("annotateEset", c("ExpressionSet","character"),
          function(object, x, ...){
    do.call(require, list(x, character.only = TRUE, quietly = TRUE))
    x <- get(x)
    annotateEset(object, x, ...)
})

##' @describeIn annotateEset Annotate an ExpressionSet using a user-supplied data.frame.
##' @export
##' @param probecol Column of the data.frame that contains the probeset IDs. Can be either
##' numeric (the column number) or character (the column header).
##' @param annocols Column(x) of the data.frame to use for annotating. Can be a vector of numbers
##' (which column numbers to use) or a character vector (vector of column names).

setMethod("annotateEset", c("ExpressionSet","data.frame"),
          function(object, x, probecol = NULL, annocols = NULL, ...){
    if(is.null(probecol) || is.null(annocols))
        stop(paste("You must specify the column containing the probeset IDs (probecol argument)",
                   "and the column(s) containing the annotation data you wish to use (annocols argument).\n"),
             call. = FALSE)
    rn <- as.character(x[,probecol])
    anno <- x[,annocols]
    row.names(anno) <- rn
    anno <- anno[featureNames(object),]
    if(!isTRUE(all.equal(row.names(anno), featureNames(object))))
        stop(paste("There appears to be a mismatch between the ExpressionSet and the",
                   "data.frame that contains the annotation data. Please make sure that",
                   "the column containing the probeset IDs matches the featureNames from",
                   "your ExpressionSet. You may need to subset one or the other to make",
                   "them conform to each other.\n"), call. = FALSE)
    andf <- dfToAnnoDF(anno)
    featureData(object) <- andf
    stopifnot(validObject(object))
    object
})

## functions to convert either data.frame or DataFrames into AnnotatedDataFrames
setGeneric("dfToAnnoDF", function(x) standardGeneric("dfToAnnoDF"))

setMethod("dfToAnnoDF", "data.frame",
          function(x) {
    adf <- AnnotatedDataFrame(data = x)
    return(adf)
})

setMethod("dfToAnnoDF", "DataFrame",
          function(x) {
    adf <- AnnotatedDataFrame(data = as(x, "data.frame"))
    return(adf)
})

## If we use 'list' or 'CnaracterList' to annotate, we may get multiple values returned for a given ID
## in which case, we need to collapse to a single vector (separated by <BR>) so it will all fit into a single
## cell in the resulting table.
listToVector <- function(object, ...){
    if(any(apply(object, 2, class) == "list")){
        ind <- apply(object, 2, class) == "list"
        for(i in seq(along = ind)) object[,i] <- sapply(object[,i], function(x) paste(x, collapse = "<BR>"))
    }
    return(object)
}


                   
