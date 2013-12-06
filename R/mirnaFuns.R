########################################
##
## Copyright 2012 James W. MacDonald
##
##
## mirnaFuns.R - Functions for miRNA data, primarily in the context of
##               existing mRNA data.
##
##
## 2012-04-09 First pass at functions
##
###########################################

convertIDs <- function(ids, from = TRUE){
    if(from) {
        ids <- gsub("_st|_x_st|_s_st|_at", "", ids)
        ids <- gsub("-star", "*", ids)
    }else{
        ids <- gsub("\\*", "-star", ids)
        ids <- paste(ids, "_st", sep = "")
    }
    ids
}

allToMat <- function(object){
    if(is.matrix(object)){
        return(object)
    } else {
        if(is(object, "ExpressionSet")){
            return(exprs(object))
        } else {
            if(is.data.frame(object)){
                return(as.matrix(object))
            } else {
                stop(paste(ls()[1], "should be a matrix, ExpressionSet or data.frame"),
                     call. = FALSE)
            }
        }
    }
}



#' A function to map miRNA to mRNA.
#' 
#' This function is intended use when there are miRNA and mRNA data for the
#' same subjects, and the goal is to detect mRNAs that appear to be targeted by
#' the miRNA.
#' 
#' This function is intended to take a vector of miRNA IDs that are
#' significantly differentially expressed in a given experiment and then map
#' those IDs to putative mRNA transcripts that the miRNAs are supposed to
#' target. The mRNA transcript IDs are then mapped to chip-specific probeset
#' IDs, which are then subsetted to only include those probesets that were also
#' significantly differentially expressed.
#' 
#' The output from this function is intended as input for
#' \code{\link{makeHmap}}.
#' 
#' @param miRNAids A character vector of miRNA IDs. Currently only supports
#' Affymetrix platform.
#' @param miRNAannot Character. The filename (including path if not in working
#' directory) for the file containing miRNA to mRNA mappings.
#' @param mRNAids A character vector of mRNA IDs. Currently only supports
#' Affymetrix platform.
#' @param orgPkg Character. The Bioconductor organism package (e.g.,
#' org.Hs.eg.db) to be used for mapping.
#' @param chipPkg Character. The Bioconductor chip-specific package (e.g.,
#' hgu133plus2.db) to be used for mapping.
#' @param sanger Boolean. Is the miRNAannot file a Sanger miRBase targets file?
#' These can be downloaded from
#' http://www.ebi.ac.uk/enright-srv/microcosm/cgi-bin/targets/v5/download.pl
#' @param miRNAcol Numeric. If using a Sanger miRBase targets file, leave
#' \code{NULL}. Otherwise, use this to indicate which column of the miRNAannot
#' file contains miRNA IDs.
#' @param mRNAcol Numeric. If using Sanger miRBase targets file, leave
#' \code{NULL}. Otherwise, use this to indicate which column of the miRNAannot
#' file contains mRNA IDs.
#' @param transType Character. Designates the type of transcript ID for mRNA
#' supplied by the miRNAannot file. If using the Sanger miRBase files, this is
#' ensembl. Other choices include refseq and accnum.
#' @return A list with names that correspond to each significant miRNA, and the
#' mRNA probeset IDs that are targeted by that miRNA.
#' @author James W. MacDonald
#' @seealso makeHmap
#' @keywords manip
#' @export mirna2mrna
mirna2mrna <- function(miRNAids, miRNAannot, mRNAids, orgPkg, chipPkg, sanger = TRUE,
                       miRNAcol = NULL, mRNAcol = NULL, transType = "ensembl"){
    transType <- match.arg(transType, c("ensembl","refseq","accnum"))
    if(!sanger && (is.null(miRNAcol) || is.null(mRNAcol)))
        stop(paste("If the annotation data are not Sanger microcosm target data\n",
                   "then you must supply the correct column number for miRNA and mRNA data.\n\n",
                   sep = ""), call. = FALSE)
    if(sanger){
        if(is.null(miRNAcol)) miRNAcol <- 2
        if(is.null(mRNAcol)) mRNAcol <- 12
    }
    tab <- switch(transType,
                  ensembl = "ensembl_trans",
                  "refseq")

    miRNAids <- convertIDs(miRNAids)
    annot <- read.table(miRNAannot, sep = "\t", stringsAsFactors = FALSE)
    
    annotlst <- tapply(1:nrow(annot), annot[,miRNAcol], function(x) annot[x,mRNAcol])
    ## for each miRNA get corresponding transcript IDs
    prblst <- annotlst[names(annotlst) %in% miRNAids]

    ## convert to probeIDs
    require(chipPkg, quietly = TRUE, character.only = TRUE)
    require(orgPkg, quietly = TRUE, character.only = TRUE)
    orgDBloc <- system.file("extdata", sub("db","sqlite", orgPkg),
                            package = orgPkg)
    attachSQL <- paste("attach '", orgDBloc, "' as orgDB;", sep = "")
    con <- dbconn(get(paste(sub("\\.db", "", chipPkg), "ENTREZID", sep = "")))
    dbGetQuery(con, attachSQL)

    makeSql <- function(ids, tab){
        paste("select distinct probe_id from probes ",
              "inner join orgDB.genes as gi on probes.gene_id=gi.gene_id ",
              "inner join orgDB.", tab, " as et on gi._id=et._id ",
              "where et.trans_id in ('", paste(ids, collapse = "','"), "');", sep = "")
    }
    
    mrnaprbs <- lapply(prblst, function(x) dbGetQuery(con, makeSql(x, tab))[,1])
    intprbs <- lapply(mrnaprbs, function(x) mRNAids[mRNAids %in% x])
    intprbs <- intprbs[sapply(intprbs, length) > 0]
    dbGetQuery(con, "detach orgDB;")
    intprbs
}



#' A function to create a heatmap-like object or matrix of correlations between
#' miRNA and mRNA data.
#' 
#' This function is intended for use when both miRNA and mRNA data are
#' available for the same samples. In this situation it may be advantageous to
#' compute correlations between the two RNA types, in order to detect mRNA
#' transcripts that are targeted by miRNA.
#' 
#' As noted above, this function is intended to generate output from
#' simultaneous analyses of miRNA and mRNA data for the same samples, the goal
#' being either a heatmap like plot of correlations, or the data (or both).
#' 
#' If creating a plot, note that if the number of significant mRNA probes is
#' large, the resulting heatmap will have many rows and will not plot correctly
#' on the usual graphics device within R. In order to visualize, it is almost
#' always better to output as a pdf. In addition, the dimensions of this pdf
#' will have to be adjusted so the row names for the heatmap will be legible.
#' As an example, a heatmap with 10 miRNA transcripts and 100 mRNA transcripts
#' will likely need a pdf with a width argument of 6 and a height argument of
#' 25 or 30. It may require some experimentation to get the correct arguments
#' to the \code{pdf} function.
#' 
#' Also please note that this function by necessity outputs rectangular data.
#' However, there will be many instances in which a given miRNA isn't thought
#' to target a particular mRNA. Whenever this occurs, the heatmap will have a
#' white cell, and the output data for that combination will be NA.
#' 
#' @param mRNAdat An \code{ExpressionSet}, \code{data.frame} or \code{matrix}
#' of mRNA expression values. The row.names for these data should correspond to
#' the manufacturer's probe ID. Currently, the only manufacturer supported is
#' Affymetrix.
#' @param miRNAdat An \code{ExpressionSet}, \code{data.frame} or \code{matrix}
#' of mRNA expression values. The row.names for these data should correspond to
#' the manufacturer's probe ID. Currently, the only manufacturer supported is
#' Affymetrix.
#' @param mRNAlst A \code{list} of mRNA probe IDs where the names of each list
#' item are mirBase miRNA IDs. Usually this will be the output from
#' \code{\link{mirna2mrna}}.
#' @param mRNAvec A numeric vector used to subset or reorder the mRNA data, by
#' column. If \code{NULL}, this will simply be 1:ncol(mRNAdat).
#' @param miRNAvec A numeric vector used to subset or reorder the miRNA data,
#' by column. If \code{NULL}, this will simply be 1:ncol(miRNAdat).
#' @param chipPkg Character. The name of the chip-specific annotation package
#' (e.g., "hgu133plus2.db").
#' @param header Character. The plot title if a heatmap is output.
#' @param plot Boolean. Should a heatmap be generated?
#' @param out Boolean. Should the matrix of correlation coefficients be output?
#' @return This function will output a numeric matrix if the 'out' argument is
#' \code{TRUE}.
#' @author James W. MacDonald
#' @seealso mirna2mrna
#' @keywords manip
#' @export makeHmap
makeHmap <- function(mRNAdat, miRNAdat, mRNAlst, mRNAvec = NULL, miRNAvec = NULL, chipPkg,
                     header, plot = TRUE, out = TRUE){
    mRNAdat <- allToMat(mRNAdat)
    miRNAdat <- allToMat(miRNAdat)
    if(is.null(mRNAvec)) mRNAvec <- seq_len(ncol(mRNAdat))
    if(is.null(miRNAvec)) miRNAvec <- seq_len(ncol(miRNAdat))
    rn <- unique(do.call("c", mRNAlst))
    cn <- gsub("\\.probe_id", "", names(mRNAlst))
    if(!all(rn %in% row.names(mRNAdat)))
        warning(paste("Not all mRNA probes being tested are found in\n",
                      "the row.names of the mRNA data.\n\n", sep = ""), call. = FALSE)
    
    mat <- matrix(NA, nrow = length(rn), ncol = length(cn))
    row.names(mat) <- rn
    colnames(mat) <- cn
    cn2 <- convertIDs(cn, FALSE)
    ## have to account for probesets that start with v11_
    bads <- cn2[!cn2 %in% row.names(miRNAdat)]
    if(length(bads) > 0){
        ind <- which(cn2 %in% bads)
        for(i in seq(along = bads))
            cn2[ind[i]] <- grep(cn2[ind[i]], row.names(miRNAdat), value = TRUE)
    }
    if(!all(cn2 %in% row.names(miRNAdat)))
        warning(paste("Not all miRNA probes being tested are found in\n",
                      "the row.names of the miRNA data.\n\n", sep = ""), call. = FALSE)
    
    getCor <- function(miRNA, mRNAlst.itm){
        sapply(seq(along = mRNAlst.itm), function(x) cor(miRNAdat[miRNA,miRNAvec],
                   mRNAdat[mRNAlst.itm[x],mRNAvec]))
    }
    for(i in seq(along = cn2))
        mat[mRNAlst[[i]], cn[i]] <- getCor(cn2[i], mRNAlst[[i]])
    
    rn2 <- sapply(AnnotationDbi::mget(rn, get(paste(sub("\\.db", "", chipPkg), "SYMBOL", sep = "")),
                       ifnotfound=NA), function(x) x[1])
    rn2 <- ifelse(is.na(rn2), rn, rn2)
    row.names(mat) <- rn2
    ## alphabetize matrix
    ord <- order(row.names(mat))
    mat <- mat[ord,,drop = FALSE]
    if(plot){
        heatmap.2(mat, Rowv = FALSE, Colv = FALSE, dendrogram = "none",
                  trace = "none", density.info = "none", lhei = c(0.05, 0.95),
                  cexCol = 0.6, colsep = 1:ncol(mat), rowsep = 1:nrow(mat),
                  sepcolor = "lightgrey", main = header, margins = c(5,8))
    }
    if(out)
        mat <- data.frame(mRNA.ID = rn[ord], mat, row.names = NULL)
        return(mat)
}

  
