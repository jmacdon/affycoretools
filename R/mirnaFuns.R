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
        ids <- sapply(strsplit(ids, "_"), function(x) if(length(x) > 2) x[2] else x[1])
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
        return(mat)
}

  
