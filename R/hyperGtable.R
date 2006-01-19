##########################################
##
##  Copyright 2005 James W. MacDonald 
##
##  A function to output a table from running
##  GOHyperG() on a set of Affymetrix probe IDs
##
##
##  Modified 6-9-05 to output the table
##  Added hyperG2Affy 10-10-2005
##   - This function takes GOIDs form hyperGtable and outputs
##     a list giving the Affy IDs associated with each GOID
##
##  Added hyperG2annaffy 12-12-05
##    -This function takes the output of hyperGtable and outputs
##     HTML tables from the Affy IDs associated with each GOID
###########################################

hyperGtable <- function(probids, lib, type="MF", pvalue=0.05,
                        min.count=10, save = FALSE, output = TRUE,
                        filename = NULL){
  require("GOstats", quietly = TRUE) || stop("The GOstats package is required")
  require(lib, quietly = TRUE, character.only = TRUE) || stop(paste("The ", lib, " package is required"))
  lls <- getLL(probids, lib)
  lls <- unique(lls)
  lls <- lls[!is.na(lls)]
  tmp <- GOHyperG(lls, lib, type)
  index <- tmp$pvalues < pvalue & tmp$goCounts > min.count
  wh <- mget(names(tmp$pvalues[index]), GOTERM)
  tmp.terms <- sapply(wh, Term)
  out <- data.frame(names(tmp$pvalues[index]), tmp.terms, round(tmp$pvalues[index], 3), 
                       tmp$intCounts[index], tmp$goCounts[index])
  names(out) <- c("GO Term","Description","p-values",
                     "Number significant","Number on chip")
  if(output){
    if(!is.null(filename)){
      write.table(out, paste(filename, ".txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
    }else{
      stop("You need to give a filename in order to output this table!\n")
    }
  }
  if(save)
    out
}

hyperG2Affy <- function(probids, lib, type="MF", pvalue=0.05,
                        min.count=10){
  require("GOstats", quietly = TRUE) || stop("The GOstats package is required")
  require(lib, quietly = TRUE, character.only = TRUE) || stop(paste("The ", lib, " package is required"))
  lls <- getLL(probids, lib)
  lls <- unique(lls)
  lls <- lls[!is.na(lls)]
  tmp <- GOHyperG(lls, lib, type)
  index <- tmp$pvalues < pvalue & tmp$goCounts > min.count
  wh <- mget(names(tmp$pvalues[index]), GOTERM)
  tmp.terms <- sapply(wh, Term)
  index2 <- match(names(tmp.terms), names(tmp$go2Affy))
  out <- lapply(tmp$go2Affy[index2], function(x) x[x %in% probids])
  names(out) <- names(tmp.terms)
  out
}
 
  
hyperG2annaffy <- function(probids, lib, eset, fit = NULL, subset = NULL, comp = 1,
                           type="MF", pvalue = 0.05, min.count = 10){
  tab <- hyperGtable(probids = probids, lib = lib, type = type, pvalue = pvalue,
                     min.count = min.count, save = TRUE, output = FALSE)
  if(!is.null(subset)) tab <- tab[subset,]
  prbs <- vector("list", dim(tab)[1])
  for(i in seq(along = prbs)) 
    prbs[[i]] <- probids[probids %in% get(as.character(tab[i,1]), 
                                          get(paste(lib, "GO2ALLPROBES", sep = "")))]
  for(i in seq(along = prbs)){
    class(prbs[[i]]) <- "character"
    index <- geneNames(eset) %in% prbs[[i]]
    if(!is.null(fit)){
      ord <- order(abs(fit$t)[index], decreasing = TRUE)
      otherdata <- list("t-statistic" = round(fit$t[index][ord], 2),
                        "p-value" = round(p.adjust(fit$p.value, "fdr")[index][ord], 3),
                        "Fold Change" = round(fit$coef[index][ord],2))
      probes2table(eset, prbs[[i]][ord], annotation(eset), otherdata = otherdata,
                   filename = paste(as.character(tab[i,2]), "genes", sep=" "))
    }else{
      probes2table(eset, prbs[[i]], annotation(eset),
                   filename = paste(as.character(tab[i,2]), "genes", sep = " "))
    }
  }
}
  
