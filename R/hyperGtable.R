##########################################
##
##  Copyright 2005 James W. MacDonald 
##
##  A function to output a table from running
##  GOHyperG() on a set of Affymetrix probe IDs
##
##
##  Modified 6-9-05 to output the table
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
