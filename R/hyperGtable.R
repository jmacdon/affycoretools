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
##
##  Moved hyperGtable and hyperG2Affy to GOstats 2-23-06
###########################################

## FIXME: this needs to change to use probeSetSummary() instead of hyperG2Affy().

hyperG2annaffy <- function(probids, lib, eset, fit = NULL, subset = NULL, comp = 1,
                           type = "MF", pvalue = 0.05, min.count = 10){
  tab <- hyperG2Affy(probids, lib, type, pvalue, min.count)
  if(!is.null(subset))
    tab <- tab[subset]
  tab <- lapply(tab, unique)
  for(i in seq(along = tab)){
    if(!is.null(fit)){
      index <- geneNames(eset) %in% tab[[i]]
      ord <- order(abs(fit$t)[index], decreasing = TRUE)
      otherdata <- list("t-statistic" = round(fit$t[,comp][index][ord], 2),
                        "p-value" = round(p.adjust(fit$p.value[,comp], "fdr")[index][ord], 3),
                        "Fold Change" = round(fit$coef[,comp][index][ord],2))
      probes2table(eset, tab[[i]][ord], annotation(eset), otherdata = otherdata,
                   filename = paste(Term(get(names(tab)[i], GOTERM)), "genes", sep=" "))
    }else{
      probes2table(eset, tab[[i]], annotation(eset),
                   filename = paste(Term(get(names(tab)[i], GOTERM)), "genes", sep = " "))
    }
  }
}
