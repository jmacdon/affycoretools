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
##
##  hyperG2annaffy defunct as of 12-13-08
###########################################


hyperG2annaffy <- function(probids, lib, eset, fit = NULL, subset = NULL, comp = 1,
                           type = "MF", pvalue = 0.05, min.count = 10){
  .Defunct("hyperGoutput", msg = "hyperG2annaffy is defunct. Use hyperGoutput instead.\n")
##  tab <- hyperG2Affy(probids, lib, type, pvalue, min.count)
##   if(!is.null(subset))
##     tab <- tab[subset]
##   tab <- lapply(tab, unique)
##   for(i in seq(along = tab)){
##     if(!is.null(fit)){
##       index <- featureNames(eset) %in% tab[[i]]
##       ord <- order(abs(fit$t)[index], decreasing = TRUE)
##       otherdata <- list("t-statistic" = round(fit$t[,comp][index][ord], 2),
##                         "p-value" = round(p.adjust(fit$p.value[,comp], "fdr")[index][ord], 3),
##                         "Fold Change" = round(fit$coef[,comp][index][ord],2))
##       probes2table(eset, tab[[i]][ord], annotation(eset), otherdata = otherdata,
##                    filename = paste(Term(get(names(tab)[i], GOTERM)), "genes", sep=" "))
##     }else{
##       probes2table(eset, tab[[i]], annotation(eset),
##                    filename = paste(Term(get(names(tab)[i], GOTERM)), "genes", sep = " "))
##     }
##  }
}

hyperGoutput <- function(hyptObj, eset, pvalue, categorySize, sigProbesets, fit = NULL,
                         subset = NULL, comp = 1, output = c("significant", "all", "split"),
                         statistics = c("tstat","pval","FC"), html = TRUE, text = TRUE, ...){
  require(paste(annotation(hyptObj), "db", sep = "."), character.only = TRUE, quietly = TRUE)
  if (!is(hyptObj, "GOHyperGResult")) 
    stop("result must be a GOHyperGResult instance (or subclass)")
  if(!all(output %in% c("significant", "all", "split")))
    stop(paste(output, "is not a valid choice!",
               "Please choose from 'significant', 'all', and 'split'.", call. = FALSE))
  tmp <- probeSetSummary(result=hyptObj, pvalue=pvalue, categorySize=categorySize,
                         sigProbesets=sigProbesets)
  if(!is.null(subset))
    tmp <- tmp[subset]
  output <- match.arg(output)
  for(i in seq(along = tmp)){
    switch(output,
           significant = houtSel(tmp[i], eset, fit, comp, statistics, html, text),
           all = houtAll(tmp[i], eset, fit, comp, statistics, html, text),
           split = houtSplit(tmp[i], eset, fit, comp, statistics, html, text))
  }
}

houtSel <- function(tab, eset, fit, comp, statistics, html, text){
  nam <- names(tab)
  tab <- tab[[1]]
  index <- tab[,3] == 1
  prbs <- tab[index,2]
  if(!is.null(fit)){
    ord <- order(abs(fit$t[prbs,comp]), decreasing = TRUE)
    otherdata <- extractStats(prbs[ord], fit, comp, statistics)
    probes2table(eset, prbs[ord], annotation(eset), otherdata, html = html, text = text,
                 filename = paste("Probesets annotated to", Term(get(nam, GOTERM))))
  }else{
    probes2table(eset, prbs, annotation(eset), html = html, text = text,
                 filename =  paste("Probesets annotated to", Term(get(nam, GOTERM))))
  }
}

houtAll <- function(tab, eset, fit, comp, statistics, html, text){
  nam <- names(tab)
  tab <- tab[[1]]
  ## subset tab to those under consideration
  index <- tab[,2] %in% featureNames(eset)
  tab <- tab[index,]
  tmpprbs <- tab[tab[,3] == 1,2]
  if(!is.null(fit)){
    ## First order the selected probesets, then the unselected ones
    ord <- order(abs(fit$t[tmpprbs,comp]), decreasing = TRUE)
    oprbs <- tmpprbs[ord]
    oegids <- tab[match(oprbs, tab[,2]),1]
    smtab <- tab[tab[,3] == 0,]
    prbs <- vector()
    for(i in seq(along = oegids)){
      index <- smtab[,1] %in% oegids[i]
      ord <- order(abs(fit$t[smtab[index,2],comp]), decreasing = TRUE)
      prbs <- c(prbs, oprbs[i], smtab[index,2][ord])
    }
    otherdata <- extractStats(prbs, fit, comp, statistics)
    probes2table(eset, prbs, annotation(eset), otherdata, html = html, text = text,
                 filename = paste("Probesets annotated to", Term(get(nam, GOTERM))))
  }else{
    ord <- order(tab[,1], -tab[,3])
    prbs <- tab[,2][ord]
    probes2table(eset, prbs, annotation(eset), html = html, text = text,
                 filename =  paste("Probesets annotated to", Term(get(nam, GOTERM))))
  }
}

houtSplit <- function(tab, eset, fit, comp, statistics, html, text){
  nam <- names(tab)
  tab <- tab[[1]]
  ## subset tab to those under consideration
  index <- tab[,2] %in% featureNames(eset)
  tab <- tab[index,]
  tmpprbs <- tab[tab[,3] == 1,2]
  if(!is.null(fit)){
    ord <- order(abs(fit$t[tmpprbs,comp]), decreasing = TRUE)
    oprbs <- tmpprbs[ord]
    oegids <- tab[match(oprbs, tab[,2]),1]
    smtab <- tab[tab[,3] == 0,]
    prbs <- oprbs
    for(i in seq(along = oegids)){
      index <- smtab[,1] %in% oegids[i]
      ord <- order(abs(fit$t[smtab[index,2],comp]), decreasing = TRUE)
      prbs <- c(prbs, smtab[index,2][ord])
    }
    otherdata <- extractStats(prbs, fit, comp, statistics)
    probes2table(eset, prbs, annotation(eset), otherdata, html = html, text = text,
                 filename = paste("Probesets annotated to", Term(get(nam, GOTERM))))
  }else{
    ord <- order(-tab[,3], tab[,1])
    prbs <- tab[,2][ord]
    probes2table(eset, prbs, annotation(eset), html = html, text = text,
                 filename =  paste("Probesets annotated to", Term(get(nam, GOTERM))))
  }
}

extractStats <- function(prbs, fit, comp, statistics){
  if(!all(statistics %in% c("tstat","pval","FC")))
    stop(paste("The 'statistics' supplied (", statistics,") are not correct\n",
               "valid values are 'tstat', 'pval', 'FC'.", call. = FALSE))
  statistics <- as.list(statistics)
  out <- lapply(statistics, function(x) {
    switch(x,
           tstat = round(fit$t[prbs,comp], 2),
           pval = round(fit$p.value[prbs,comp],3),
           FC = round(fit$coef[prbs,comp],2))})
  nam <- vector()
  for(i in seq(along = statistics))
    nam[i] <- switch(statistics[[i]],
                     tstat = "t-statistic",
                     pval = "p-value",
                     FC = "Fold change")
  names(out) <- nam
  out
}
  
getUniqueLL <- function(probes, annot){

  out <- getEG(probes, annot)
  out <- out[match(unique(out), out)]
  out <- out[!is.na(out)]
  out
}
