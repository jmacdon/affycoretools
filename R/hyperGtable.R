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




#' HTML tables from GOIDs
#' 
#' Output HTML tables containing the 'enriched' genes for each GO term
#' resulting from a call to \code{hyperGtable}.
#' 
#' This function is used to create HTML tables based on the output of
#' \code{hyperGtable}. The basic idea is as follows; as part of an analysis,
#' say \code{hyperGtable} was used to create a table of 'enriched' GO terms.
#' Unfortunately, the table only lists GO terms and the number of probesets
#' that are annotated to those GO terms, and the client may be interested in
#' knowing what probesets are enriched for each (or some) GO term.
#' 
#' The default behaviour is to output an HTML table for each GO term,
#' containing the probesets that are annotated at that term (and that are in
#' the set of significant genes). If only some of the GO terms are of interest,
#' one may use the \code{subset} argument to select only particular rows. In
#' addition, if the relevant t-statistics, p-values and fold changes are of
#' interest, one can also use the \code{fit} argument to point to an
#' \code{\link[limma]{lmFit}} object that contains these data, as well as the
#' \code{comp} argument to indicate which parameter or contrast to use. Note
#' that the \code{comp} argument defaults to 1, so the first parameter or
#' contrast will be extracted by default.
#' 
#' @param probids A vector of Affymetrix probe IDs
#' @param lib An annotation package (e.g., \code{hgu95av2})
#' @param eset An \code{ExpressionSet}
#' @param fit An \code{\link[limma]{lmFit}} object. Only necessary if
#' statistics are desired in the resulting table. Defaults to \code{NULL}.
#' @param subset A numeric vector used to select GO terms to output (see
#' description for more information). Defaults to \code{NULL}
#' @param comp Which contrast/parameter estimate should be used to extract the
#' relevant statistics? Only used if \code{fit} is not \code{NULL}. See
#' description for more information.
#' @param type One of "MF", "CC", "BP", indicating molecular function, cellular
#' component, or biological process, respectively.
#' @param pvalue The significance level used to choose GO terms
#' @param min.count The minimum number of a given GO term that must be on the
#' chip in order to choose that GO term. This protects against very low
#' p-values that result from the situation where there are very few genes with
#' a given GO term on the chip, but one or two are found in the set of
#' significant genes.
#' @return This function is used only for the side effect of creating HTML
#' tables.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @keywords manip htest
#' @export hyperG2annaffy
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



#' Output Tables Based on Hypergeometric Test
#' 
#' This function will output various tables containing probesets that are
#' annotated to a particular GO, KEGG, or PFAM term. The tables are based on
#' the results from a call to hyperGtest.
#' 
#' This function is designed to be used to output the results from a
#' hypergeometric test for over-represented terms. This function would be used
#' at the end of an analysis such as:
#' 
#' 1.) Compute expression values 2.) Fit a model using \code{limma} 3.) Output
#' significant probesets using \code{limma2annaffy} 4.) Perform hypergeometric
#' test using \code{\link[Category:HyperGResult-accessors]{hyperGTest}}
#' 
#' At step 4, one can output a list of the over-represented terms using
#' \code{\link[Category:HyperGResult-accessors]{htmlReport}}. One might then be
#' interested in knowing which probesets contributed to the significance of a
#' particular term, which is what this function is designed to do.
#' 
#' One argument that can be passed to
#' \code{\link[Category:HyperGResult-accessors]{htmlReport}} (and also to
#' \code{hyperGoutput}) is \code{categorySize}, which gives a lower bound for
#' the number of probesets with a particular term in the universe. In other
#' words, assume that a particular GO term is annotated to three probesets on a
#' given chip. If, after doing a t-test to detect differentially expressed
#' probesets, one of those probesets were found to be significantly
#' differentially expressed and was then used to do a hypergeometric test, that
#' GO term would be significant, with a small p-value. However, this is
#' probably not very strong evidence that the GO term is actually
#' over-represented, since there were only three to begin with. By setting
#' \code{categorySize} to a sensible value (such as 10), this situation can be
#' avoided.
#' 
#' This function will output HTML and/or text tables containing annotation
#' information about each probeset as well as the expression values. In
#' addition, if limma were used to fit the model, the relevant statistics
#' (t-statistic, p-value, fold change) can also be output in the table by
#' passing the \code{\link[limma:marraylm]{MArrayLM}} object that resulted from
#' a
#' 
#' call to \code{\link[limma:ebayes]{eBayes}}. The \code{statistics} argument
#' can
#' 
#' be used to control which statistics are output.
#' 
#' By default \code{hyperGoutput} will output tables for all significant terms,
#' which may end up being quite a few tables. Usually only a few terms are of
#' interest, so there is a \code{subset} argument that can be used to select
#' only those terms. This argument follows directly from the order of the table
#' output by \code{\link[Category:HyperGResult-accessors]{htmlReport}} or
#' \code{\link[GOstats:GOHyperGResult-class]{summary}}. For instance, if the
#' first, third and fifth terms in the HTML table output by
#' \code{\link[Category:HyperGResult-accessors]{htmlReport}} were of interest,
#' one would use subset=c(1,3,5).
#' 
#' One critical step prior to the hypergeometric test is to subset the
#' probesets to unique Entrez Gene IDs. It should be noted however, that the
#' functions used by \code{hypergOutput} will output all the probesets
#' annotated to a particular term. The \code{output} argument is used to
#' control this behavior. If output = "significant" (the default), then only
#' those probesets that correspond to the original subsetting will be output.
#' If output = "all", then all probesets will be output (grouped by Entrez ID),
#' with the 'significant' probeset first. If output = "split", then all the
#' probesets will be output, with all the 'significant' probesets first,
#' followed by the other probesets, grouped by Entrez ID.
#' 
#' Note that the 'significant' probesets come from one of two sources. First,
#' one can pass a character vector of probeset IDs corresponding to those that
#' were significant in the original analysis (recommended). Second, if the
#' \code{geneIds} slot of the \code{GOHyperGParams} object containes a named
#' vector of Entrez Gene IDs, then the names from that vector will be used.
#' This can be accomplished by using either
#' \code{\link[genefilter]{findLargest}} or \code{getUniqueLL}.
#' 
#' Since the \code{geneIds} are by definition a unique set of Entrez Gene IDs,
#' any duplicate probeset IDs will have been removed, so the first method is to
#' be preferred for accuracy.
#' 
#' @param hyptObj A \code{HyperGResult} object, usually produced by a call to
#' \code{\link[Category]{hyperGTest}}
#' @param eset An \code{ExpressionSet} object
#' @param pvalue The p-value cutoff used for selecting significant GO terms. If
#' not specified, it will be extracted from the \code{HyperGResult} object
#' @param categorySize Number of terms in the universe required for a term to
#' be significant. See details for more information
#' @param sigProbesets Vector of probeset IDs that were significant in the
#' original analysis.
#' @param fit An \code{\link[limma:marraylm]{MArrayLM}} object, produced from a
#' call to \code{\link[limma:ebayes]{eBayes}}
#' @param subset Numeric vector used to select particular tables to output. The
#' default is to output tables for all terms. See details for more information
#' @param comp Numeric vector of length one, used to indicate which comparison
#' in the \code{\link[limma:marraylm]{MArrayLM}} object to use for extracting
#' relevant statistics. See details for more information
#' @param output One of 'selected', 'all', or 'split'. See details for more
#' information
#' @param statistics Which statistics to output in the resulting tables.
#' Choices include 'tstat', 'pval', or 'FC', corresponding to t-statistics,
#' p-values, and fold change, respectively
#' @param html Boolean. Output HTML tables? Defaults to \code{TRUE}
#' @param text Boolean. Output text tables? Defaults to \code{TRUE}
#' @param \dots Allows end user to pass further arguments. The most notable
#' would be an \code{anncols} argument, passed to \code{probes2table} to
#' control the hyperlinked annotation columns. See
#' \code{\link[annaffy]{aaf.handler}} for more information
#' @return This function returns no value, and is called solely for the side
#' effect of outputting HTML and/or text tables.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @seealso \code{\link[Category]{hyperGTest}},
#' \code{\link[Category:HyperGResult-accessors]{htmlReport}},
#' \code{\link[GOstats]{probeSetSummary}}
#' @keywords manip
#' @export hyperGoutput
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



#' Internal functions for hyperGoutput
#' 
#' These functions are internal functions for \code{hyperGoutput} and are not
#' intended to be called directly by the end user.
#' 
#' These functions control the output and ordering of the probesets for the
#' \code{hyperGoutput} function.
#' 
#' @aliases houtSel houtAll houtSplit
#' @param tab A \code{data.frame}, resulting from a call to
#' \code{\link[GOstats]{probeSetSummary}}
#' @param eset An \code{link[Biobase:class.ExpressionSet]{ExpressionSet}}
#' object
#' @param fit An \code{\link[limma:marraylm]{MArrayLM}} object, resulting from
#' a call to \code{\link[limma:ebayes]{eBayes}}
#' @param comp Numeric vector of length one, designating which comparison to
#' select from the \code{\link[limma:marraylm]{MArrayLM}} object
#' @param statistics Any or all of 'tstat', 'pval', 'FC', used to control the
#' output of statistics from the \code{\link[limma:marraylm]{MArrayLM}} object
#' @param html Boolean. Output HTML tables? Defaults to \code{TRUE}
#' @param text Boolean. Output text tables? Defaults to \code{TRUE}
#' @return Nothing is returned. These functions are called solely for their
#' side effect, the output of HTML and/or text tables.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @keywords internal
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



#' Extract Statistics from an MArrayLM object
#' 
#' This is a utility function that can be used to extract the statistics for a
#' set of probesets from an \code{\link[limma:marraylm]{MArrayLM}} object. This
#' may be useful for others, but currently is considered an internal function
#' and is not intended to be called by end users.
#' 
#' This is an internal function used to return a named list of statistics,
#' which can then be passed to \code{probes2table} for creating HTML and/or
#' text tables.
#' 
#' @param prbs A vector of probe IDs
#' @param fit An \code{\link[limma:marraylm]{MArrayLM}} object
#' @param comp Numeric vector of length one, designating which comparison to
#' select from the \code{\link[limma:marraylm]{MArrayLM}} object
#' @param statistics Any or all of 'tstat', 'pval', 'FC', used to control the
#' output of statistics from the \code{\link[limma:marraylm]{MArrayLM}} object
#' @return A named list that can have up to three vectors of statistics.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @keywords internal
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
  


#' Subset a Vector of Probesets
#' 
#' This function will take a vector of Affy IDs and return a vector of Entrez
#' IDs that have replicated IDs removed. The resulting vector will still have
#' the corresponding Affy IDs appended as names, which is important for some
#' functions.
#' 
#' Subsetting a set of Affy IDs to unique Entrez Gene IDs is a common thing to
#' do prior to doing a hypergeometric test. Functions such as
#' \code{\link[Category]{hyperGTest}} can use un-named vectors of Entrez IDs
#' (e.g., unique(getLL(probeIDs, annot))), but there is some functionality that
#' requires the Entrez Gene IDs to be in a named vector, with the names being
#' the associated Probeset IDs.
#' 
#' As an example, \code{hyperGoutput} will only work correctly if the input
#' Entrez ID vector is named with the associated Probeset IDs.
#' 
#' @param probes A vector of probe IDs
#' @param annot The annotation package for the chip used
#' @return A named vector of unique Entrez IDs
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @keywords manip
#' @export getUniqueLL
getUniqueLL <- function(probes, annot){

  out <- getEG(probes, annot)
  out <- out[match(unique(out), out)]
  out <- out[!is.na(out)]
  out
}
