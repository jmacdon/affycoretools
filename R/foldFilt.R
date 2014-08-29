##################################################
##
##  Copyright 2005 James W. MacDonald
##
##  foldFilt - functions to output HTML and text tables
##             based on fold change
##
##  2-14-2006 Modified to allow a filterfun from genefilter to be passed
##            to the function, allowing for filtering of each comparison
##
###################################################



#' Output Fold Change Data
#' 
#' This function is designed to take an \code{ExpressionSet} and some
#' comparisons and output either HTML tables, text files, or both.
#' 
#' This function is useful for outputting annotated gene lists for multiple
#' fold change comparisons. The genes will be ordered by the absolute fold
#' change. Note that this function is essentially a wrapper to call
#' \code{annaffy}, so is only useful for Affymetrix GeneChips for which there
#' is an annotation package.
#' 
#' Without attaching a data file to this package, it is not possible to give a
#' working example. Instead, here is a 'for instance'.
#' 
#' Say you have an \code{ExpressionSet} containing four Affy HG-U133Plus2
#' chips. There is no replication, and you simply want to output genes with a
#' two-fold or greater difference between the first chip and each of the last
#' three (the first chip is the control, and the other three are
#' experimentals). The \code{ExpressionSet} is called eset.
#' 
#' Additionally, say we don't want any genes called significant if both of the
#' samples have very low expression. We can set up a filter using the
#' \pkg{genefilter} package.
#' 
#' f1 <- kOverA(1,6)
#' 
#' filt <- filterfun(f1)
#' 
#' foldFilt(eset, groups=1:4, comps=list(c(2, 1), c(3, 1), c(4, 1)),
#' compnames=c("Expt1-Cont","Expt2-Cont","Expt3-Cont"), filterfun = filt)
#' 
#' This will output three HTML tables called 'Expt1-Cont.html', etc., each
#' containing sorted genes that have two-fold or greater differences between
#' the two samples.
#' 
#' @param object An \code{ExpressionSet} object
#' @param fold The log fold change cutoff to use. Note that this is log base
#' two.
#' @param groups A vector of group identifiers. Probably easiest to use a
#' numeric vector
#' @param comps A list containing all the comparisons to be made. Each list
#' item should be a vector of length two. See details for more information.
#' @param compnames A character vector of the names for each of the comparisons
#' to be made. This will be the name of the resulting HTML or text file.
#' @param save Boolean. If \code{TRUE}, a list will be returned. The first item
#' in the list will be a vector showing the number of 'significant' genes for
#' each comparison. The second item will be a matrix of -1's, 0's and 1's
#' indicating a significant difference, and the direction of the difference.
#' The first item is useful for creating Sweave - based reports and the second
#' is useful for making Vennn diagrams using the \code{vennDiagram} from the
#' limma package.
#' @param html Boolean - if \code{TRUE}, output HTML tables
#' @param text Boolean - if \code{TRUE}, output text tables
#' @param filterfun A filtering function, created by
#' \code{\link[genefilter]{genefilter}} to filter the data using additional
#' criteria. See details for more information
#' @return Returns a list; see above for the elements of the list. This
#' function is mainly called for the side effect of outputting HTML or text
#' files containing annotated 'significant' gene lists.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @keywords manip
#' @export foldFilt
foldFilt <-function(object, fold = 1, groups, comps, compnames,
                    save = FALSE, text = TRUE, html = TRUE, filterfun = NULL){
  .Deprecated(new = "", msg = paste("foldFilt is being deprecated. Please see the RefactoredAffycoretools",
                        "vignette for more up to date ways to annotate results."))
  if(is(object, "ExpressionSet"))
    x  <- exprs(object)
  if(length(unique(groups)) != length(groups)){
    gps <- matrix(NA, ncol = length(unique(groups)), nrow = dim(x)[1])
    for(i in unique(groups)){
      if(length(groups[which(groups == i)]) != 1)
        gps[,i] <- rowMeans(x[,which(groups == i)])
      else
        gps[,i] <- x[,which(groups == i)]
    }
  }else{
    gps <- x
  }
  colnames(gps) <- unique(groups)
  flds <- lapply(comps, function(y) gps[,y[1]] - gps[,y[2]])
  indices <- lapply(flds, function(y) abs(y) > fold)
  
  if(!is.null(filterfun)){
    filt.ind <- lapply(comps, function(y) genefilter(cbind(gps[,y[1]], gps[,y[2]]), filterfun))
    indices <- mapply(function(x, y) x * y, indices, filt.ind, SIMPLIFY = FALSE)
    indices <- lapply(indices, as.logical)
  }
  
  probes <- lapply(indices, function(y) row.names(x)[y])
  FCs <- vector("list", length(unique(groups)))
  for(i in seq(along=flds)){
    FCs[[i]] <- flds[[i]][indices[[i]]]
    if(length(FCs[[i]]) > 0){
      ord <- order(abs(FCs[[i]]), decreasing = TRUE)
      idx <- getIndex(comps[[i]], groups)
      if(html || text)
        probes2table(object[,idx], probes[[i]][ord], object@annotation,
                     otherdata = list("Fold Change" = FCs[[i]][ord]),
                     text = text, html = html, filename = compnames[i])
    }
  }
  direct <- matrix(NA, ncol = length(comps), nrow = dim(gps)[1],
                   dimnames = list(featureNames(object), compnames))
  for(i in seq(along = flds)){
    direct[,i] <- sign(flds[[i]] * indices[[i]])
  }
  if(save)
   return(list(sums =  sapply(indices, sum), dirs = direct))
}
        


#' Get Indices
#' 
#' This is an internal function and not intended to be called by end users.
#' 
#' 
#' @param comps A list of comparisons
#' @param groups A vector of groupnames
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @keywords internal
getIndex <- function(comps, groups){
  first <- which(groups == comps[1])
  second <- which(groups == comps[2])
  c(first, second)
}



#' Select and Output Gene Lists Based on Venn Diagrams
#' 
#' This function is designed to output text and/or HTML tables based on the
#' results of a call to \code{foldFilt}. The general idea being that one might
#' want to create a Venn diagram showing probesets that are unique to
#' particular comparisons, or consistent between comparisons, and then might
#' want to output the probesets that are contained in each cell of the Venn
#' diagram.
#' 
#' The purpose of this function is to output the probesets listed in a Venn
#' diagram that has been produced by a call to \code{foldFilt}. A small example
#' would be as follows:
#' 
#' Assume an \code{ExpressionSet} exists that contains expression values for
#' three Affymetrix chips, say a control, and two experimentals. One might want
#' to know what probesets are different between each of the experimentals and
#' the control, and those that are different between both of the experimentals
#' and the control. We first make the comparisons, based on a fold change of 2
#' (or a difference of 1 on the log scale).
#' 
#' comps <- list(c(1,2), c(1,3))
#' 
#' This list indicates what comparisons we want. In this case 1vs2 and 1vs3.
#' 
#' out <- foldFilt(eset, fold = 1, groups = 1:3, comps = comps,
#' compnames=c("Control vs experimental1", "Control vs experimental2"), save =
#' TRUE)
#' 
#' By setting save = TRUE, we are saving a list, the first item being a vector
#' of the number of probesets in each comparison, the second item being an
#' indicator matrix showing up or down regulation based on a two-fold
#' difference. We could make a Venn diagram using this matrix with
#' \code{vennCounts2} and \code{\link[limma:venn]{vennDiagram}}. If we then
#' wanted to output the probesets in each cell of that Venn diagram, we could
#' use \code{vennSelectFC} as follows:
#' 
#' vennSelectFC(eset, out[[2]], comps)
#' 
#' One thing to note here is that the names of the resulting tables as well as
#' the columns containing the fold change values will be extracted from the
#' column names of the indicator matrix. This matrix will get its column names
#' from the 'compnames' argument to \code{foldFilt}, so it is best to use
#' reasonable names here. Also note that any character used in the 'compnames'
#' argument that is not a valid character for a file name will be stripped out.
#' 
#' @param eset A \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}
#' object.
#' @param x An indicator matrix showing up or down regulation based on fold
#' change, usually from a call to \code{foldFilt}. See details for more
#' information.
#' @param comps A list containing all the comparisons to be made. Each list
#' item should be a vector of length two. This should be identical to the
#' 'comps' argument used in the call to \code{foldFilt} See details for more
#' information.
#' @param order.by One of 'sum', 'max', 'median', or 'mean'. This orders the
#' output for those tables that have multiple fold change values based on the
#' summary statistic chosen. Defaults to 'sum'.
#' @param method One of "same", "both", "up", "down", "sameup", or "samedown".
#' See details for more information.
#' @param text Boolean. Output text tables? Defaults to \code{TRUE}
#' @param html Boolean. Output HTML tables? Defaults to \code{TRUE}
#' @param ... Used to pass other arguments to \code{probes2table}, in
#' particular, to change the argument to \code{anncols} which controls the
#' columns of hyperlinks to online databases (e.g., Entrez Gene, etc.). See
#' \code{\link[annaffy]{aaf.handler}} for more information.
#' @return Called only for the side effect of outputting HTML and/or text
#' tables.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @keywords manip
#' @export vennSelectFC
vennSelectFC <- function(eset, x, comps, order.by = "sum", method = "same", text = TRUE,
                         html = TRUE, ...){
  .Deprecated(new = "makeVenn", msg = paste("vennSelectFC is being deprecated. Please see the RefactoredAffycoretools",
              "vignette for more current ways to create Venn diagrams."))
  ncomps <- dim(x)[2]
  if(ncomps < 2 || ncomps > 3)
    stop("This function only works for two or three comparisons at a time.\n",
         call. = FALSE)
  if (ncomps == 2)
    name <- c(paste("Genes unique to", colnames(x)),
              paste("Genes in intersection of",  intNames(x)))
  else
    if (ncomps == 3) 
      name <- c(paste("Genes unique to", colnames(x)),
                paste("Genes in intersection of", intNames(x)),
                "Genes common to all comparisons")
  indices <- makeIndices(x, method)
  tmp <- makeFCList(x, indices, comps, exprs(eset))
  fcList <- tmp$out
  ord <- getOrder(fcList, order.by)
  ## order the fold change list
  for(i in seq(along = fcList)){
    if(length(fcList[[i]]) == 1) fcList[[i]][[1]] <- fcList[[i]][[1]][ord[[i]]]
    else
      for(j in seq(along = fcList[[i]])) fcList[[i]][[j]] <- fcList[[i]][[j]][ord[[i]]]
  }
  comps <- lapply(tmp$comps, function(x) if(is.list(x)) unique(unlist(x)) else x)
  for(i in seq(along = fcList)){
    ## skip empty cells  
    if(length(fcList[[i]][[1]]) == 0)
        next
    tmp <- featureNames(eset)[indices[[i]]][ord[[i]]]
    otherdata <- fcList[[i]]
    probes2table(eset[,comps[[i]]], tmp, annotation(eset), otherdata, html = html,
                 text = text, express = TRUE, filename = name[i])
  }
}




#' Make a Fold Change List
#' 
#' This function is used to create a \code{list} containing named lists of fold
#' change values suitable to be passed to \code{probes2table}, and ultimately
#' \code{\link[annaffy]{aafTableAnn}} in the annaffy package. This is an
#' internal function and is not intended to be called by end users.
#' 
#' The purpose of this function is to create a list of named lists that will be
#' suitable to be passed to \code{probes2table}, which will result in a column
#' of fold change values with a reasonable name when output as an HTML or text
#' table using annaffy.
#' 
#' @param x An indicator matrix showing up or down regulation based on fold
#' change, usually from a call to \code{foldFilt}. See details for more
#' information.
#' @param indices A \code{list} created by a call to \code{makeIndices}.
#' @param comps A list containing all the comparisons to be made. Each list
#' item should be a vector of length two. This should be identical to the
#' 'comps' argument used in the call to \code{foldFilt}
#' @param dat A \code{matrix} containing the expression values, usually
#' extracted using the \code{exprs} accessor.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @keywords internal
makeFCList <- function(x, indices, comps, dat){
  ## Make list of lists containing fold change values
  if(length(comps) == 2) comps <- c(comps, list(comps))
  else
    if(length(comps) == 3) comps <- c(comps, list(comps[1:2], comps[c(1,3)], comps[2:3]),
               list(comps))
  out <- mapply(function(x, y){
    if(is.vector(x) && !is.list(x)) list(dat[y, x[1]] - dat[y, x[2]])
    else
      lapply(x, function(x) dat[y, x[1]] - dat[y, x[2]])}, comps, indices)
  ## Make list of names to use as names for output list
  nam <- paste("Fold change", colnames(x))
  if(length(nam) == 2)
    tmp <- c(as.list(nam), list(nam))
  if(length(nam) == 3)
    tmp <- c(as.list(nam), list(nam[1:2], nam[c(1,3)], nam[2:3]), list(nam))
  for(i in seq(along = out)) names(out[[i]]) <- tmp[[i]]
  return(list(out = out, comps = comps))
}




#' Order Probesets Based on Fold Change Values
#' 
#' This function is used to create a \code{list} containing the correct order
#' to output a set of probesets based on the fold change values for those
#' probesets. This is an internal function and is not intended to be called by
#' end users.
#' 
#' The purpose of this function is to come up with some ordering for probesets
#' that will be ouput by a call to \code{vennSelectFC}. The ordering for
#' probesets unique to a given comparison is simple, since there is only one
#' fold change to be considered. However, when there are more than one
#' comparison, we have to decide how the probesets should be ordered. The
#' choices available include those that seemed reasonable to me when I wrote
#' this function. Presumably there are other reasonable choices that may need
#' to be included.
#' 
#' @param fcList A \code{list} of named lists, generally produced by a call to
#' \code{makeFCList}
#' @param order.by One of 'sum', 'max', 'median', or 'mean'. This orders the
#' output for those tables that have multiple fold change values based on the
#' summary statistic chosen. Defaults to 'sum'.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @keywords internal
getOrder <- function(fcList, order.by){
  ord <- lapply(fcList, function(x){
    if(length(x) == 1) order(abs(x[[1]]), decreasing = TRUE)
    else
      if(length(x) > 1) switch(order.by,
                               sum = order(abs(rowSums(data.frame(x))), decreasing = TRUE),
                               max = order(abs(apply(data.frame(x), 1, max)),
                                 decreasing = TRUE),
                               mean = order(abs(rowMeans(data.frame(x))), decreasing = TRUE),
                               median = order(abs(apply(data.frame(x), 1, median)),
                                 decreasing = TRUE))})
}
