################################################
##
##  Copyright 2005 James W. MacDonald
##
##  Functions to both create Venn diagrams
##  and to output the subsets from the Venn
##  diagrams in HTML and text tables
##
################################################




#' Compute Counts for Venn Diagram
#' 
#' This function is designed to compute counts for a Venn diagram. It is
#' slightly different from \code{\link[limma:venn]{vennCounts}} in the
#' additional ability to compute counts for genes that are differentially
#' expressed in the same direction.
#' 
#' The function \code{\link[limma:venn]{vennCounts}} will return identical
#' results except for the "same" method. This will only select those genes that
#' both pass the criteria of \code{\link[limma]{decideTests}} as well as being
#' differentially expressed in the same direction. Note that this is different
#' from the "both" method, which simply requires that a given gene be
#' differentially expressed in e.g., two different comparisons without any
#' requirement that the direction be the same.
#' 
#' @param x A \code{\link[limma]{TestResults}} object, produced by a call to
#' \code{\link[limma]{decideTests}} or \code{foldFilt}.
#' @param method One of "same", "both", "up", "down". See details for more
#' information.
#' @param fit An \code{\link[limma:marraylm]{MArrayLM}} object, produced by a
#' call
#' 
#' to \code{\link[limma]{lmFit}} and \code{\link[limma:ebayes]{eBayes}}. Only
#' necessary if 'foldFilt' = \code{TRUE}.
#' @param foldFilt A fold change to filter samples. This is primarily here for
#' consistency with the corresponding argument in \code{vennSelect}.
#' @return A \code{\link[limma:venn]{VennCounts}} object.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @keywords manip hplot
#' @examples
#' 
#' library("limma")
#' tstat <- matrix(rt(300,df=10),100,3)
#' tstat[1:33,] <- tstat[1:33,]+2
#' clas <- classifyTestsF(tstat,df=10,p.value=0.05)
#' a <- vennCounts2(clas)
#' print(a)
#' vennDiagram(a)
#' 
#' @export vennCounts2
vennCounts2 <- function(x, method = "same", fit = NULL,
                        foldFilt = NULL){

  ## x is a TestResults object from a call to
  ## decideTests()
  ## alternatively it can be the 'dirs' matrix from
  ## a call to foldFilt()
  ## fit is an MArrayLM object and foldFilt is a fold change to
  ## filter on. Note that fit is only required if foldFilt != NULL

  if(!is.null(foldFilt) && is.null(fit))
    stop("If you want to filter by fold change, you must also pass an MArrayLM object\n",
         "that was produced from a call to lmFit() and then eBayes().", call. = FALSE)
  if(!is.null(foldFilt)){
        idx <- abs(fit$coefficients) > foldFilt
        x <- as.matrix(x) * idx
    }
  tmp <- vennCounts(x)
  newcnts <- sapply(makeIndices(x, method), sum)
  all <- sum(tmp[,dim(tmp)[2]])
  if(dim(x)[2] == 2)
    ord <- c(2, 1, 3)
  if(dim(x)[2] == 3)
    ord <- c(3, 2, 6, 1, 5, 4, 7)
  tmp[,dim(tmp)[2]] <- c(all - sum(newcnts), newcnts[ord])
  tmp
}



#' Create Names for Venn Diagram Intersections
#' 
#' This function is designed to create the names for all the intersections of a
#' Venn diagram based on the names of all the single comparisons. This is an
#' internal function and is not intended to be called by the end user.
#' 
#' 
#' @param x A \code{TestResults} object from a call to \code{decideTests}.
#' @return A vector of names.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @keywords internal
intNames <- function(x){
  tmp <- colnames(x)
  if(is.null(dim(x)) || dim(x)[2] > 3)
    stop("This function only works for two or three comparisons at a time\n",
         call. = FALSE)
  if(dim(x)[2] == 2)
    return(paste(tmp[1], "and", tmp[2], sep = " "))
  if(dim(x)[2] == 3)
    return(paste(c(tmp[1], tmp[1], tmp[2]), "and", c(tmp[2], tmp[3], tmp[3])))
}



#' Create Indices for Venn Diagrams
#' 
#' This function is used to create indices for making Venn diagrams or
#' \code{vennCounts} objects. This is an internal method and is not intended to
#' be called by the end user.
#' 
#' 
#' @param x A \code{TestResults} object from a call to \code{decideTests}.
#' @param method One of "same", "up", "down", "both", indicating direction of
#' differential expression. See details of \code{vennSelect} for more
#' information.
#' @return A list containing \code{TRUE} and \code{FALSE} values, to be used by
#' \code{vennSelect}.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @keywords internal
makeIndices <- function(x, method = "same"){
  method <- match.arg(method, c("same", "up", "down", "both", "sameup", "samedown"))
  indices <- list()
  x <- sign(switch(method, both = abs(x), up = x > 0, down = x < 0, same = x,
                   sameup = x, samedown = x))
  if(ncol(x) == 2){
    if(method != "sameup" && method != "samedown"){
      indices[[1]] <- abs(x[,1]) == 1 & x[,2] != x[,1]
      indices[[2]] <- x[,1] != x[,2] & abs(x[,2]) == 1
      indices[[3]] <- abs(rowSums(x)) == 2
    }
    if(method == "sameup"){
      indices[[1]] <- x[,1] == 1 & x[,2] != x[,1]
      indices[[2]] <- x[,1] != x[,2] & x[,2] == 1
      indices[[3]] <- rowSums(x) == 2
    }
    if(method == "samedown"){
      indices[[1]] <- x[,1] == -1 & x[,2] != x[,1]
      indices[[2]] <- x[,1] != x[,2] & x[,2] == -1
      indices[[3]] <- rowSums(x) == -2
    }
  }
  if(ncol(x) == 3){
    if(method != "sameup" && method != "samedown"){
      indices[[1]] <- abs(x[,1]) == 1 & x[,2] != x[,1] & x[,3] != x[,1]
      indices[[2]] <- x[,1] != x[,2] & abs(x[,2]) == 1 & x[,3] != x[,2]
      indices[[3]] <- x[,1] != x[,3] & x[,2] != x[,3] & abs(x[,3]) == 1
      indices[[4]] <- abs(rowSums(x[,1:2])) == 2 & x[,3] != x[,1]
      indices[[5]] <- abs(rowSums(x[,c(1,3)])) == 2 & x[,2] != x[,1]
      indices[[6]] <- abs(rowSums(x[,2:3])) == 2 & x[,1] != x[,2]
      indices[[7]] <- abs(rowSums(x)) == 3
    }
    if(method == "sameup"){
      indices[[1]] <- x[,1] == 1 & x[,2] != x[,1] & x[,3] != x[,1]
      indices[[2]] <- x[,1] != x[,2] & x[,2] == 1 & x[,3] != x[,2]
      indices[[3]] <- x[,1] != x[,3] & x[,2] != x[,3] & x[,3] == 1
      indices[[4]] <- rowSums(x[,1:2]) == 2 & x[,3] != x[,1]
      indices[[5]] <- rowSums(x[,c(1,3)]) == 2 & x[,2] != x[,1]
      indices[[6]] <- rowSums(x[,2:3]) == 2 & x[,1] != x[,2]
      indices[[7]] <- rowSums(x) == 3
    }
    if(method == "samedown"){
      indices[[1]] <- x[,1] == -1 & x[,2] != x[,1] & x[,3] != x[,1]
      indices[[2]] <- x[,1] != x[,2] & x[,2] == -1 & x[,3] != x[,2]
      indices[[3]] <- x[,1] != x[,3] & x[,2] != x[,3] & x[,3] == -1
      indices[[4]] <- rowSums(x[,1:2]) == -2 & x[,3] != x[,1]
      indices[[5]] <- rowSums(x[,c(1,3)]) == -2 & x[,2] != x[,1]
      indices[[6]] <- rowSums(x[,2:3]) == -2 & x[,1] != x[,2]
      indices[[7]] <- rowSums(x) == -3
    }
  }
  indices
}
      



#' Correct Ordering of Contrasts
#' 
#' A function to determine the correct ordering of contrasts for
#' \code{vennSelect}. This is an internal function and should not be called by
#' the end user
#' 
#' 
#' @param design A design matrix, usually produced by a call to
#' \code{model.matrix}.
#' @param contrast A contrasts matrix, either produced by hand or as a result
#' of a call to \code{\link[limma]{makeContrasts}}
#' @return This function returns the correct order for the contrasts.
#' @author James W. MacDonald
#' @keywords internal
getCols <- function(design, contrast){
  ## A function to get the correct columns of the ExpressionSet
  ## to output based on the design and contrasts matrices
  ## This is sort of a kludge and could probably be improved
  ncontrasts <- ncol(contrast)
  ids <- 1:dim(design)[1]
  tmp <- function(design, contrast, x){
    a <- design[,contrast[,x] > 0]
    if(is.matrix(a))
      a <- rowSums(abs(a)) > 0
    b <- design[,contrast[,x] < 0]
    if(is.matrix(b))
      b <- rowSums(abs(b)) > 0
    c(ids[as.logical(a)], ids[as.logical(b)])
  }
  cols <- list()
  if(ncontrasts == 2){
    for(i in 1:2)
      cols[[i]] <- tmp(design, contrast, i)
    cols[[3]] <- unique(unlist(cols[1:2]))
  }else{
    for(i in 1:3)
      cols[[i]] <- tmp(design, contrast, i)
    comb <- list(c(1, 2), c(1, 3), c(2, 3), 1:3)
    cols[4:7] <- lapply(comb, function(x) unique(unlist(cols[x])))
  }
  cols
}

getStat <- function(stat, ncontrasts, num, fit, contrast, index, adj.meth){

  ## function to get the correct statistics based on the cell of the Venn diagram

  if(ncontrasts == 2)
    tmp <- list(1, 2, 1:2)
  else
    tmp <- list(1, 2, 3, 1:2, c(1, 3), 2:3, 1:3)
  if(stat == "tstat")
    if(length(tmp[[num]]) == 1){
      out <- list("t-statistic" = fit$t[index,tmp[[num]]])
    }else{
      out <- fit$t[index,tmp[[num]]]
      if(is.matrix(out)){
          out <- as.list(as.data.frame(out))
      }else{
          out <- as.list(out)
      }
      names(out) <- paste("t-statistic for ", colnames(contrast)[tmp[[num]]], sep = "")
    }
  if(stat == "pval")
    if(length(tmp[[num]]) == 1){
      out <- list("p-value" = p.adjust(fit$p.value[,tmp[[num]]], adj.meth)[index])
    }else{
      out <- fit$p.value[,tmp[[num]]]
      if(is.matrix(out)){
          out <- as.list(as.data.frame(out))
      }else{
          out <- as.list(out)
      }
      out <- lapply(out, p.adjust, method = adj.meth)
      out <- lapply(out, function(x) x[index])
      names(out) <- paste("p-value for ", colnames(contrast)[tmp[[num]]], sep = "")
    }
  if(stat == "FC")
    if(length(tmp[[num]]) == 1){
      out <- list("log2 fold change" = fit$coefficients[index,tmp[[num]]])
    }else{
      out <- fit$coefficients[index,tmp[[num]], drop = FALSE]
      if(is.matrix(out)){
          out <- as.list(as.data.frame(out))
      }else{
          out <- as.list(out)
      }
      names(out) <- paste("log2 fold change for ", colnames(contrast)[tmp[[num]]], sep = "")
    }
  
  out
}


f.stat.dat <- function(fit, index, contrast, adj.meth = "BH", stats,
                       order.by = "pval", ncontrasts, num){

  ## A function to get the f-statistics, fold change, and p-value (or some subset)
  ## for each cell of the Venn diagram
 
  if(!order.by %in% stats)
    stop(paste("If you choose to order by ", order.by, " you must also include that as ",
               "one of the statistics to output.",sep = ""), call. = FALSE)
  
  if("fstat" %in% stats)
    f.stat <- list("F-statistic" = fit$F[index])
  if("pval" %in% stats){
    p.val <- list("p-value" = p.adjust(fit$F.p.value, method = adj.meth)[index])
    if(exists("f.stat")) out <- c(f.stat, p.val) else out <- p.val
  }
  if("FC" %in% stats){
    fc <- getStat("FC", ncontrasts, num, fit, contrast, index)
    if(exists("out")) out <- c(out, fc) else out <- fc
  }
  
  ord <- switch(order.by, pval = order(p.val[[1]]),
                fstat = order(f.stat[[1]], decreasing = TRUE),
                fc = order(fc[[1]], decreasing = TRUE))
  out <- lapply(out, function(x) x[ord])
  out <- lapply(out, function(x)
                ifelse(abs(x) < 1e-3, sprintf("%0.3e", x), sprintf("%0.3f", x)))
  return(list(out = out, ord = ord))
}

t.stat.dat <- function(fit, index, contrast, adj.meth = "BH", stats,
                       order.by = "pval", ncontrasts, num){
  ## A function to get the t-statistics, fold change, and p-value (or some subset)
  ## for each cell of the Venn diagram
 
  if(!order.by %in% stats)
    stop(paste("If you choose to order by ", order.by, " you must also include that as ",
               "one of the statistics to output.",sep = ""), call. = FALSE)
  
  if("tstat" %in% stats)
    t.stat <- getStat("tstat", ncontrasts, num, fit, contrast, index)
  if("pval" %in% stats){
    p.val <- getStat("pval", ncontrasts, num, fit, contrast, index, adj.meth)
    if(exists("t.stat")) out <- c(t.stat, p.val) else out <- p.val
  }
  if("FC" %in% stats){
    fc <- getStat("FC", ncontrasts, num, fit, contrast, index)
    if(exists("out")) out <- c(out, fc) else out <- fc
  }
  
  ord <- switch(order.by, pval = order(p.val[[1]]),
                tstat = order(t.stat[[1]], decreasing = TRUE),
                fc = order(fc[[1]], decreasing = TRUE))
  out <- lapply(out, function(x) x[ord])
  out <- lapply(out, function(x)
                ifelse(abs(x) < 1e-3, sprintf("%0.3e", x), sprintf("%0.3f", x)))
  return(list(out = out, ord = ord))
}




#' Select and Output Genelists Based on Venn Diagrams
#' 
#' This function is designed to output text and/or HTML tables based on the
#' results of a call to \code{\link[limma]{decideTests}}.
#' 
#' The purpose of this function is to output HTML and text tables with lists of
#' genes that fulfill the criteria of a call to
#' \code{\link[limma]{decideTests}} as well as the direction of differential
#' expression.
#' 
#' Some important things to note: First, the names of the HTML and text tables
#' are extracted from the \code{colnames} of the \code{TestResults} object,
#' which come from the contrasts matrix, so it is important to use something
#' descriptive. Second, the method argument is analogous to the \code{include}
#' argument from \code{\link[limma:venn]{vennCounts}} or
#' \code{\link[limma:venn]{vennDiagram}}. Choosing "both" will select genes
#' that are differentially expressed in one or more comparisons, regardless of
#' direction. Choosing "up" or "down" will select genes that are only
#' differentially expressed in one direction. Choosing "same" will select genes
#' that are differentially expressed in the same direction. Choosing "sameup"
#' or "samedown" will select genes that are differentially expressed in the
#' same direction as well as 'up' or 'down'.
#' 
#' Note that this is different than sequentially choosing "up" and then "down".
#' For instance, a gene that is upregulated in one comparison and downregulated
#' in another comparison will be listed in the intersection of those two
#' comparisons if "both" is chosen, it will be listed in only one comparison
#' for both the "up" and "down" methods, and it will be listed in the union
#' (e.g., not selected) if "same" is chosen.
#' 
#' Calling the function normally will result in the output of HTML and text
#' tables:
#' 
#' vennSelect(eset, fit, design, x)
#' 
#' Calling the function with save set to \code{TRUE} will output both HTML and
#' text tables as well as a vector of counts for each comparison. This is
#' useful when using the function programmatically (e.g., when making reports
#' using Sweave).
#' 
#' out <- vennSelect(eset, fit, design, x, save = TRUE)
#' 
#' An alternative would be to use \code{vennCounts2} and
#' \code{\link[limma:venn]{vennDiagram}} to output a Venn diagram, which is
#' probably more reasonable since the tables being output are supposed to be
#' based on a Venn diagram.
#' 
#' @param eset An \code{ExpressionSet} object.
#' @param design A design matrix.
#' @param x A \code{\link[limma]{TestResults}} object, usually from a call to
#' \code{\link[limma]{decideTests}}.
#' @param contrast A contrasts matrix, produced either by hand, or by a call to
#' \code{\link[limma]{makeContrasts}}
#' @param fit An \code{\link[limma:marraylm]{MArrayLM}} object, from a call to
#' \code{\link[limma:ebayes]{eBayes}}.
#' @param method One of "same", "both", "up", "down", "sameup", or "samedown".
#' See details for more information.
#' @param adj.meth Method to use for adjusting p-values. Default is 'BH', which
#' corresponds to 'fdr'. Ideally one would set this value to be the same as was
#' used for \code{\link[limma]{decideTests}}.
#' @param stat The statistic to report in the resulting HTML tables. Choices
#' are 'fstat', 'tstat', and \code{NULL}. Ideally, the statistic chosen would
#' correspond to the method used in \code{\link[limma]{decideTests}}. In other
#' words, if one used methods such as 'separate' or 'hierarchical', which are
#' based on a t-statistic, one should choose 'tstat', however, if one used
#' 'nestedF', the logical choice would be 'fstat'.
#' @param otherstats Other statistics to be included in the HTML tables.
#' Choices include 'pval' and 'FC'.
#' @param order.by Which statistic should be used to order the probesets?
#' Choices include 'fstat', 'tstat', 'pval', and 'FC'. Note that if 'FC' is
#' chosen and there are more than one set of fold changes, the first will be
#' used.
#' @param foldFilt A log fold change to filter results.
#' @param save Boolean. Save the results for further processing?
#' @param titleadd Additional text to add to the title of the HTML tables.
#' Default is NULL, in which case the title of the table will be the same as
#' the filename.
#' @param ... Used to pass other arguments to \code{probes2table}, in
#' particular, to change the argument to \code{anncols} which controls the
#' columns of hyperlinks to online databases (e.g., Entrez Gene, etc.). See
#' \code{\link[annaffy]{aaf.handler}} for more information.
#' @return Normally called only for the side effect of producing HTML and text
#' tables. However, setting save to \code{TRUE} will output a vector of counts
#' that can be used for making Sweave-style reports.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @keywords manip
#' @export vennSelect
vennSelect <- function(eset, design, x, contrast, fit, method = "same", adj.meth = "BH",
                       stat = "fstat", otherstats = c("pval", "FC"), order.by = "pval",
                       foldFilt = NULL, save = FALSE, titleadd = NULL, ...){
  ## eset is an ExpressionSet containing data used for comparisons
  ## design is a design matrix from limma
  ## x is a TestResults object from a call to decideTests()
  ## output is a list containing the probe IDs of genes from each comparison

  if(!stat %in% c("fstat", "tstat")) stop(paste(stat, " is not an acceptable argument for the ",
                                               "'stat' argument.\nPlease use either ",
                                               "'fstat' or 'tstat'.", sep = ""), call. = FALSE)

  if(!is.null(foldFilt)){
    idx <- abs(fit$coefficients) > foldFilt
    x <- as.matrix(x) * idx
  }
  ncontrasts <- ncol(x)
  if(ncontrasts < 2 || ncontrasts > 3)
    stop("This function only works for two or three comparisons at a time.\n",
         call. = FALSE)
  if(ncontrasts == 2)
    name <- c(paste("Genes unique to", colnames(x)),
              paste("Genes in intersection of", intNames(x)))
  if(ncontrasts == 3)
    name <- c(paste("Genes unique to", colnames(x)),
              paste("Genes in intersection of", intNames(x)),
              "Genes common to all comparisons")
  
   ## Remove illegal characters from filenames
  if(length(grep("[/|\\|?|*|:|<|>|\"|\\|]", name)) > 0)
    warning(paste("Some illegal characters have been removed from the filenames",
                  name, sep = " "), call. = FALSE)
  name <- gsub("[/|\\|?|*|:|<|>|\"|\\|]", "", name)
  
  indices <- makeIndices(x, method = method)
  cols <- getCols(design, contrast)
  for(i in seq(along = indices)){
    tmp <- featureNames(eset)[indices[[i]]]
    if(stat == "fstat")
      otherdata <- f.stat.dat(fit, indices[[i]], contrast, adj.meth, c(stat, otherstats),
                              order.by, ncontrasts, i)
    if(stat == "tstat")
      otherdata <- t.stat.dat(fit, indices[[i]], contrast, adj.meth, c(stat, otherstats),
                              order.by, ncontrasts, i)
    if(is.null(stat)) otherdata <- NULL
    if(!is.null(titleadd))
        title <- paste(name[i], titleadd)
    else
        title <- name[i]
    if(length(tmp) == 0) next
    else
      probes2table(eset[,cols[[i]]], tmp[otherdata$ord], annotation(eset), text = TRUE,
                   filename = name[i], otherdata = otherdata$out, title = title, ...)
  }
  if(save)
    sapply(indices, sum)
}


