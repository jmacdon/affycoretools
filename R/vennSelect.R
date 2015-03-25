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

