################################################
##
##  Copyright 2005 James W. MacDonald
##
##  Functions to both create Venn diagrams
##  and to output the subsets from the Venn
##  diagrams in HTML and text tables
##
################################################


vennCounts2 <- function(fit, x, method = "same", foldFilt = 0
                       ){
  ## x is a TestResults object from a call to
  ## decideTests()
  ## output is the counts of genes where the requirement is
  ## that they all go the same direction
  require("limma", quietly = TRUE)
  tmp <- vennCounts(x)
  newcnts <- sapply(vennSelect(fit = fit, x = x, method = method,
                               foldFilt = foldFilt, indices.only = TRUE), sum)
  all <- sum(tmp[,dim(tmp)[2]])
  if(dim(x)[2] == 2)
    ord <- c(2, 1, 3)
  if(dim(x)[2] == 3)
    ord <- c(3, 2, 6, 1, 5, 4, 7)
  tmp[,dim(tmp)[2]] <- c(all - sum(newcnts), newcnts[ord])
  tmp
}

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

makeIndices <- function(x, method = "same"){
  method <- match.arg(method, c("same", "up", "down", "both", "sameup", "samedown"))
  indices <- list()
  x <- sign(switch(method, both = abs(x), up = x > 0, down = x < 0, same = x,
                   sameup = x, samedown = x))
  if(ncol(x) == 2){
    if(method != "sameup" && method != "samedown"){
      indices[[1]] <- abs(x[,1]) == 1 & x[,2] == 0
      indices[[2]] <- x[,1] == 0 & abs(x[,2]) == 1
      indices[[3]] <- abs(rowSums(x)) == 2
    }
    if(method == "sameup"){
      indices[[1]] <- x[,1] == 1 & x[,2] == 0
      indices[[2]] <- x[,1] == 0 & x[,2] == 1
      indices[[3]] <- rowSums(x) == 2
    }
    if(method == "samedown"){
      indices[[1]] <- x[,1] == -1 & x[,2] == 0
      indices[[2]] <- x[,1] == 0 & x[,2] == -1
      indices[[3]] <- rowSums(x) == -2
    }
  }
  if(ncol(x) == 3){
    if(method != "sameup" && method != "samedown"){
      indices[[1]] <- abs(x[,1]) == 1 & x[,2] == 0 & x[,3] == 0
      indices[[2]] <- x[,1] == 0 & abs(x[,2]) == 1 & x[,3] == 0
      indices[[3]] <- x[,1] == 0 & x[,2] == 0 & abs(x[,3]) == 1
      indices[[4]] <- abs(rowSums(x[,1:2])) == 2 & x[,3] == 0
      indices[[5]] <- abs(rowSums(x[,c(1,3)])) == 2 & x[,2] == 0
      indices[[6]] <- abs(rowSums(x[,2:3])) == 2 & x[,1] == 0
      indices[[7]] <- abs(rowSums(x)) == 3
    }
    if(method == "sameup"){
      indices[[1]] <- x[,1] == 1 & x[,2] == 0 & x[,3] == 0
      indices[[2]] <- x[,1] == 0 & x[,2] == 1 & x[,3] == 0
      indices[[3]] <- x[,1] == 0 & x[,2] == 0 & x[,3] == 1
      indices[[4]] <- rowSums(x[,1:2]) == 2 & x[,3] == 0
      indices[[5]] <- rowSums(x[,c(1,3)]) == 2 & x[,2] == 0
      indices[[6]] <- rowSums(x[,2:3]) == 2 & x[,1] == 0
      indices[[7]] <- rowSums(x) == 3
    }
    if(method == "samedown"){
      indices[[1]] <- x[,1] == -1 & x[,2] == 0 & x[,3] == 0
      indices[[2]] <- x[,1] == 0 & x[,2] == -1 & x[,3] == 0
      indices[[3]] <- x[,1] == 0 & x[,2] == 0 & x[,3] == -1
      indices[[4]] <- rowSums(x[,1:2]) == -2 & x[,3] == 0
      indices[[5]] <- rowSums(x[,c(1,3)]) == -2 & x[,2] == 0
      indices[[6]] <- rowSums(x[,2:3]) == -2 & x[,1] == 0
      indices[[7]] <- rowSums(x) == -3
    }
  }
  indices
}
      

getOrd <- function(x, design){
  ord <- vector()
  for(i in seq(along = colnames(design)))
    ord <- c(ord, which(design[,i] == 1))
  ord
}
    
  
vennSelect <- function(eset, fit, design, x, method = "same", foldFilt = 0,
                       indices.only = FALSE, save = FALSE, ...){
  ## eset is exprSet containing data used for comparisons
  ## design is a design matrix from limma
  ## x is a TestResults object from a call to decideTests()
  ## output is a list containing the probe IDs of genes from each comparison

  require("limma", quietly = TRUE)
  idx <- abs(fit$coefficients) > foldFilt
  x <- as.matrix(x) * idx
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
  indices <- makeIndices(x, method = method)
  if(indices.only)
    return(indices)
  ## ugly hack to get some ordering for genes
  if(ncontrasts == 2)
    cont.ind <- c(1,2,1)
  else
    cont.ind <- c(1,2,3,1,1,2,1)
  ord <- getOrd(x, design)
  for(i in seq(along = indices)){
    tmp <- geneNames(eset)[indices[[i]]]
    gn.ord <- order(abs(fit$coefficients[tmp, colnames(x)[cont.ind[i]]]), decreasing = TRUE)
    if(length(tmp) == 0) next
    else
      probes2table(eset[,ord], tmp[gn.ord], annotation(eset), text = TRUE,
                   filename = name[i], ...)
  }
  if(save)
    sapply(indices, sum)
}
    
  
  
