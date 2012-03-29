################################################
##
##  Copyright 2005 James W. MacDonald
##
##  Functions to both create Venn diagrams
##  and to output the subsets from the Venn
##  diagrams in HTML and text tables
##
################################################


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
      

getCols <- function(design, contrast){
  ## A function to get the correct columns of the ExpressionSet
  ## to output based on the design and contrasts matrices
  ## This is sort of a kludge and could probably be improved
  ncontrasts <- ncol(contrast)
  ids <- 1:dim(design)[1]
  tmp <- function(design, contrast, x){
    a <- design[,contrast[,x] > 0]
    if(is.matrix(a))
      a <- apply(a, 1, any)
    b <- design[,contrast[,x] < 0]
    if(is.matrix(b))
      b <- apply(b, 1, any)
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

