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

foldFilt <-function(object, fold = 1, groups, comps, compnames,
                    save = FALSE, text = TRUE, html = TRUE, filterfun = NULL){
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
        
getIndex <- function(comps, groups){
  first <- which(groups == comps[1])
  second <- which(groups == comps[2])
  c(first, second)
}

vennSelectFC <- function(eset, x, comps, order.by = "sum", method = "same", text = TRUE,
                         html = TRUE, ...){
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
