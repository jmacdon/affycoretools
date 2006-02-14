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
  if(is(object, "exprSet"))
    x  <- exprs(object)
  if(length(unique(groups)) != length(groups)){
    gps <- matrix(NA, nc = length(unique(groups)), nr = dim(x)[1])
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
  FCs <- vector("list", unique(groups))
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
  direct <- matrix(NA, nc = length(comps), nr = dim(gps)[1])
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
