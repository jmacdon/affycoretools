###################################################
##
##  Copyright 2005 James W. MacDonald
##
##  MA plotting functions
##
####################################################

plotma <- function(x, y=NULL, main="", log=TRUE, out=FALSE){

## This function makes MA plots using either log2 or unlogged data. There is also the ability to output
## the results for e.g., adding a lowess line.
    .Deprecated(msg = paste("plotma is deprecated. See either plotMA in limma or mva.pairs in affy\n",
                "This function will be removed by the next release.\n", sep=""))

  if(is.matrix(x)){
    if(dim(x)[2] > 2) cat("\nThe data matrix passed to plot.ma has more than two columns.\n",
                                  "Only the first two will be used to plot.\n", sep="")
    y <- x[,2]
    x <- x[,1]
  }
  if(log){
    M <- x-y
    A <- (x+y)/2
  }else{
    M <- log2(x/y)
    A <- log2(sqrt((x/10)*(y/10)*100))
  }
  if(is.character(main)) main <- main
  ylim <- c(-max(abs(M)), max(abs(M)))
  plot(A, M, pch=".", cex=0.4, main=main, ylim=ylim, ylab="Log Ratio", xlab="Average Intensity")
  if(out) return(cbind(A,M))
}


convert.back <- function(M, A){ #Here M = maM values and A = maA values
  .Deprecated(msg = paste("convert.back is deprecated. I don't think anybody uses this function,\n",
              "so it will be removed in the next release.\n", sep=""))
  G <- (2*A-M)/2
  R <- (2*A+M)/2
  return(cbind(R,G))
}
