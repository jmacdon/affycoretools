###########################################################
##
##  Copyright 2005 James W. MacDonald
##
## affystart - functions to go from celfiles to exprSet
##             with QC plots
##
########################################################





affystart <- function(..., filenames = NULL, groups=NULL, groupnames=NULL,
                      plot=TRUE, pca=TRUE, squarepca = FALSE, plottype="pdf",
                      express=c("rma", "mas5", "gcrma"), addname=NULL,
                      phenoData = new("phenoData")){
  require(affy, quietly=TRUE)
  if(is.null(filenames)){
    filenames <- list.files(pattern="\\.[cC][eE][lL]$")
    if(length(filenames) == 0)
      stop("There are no celfiles in the current working directory!", call. = FALSE)
  }
  dat <- read.affybatch(filenames = filenames, phenoData = phenoData)
  filenames <- sub("\\.[cC][eE][lL]$", "", filenames)
  
  express <- match.arg(express, c("rma","mas5","gcrma"))
  if(express=="rma"){
    eset <- rma(dat)
  }
  if(express=="mas5"){
    eset <- mas5(dat)
    calls <- mas5calls(dat)
  }
  if(express=="gcrma"){
    require(gcrma, quietly=TRUE)
    eset <- gcrma(dat)
  }
  if(express == "rma" |express == "gcrma"){
    if(is.null(addname))
      write.exprs(eset, "Expression values.txt", col.names=NA)
    else
      write.exprs(eset, paste("Expression values", paste(" ", addname), ".txt",  sep=""), col.names=NA)
  }
  if(express == "mas5"){
    out <- round(exprs(eset), 2)
    calls1 <- exprs(calls)
    calls2 <- round(se.exprs(calls), 2)
    out.dat <- data.frame(cbind(out[,1],calls1[,1], calls2[,1]))
    if(dim(out)[2] > 1){
      for (i in 2:dim(out)[2]){
        out.dat <- data.frame(cbind(out.dat, out[,i], calls1[,i], calls2[,i]))
      }
    }
    nams <- NULL
    for(i in seq(along=filenames)) nams <- c(nams, filenames[i], "Call", "p-value")
    colnames(out.dat) <- nams
    if(is.null(addname))
       write.table(out.dat, "Expression values.txt", sep="\t", quote=FALSE, col.names=NA)
    else
      write.table(out.dat, paste("Expression values", paste(" ", addname), ".txt",  sep=""), sep="\t", quote=FALSE, col.names=NA)
  }

   ##Plots
  
  if(plot){
    plotHist(dat, filenames)
    saved.hist <- recordPlot()
    plotDeg(dat, filenames)
    saved.deg <- recordPlot()
    if(pca){
      plotPCA(eset, groups, groupnames, squarepca = squarepca)
      saved.pca <- recordPlot()
    }
    meth <- c("pdf", "postscript", "jpeg", "bmp", "png")
    if(plottype %in% meth){
      method <- get(plottype)
      ## set height and width so plots look reasonable
      if(plottype %in% meth[1:2]) height <- width <- 7.5
      if(plottype %in% meth[3:5]) height <- width <- 700
      ## Get correct extension
      if(plottype == "postscript") plottype <- "ps"
      if(plottype == "jpeg") plottype <- "jpg"
      method(paste("Density plot", plottype, sep = "."),
             height = height, width = width, ...)
      replayPlot(saved.hist)
      dev.off()
      method(paste("Digestion plot", plottype, sep = "."),
             height = height, width = width, ...)
      replayPlot(saved.deg)
      dev.off()
      if(pca){
        method(paste("PCA plot", plottype, sep = "."),
               height = height, width = width, ...)
        replayPlot(saved.pca)
        dev.off()
      }
    }else stop("Currently this function only supports outputting\n",
               "pdf, postscript, jpeg, bmp and png files.\n", call. = FALSE)
  }
  return(eset)
}

plotHist <- function(dat, filenames = NULL)
{
  if(is.null(filenames)) filenames <- sampleNames(dat)
  cl <- make.cl(filenames)
  if(length(filenames) <= 8){
    if(is(dat, "AffyBatch"))
      hist(dat, lty=1, lwd=2, col=cl)
    if(is(dat, "matrix"))
      plotDensity(log2(dat), lty = 1, lwd = 2, col = cl)
    x.ax <- legend(1,1,legend=filenames, lty=1, lwd=2, col=cl, plot=FALSE)$rect$w
    legend(par("usr")[2] - (par("usr")[2] - par("usr")[1])/100 - x.ax,
           par("usr")[4] - (par("usr")[4]-par("usr")[3])/100,
           legend=filenames, lty=1, lwd=2, col=cl)
  }else{
    if(is(dat, "AffyBatch"))
      hist(dat, lty=cl, lwd=2, col=1:length(filenames))
    if(is(dat, "matrix"))
      plotDensity(log2(dat), lty = cl, lwd = 2, col = 1:length(filenames))
    x.ax <- legend(1,1,legend=filenames, lty=1, lwd=2, col=cl, plot=FALSE)$rect$w
    y.ax <- legend(1,1,legend=filenames, lty=1, lwd=2, col=cl, plot=FALSE)$rect$h
    ydiff <- par("usr")[4] - par("usr")[3]
    ## If legend is too big, shrink to fit
    if(y.ax < ydiff){
      legend(par("usr")[2] - (par("usr")[2] - par("usr")[1])/100 - x.ax,
             par("usr")[4] - (par("usr")[4]-par("usr")[3])/100,
             legend=filenames, lty=cl, lwd=2, col=1:length(filenames))
    }else{
      cexval <- 1
      while(y.ax > ydiff){
        cexval <- cexval - 0.05
        y.ax <- legend(1,1, legend=filenames, lty=1, lwd=2, col=cl, plot=FALSE, cex=cexval)$rect$h
        x.ax <- legend(1,1,legend=filenames, lty=1, lwd=2, col=cl, plot=FALSE, cex=cexval)$rect$w
      }
      legend(par("usr")[2] - (par("usr")[2] - par("usr")[1])/100 - x.ax,
             par("usr")[4] - (par("usr")[4]-par("usr")[3])/100,
             legend=filenames, lty=cl, lwd=2, col = 1:length(filenames), cex=cexval)
    }
    
  }
}

plotDeg <- function(dat, filenames = NULL){
  if(is.null(filenames)) filenames <- sampleNames(dat)

  ## reset things when exiting
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  ## put plot on left, legend on right
  layout(matrix(1:2, nc = 2), c(3,1))
  plotAffyRNAdeg(AffyRNAdeg(dat), col=1:length(filenames))

  ## fake a plot
  par(mai = c(0,0,1.01,0))
  plot(1:10, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
  
  tmp <- legend("topleft", legend=filenames, lty=1, lwd=2,
                col=1:length(filenames), plot=FALSE)

  y.ax <- tmp$rect$h
  x.ax <- tmp$rect$w
  ydiff <- par("usr")[4] - par("usr")[3]
  xdiff <- par("usr")[2] - par("usr")[1]

  
  ## If legend is too big, shrink to fit
  if(y.ax < ydiff && x.ax < xdiff){
    legend("topleft", legend=filenames, lty=1, lwd=2, col=1:length(filenames))
  }else{
    cexval <- 1
    while(y.ax > ydiff){
      cexval <- cexval - 0.05
      tmp <- legend("topleft", legend=filenames, lty=1, lwd=2,
                    col=1:length(filenames), plot=FALSE, cex=cexval)
      y.ax <- tmp$rect$h
      x.ax <- tmp$rect$w
    }
    if(x.ax < xdiff){
      legend("topleft", legend=filenames, lty=1, lwd=2, col=1:length(filenames), cex=cexval)
    }else{
      while(x.ax > xdiff){
        cexval <- cexval - 0.05
        x.ax <- legend("topleft", legend=filenames, lty=1, lwd=2,
                       col=1:length(filenames), plot=FALSE, cex=cexval)$rect$w
      }
      legend("topleft", legend=filenames, lty=1, lwd=2, col=1:length(filenames), cex=cexval)
    }
  }
}

plotPCA <- function(eset, groups = NULL, groupnames = NULL, addtext = NULL, x.coord = NULL, y.coord = NULL,
                    screeplot = FALSE, squarepca = FALSE, pch = NULL, col = NULL, ...){
  if(is.null(groupnames)) groupnames <- sampleNames(eset)
  if(is.factor(groupnames)) groupnames <- as.character(groupnames)
  pca <- prcomp(t(exprs(eset)))
  if(screeplot){
    plot(pca, main = "Screeplot")
  }else{
    if(squarepca){
      ylim <- max(abs(range(pca$x[,1])))
      ylim <- c(-ylim, ylim)
    }else ylim <- NULL
    if(!is.null(groups)){
      if(is.null(pch)) pch <- groups
      if(is.null(col)) col <- groups
      plot(pca$x[,1:2], pch = pch, col = col, ylab="PC2", xlab="PC1",
           main="Principal Components Plot", ylim = ylim, ...)
    }else{
      if(is.null(pch)) pch <- 0:length(sampleNames(eset))
      if(is.null(col)) col <- 1:length(sampleNames(eset))
      plot(pca$x[,1:2], pch=pch, col=col,
           ylab="PC2", xlab="PC1", main="Principal Components Plot", ylim = ylim, ...)
    }
    if(is.null(addtext)){
      pca.legend(pca, groupnames, pch, col, x.coord = x.coord, y.coord = y.coord, ...)
    }else{
      smidge <-  pca.legend(pca, groupnames, pch, col, x.coord = x.coord, y.coord = y.coord,
                            saveup = TRUE, ...)
      text(pca$x[,1], pca$x[,2] + smidge, label = addtext, cex = 0.7)
    }
  }
}


make.cl <- function(filenames){
  ## A function to make a classlabel for plotting
  ## Check for number of files
  if(length(filenames) <= 8)
    cl <- 1 : length(filenames)
  if(length(filenames) > 8){
    mod <- floor(length(filenames)/8)
    cl <- NULL
    for(i in 1 : mod){
      cl <- c(cl, rep(i, 8))
    }
    rem <- length(filenames) - mod*8
    cl <- c(cl, rep(mod + 1, rem))
  }
  cl
}

pca.legend <- function(pca, groupnames, pch, col, x.coord = NULL, y.coord = NULL,
                       saveup = FALSE){
  ## A function to try to automagically place legend in a pca plot

  pch <- sort(unique(pch))
  col <- sort(unique(col))
  x.lab <- legend(1, 1, legend = groupnames, pch = pch, plot = FALSE)$rect$w
  y.lab <- legend(1, 1, legend = groupnames, pch = pch, plot = FALSE)$rect$h
  

  right <- par("usr")[2] - (par("usr")[2] - par("usr")[1])/100 - x.lab
      left <- par("usr")[1] + (par("usr")[2] - par("usr")[1])/100 + x.lab
  up <- par("usr")[4] - (par("usr")[4] - par("usr")[3])/100 - y.lab
  down <- par("usr")[3] + (par("usr")[4] - par("usr")[3])/100 + y.lab
  
  upright <- !any(pca$x[,1] > right & pca$x[,2] > up)
  upleft <- !any(pca$x[,1] < left & pca$x[,2] > up)
  downright <- !any(pca$x[,1] > right & pca$x[,2] < down)
  downleft <- !any(pca$x[,1] < left & pca$x[,2] < down)
  
  where <- match(TRUE, c(upright, upleft, downleft, downright))
  if(!is.na(where)){
    if(where == 1)
      legend(right, up + y.lab, legend=groupnames, pch=pch, col=col)
    if(where == 2)
      legend(left - x.lab, up + y.lab, legend=groupnames, pch=pch, col=col)
    if(where == 3)
      legend(left - x.lab, down, legend=groupnames, pch=pch, col=col)
    if(where == 4)
      legend(right, down, legend=groupnames, pch=pch, col=col)
  }else if(!is.null(x.coord) & !is.null(y.coord)){
    legend(x.coord, y.coord, legend = groupnames, pch = pch, col = col)
  }else{
    answer <- readline("Please give the x-coordinate for a legend.")
    x.c <- as.numeric(answer)
    answer <- readline("Please give the y-coordinate for a legend.")
    y.c <- as.numeric(answer)
    legend(x.c, y.c, legend=groupnames, pch=pch, col=col)
  }
  if(saveup)
    return((par("usr")[4] - par("usr")[3])/50)
}

