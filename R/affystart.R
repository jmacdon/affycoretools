###########################################################
##
##  Copyright 2005 James W. MacDonald
##
## affystart - functions to go from celfiles to
##             ExpressionSet with QC plots
##
########################################################





affystart <- function(..., filenames = NULL, groups=NULL, groupnames=NULL,
                      plot=TRUE, pca=TRUE, squarepca = FALSE, plottype="pdf",
                      express=c("rma", "mas5", "gcrma"), addname=NULL,
                      output = "txt", annotate = FALSE,
                      ann.vec = c("SYMBOL","GENENAME","ENTREZID","UNIGENE","REFSEQ")){

  if(is.null(filenames)){
    filenames <- list.files(pattern="\\.[cC][eE][lL]$")
    if(length(filenames) == 0)
      stop("There are no celfiles in the current working directory!", call. = FALSE)
  }
  dat <- ReadAffy(filenames = filenames )
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

    eset <- gcrma(dat)
  }
  if(express == "rma" || express == "gcrma"){
    if(annotate){
        ann <- addAnnot(eset, ann.vec = ann.vec)
        out <- data.frame(ann, exprs(eset))
        outPut(out, addname = addname, meth = output)
    }else{
        outPut(eset, addname = addname, meth = output)
    }

}
  if(express == "mas5"){
    out <- round(exprs(eset), 2)
    calls1 <- exprs(calls)
    calls2 <- round(assayDataElement(calls, "se.exprs"), 2)
    out.dat <- data.frame(cbind(out[,1],calls1[,1], calls2[,1]))
    if(dim(out)[2] > 1){
      for (i in 2:dim(out)[2]){
        out.dat <- data.frame(cbind(out.dat, out[,i], calls1[,i], calls2[,i]))
      }
    }
    nams <- NULL
    for(i in seq(along=filenames)) nams <- c(nams, filenames[i], "Call", "p-value")
    colnames(out.dat) <- nams
    if(annotate){
        ann <- addAnnot(eset, ann.vec = ann.vec)
        out <- data.frame(ann, out.dat)
        outPut(out, addname = addname, meth = output)
    }else{
        outPut(out.dat, addname = addname, meth = output)
    }
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

addAnnot <- function(eset, ann.vec = c("SYMBOL","GENENAME","ENTREZID","UNIGENE","REFSEQ")){
    require(annotation(eset), character.only = TRUE, quietly = TRUE) ||
    stop(paste("The", annotation(eset), "package needs to be installed to annotate your data"),
         call. = FALSE)
    make.name <- function(x, y)
        get(paste(x, y, sep = ""))
    make.vec <- function(chip, annot, probes){
        tmp <- mget(probes, make.name(chip, annot))
    if(any(sapply(tmp, length) > 1))
        sapply(tmp, function(x) paste(x, collapse = "//"))
    else
        unlist(tmp)
    }
    out <- data.frame(matrix(NA, ncol = length(ann.vec), nrow = length(featureNames(eset))),
                      stringsAsFactors = FALSE)
    for(i in seq(along = ann.vec))
        out[,i] <- make.vec(annotation(eset), ann.vec[i], featureNames(eset))
    names(out) <- ann.vec
    row.names(out) <- featureNames(eset)
    out
}


outPut <- function(out, addname, meth = c("txt","xls")){
    meth <- match.arg(meth, c("txt", "xls"))

    switch(meth,
           txt = {
               if(!is.null(addname)) nam <- paste("Expression values ",
                                                  addname, ".txt", sep = "")
               else nam <- "Expression values.txt"
               if(is(out, "data.frame"))
                   write.table(out, nam, quote = FALSE, sep = "\t",
                               col.names = NA)
               if(is(out, "ExpressionSet"))
                   write.exprs(out, nam, col.names = NA)
           }, ##xls = {
              ## require("xlsReadWrite", quietly = TRUE, character.only = TRUE) ||
              ## stop("The xlsReadWrite package is required to output Excel files", call.=FALSE)
              ## if(!is.null(addname)) nam <- paste("Expression values ",
              ##                                    addname, ".xls", sep = "")
              ## else nam <- "Expression values.xls"
              ## if(is(out, "data.frame"))
              ##     write.xls(out, nam)
              ## if(is(out, "ExpressionSet"))
              ##     write.xls(exprs(out), nam)
          ## }
           )
}


plotHist <- function(dat, filenames = NULL)
{
  if(is.null(filenames)) filenames <- sampleNames(dat)
  cl <- make.cl(filenames)
  if(length(filenames) <= 8){
    if(is(dat, "AffyBatch") || is(dat, "GeneFeatureSet"))
      hist(dat, lty=1, lwd=2, col=cl)
    if(is(dat, "matrix"))
      plotDensity(log2(dat), lty = 1, lwd = 2, col = cl)
    x.ax <- legend(1,1,legend=filenames, lty=1, lwd=2, col=cl, plot=FALSE)$rect$w
    legend(par("usr")[2] - (par("usr")[2] - par("usr")[1])/100 - x.ax,
           par("usr")[4] - (par("usr")[4]-par("usr")[3])/100,
           legend=filenames, lty=1, lwd=2, col=cl)
  }else{
    if(is(dat, "AffyBatch") || is(dat, "GeneFeatureSet"))
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
  layout(matrix(1:2, ncol = 2), c(3,1))
  plotAffyRNAdeg(AffyRNAdeg(dat), cols=1:length(filenames))

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

plotPCA <- function(object, groups = NULL, groupnames = NULL, addtext = NULL, x.coord = NULL, y.coord = NULL,
                    screeplot = FALSE, squarepca = FALSE, pch = NULL, col = NULL, pcs = c(1,2),
                    legend = TRUE, main = "Principal Components Plot", plot3d = FALSE, outside = FALSE, ...){
  if(is.character(groups)) stop("The groups argument should be numeric, not character!\n", call. = FALSE)  
  if(length(pcs) != 2 && !plot3d) stop("You can only plot two principal components.\n", call. = FALSE)
  if(length(pcs) != 3 && plot3d) stop("For 3D plotting, you should specify 3 principal components.\n", call. = FALSE)



  if(is(object, "ExpressionSet")){
      if(max(pcs) > dim(exprs(object))[2])
          stop(paste("There are only", dim(exprs(object))[2], "principal components to plot.\n", call. = FALSE))
      if(is.null(groupnames)) groupnames <- sampleNames(object)
      if(is.factor(groupnames)) groupnames <- as.character(groupnames)
      pca <- prcomp(t(exprs(object)))
      len <- length(sampleNames(object))
 }else{
     if(class(object) == "matrix"){
         if(max(pcs) > dim(object)[2])
             stop(paste("There are only", dim(object)[2], "principal components to plot.\n", call. = FALSE))
         if(is.null(groupnames)) groupnames <- colnames(object)
         if(is.factor(groupnames)) groupnames <- as.character(groupnames)
         pca <- prcomp(t(object))
         len <- dim(object)[2]
  }else{
      if(class(object) == "prcomp"){
          if(max(pcs) > dim(object$x)[2])
              stop(paste("There are only", dim(object$x)[2], "principal components to plot.\n", call. = FALSE))
          if(is.null(groupnames)) groupnames <- row.names(object$x)
          if(is.factor(groupnames)) groupnames <- as.character(groupnames)
          pca <- object
          len <- dim(object$x)[2]
      }else{
          stop("plotPCA currently only supports ExpressionSets, matrices and prcomp objects")
      }
  }
 }
  if(screeplot){
      plot(pca, main = "Screeplot")
  }else{
      if(plot3d){
          require("rgl", quietly = TRUE) || stop("The rgl package must be installed to do 3D plots.\n", call. = FALSE)
          plotstuff <- pcaPCH(len, groups, pch, col)
          plot3d(pca$x[,pcs], type = "s", col = plotstuff$col, size = 2)
          cat(paste("Sometimes rgl doesn't plot the first time.\nIf there",
                    "isn't anything in the plotting window, close it and re-run plotPCA().\n"))
      }else{
          if(legend && outside){
              opar <- par(no.readonly = TRUE)
              par(xpd = TRUE, mar = par()$mar + c(0,0,0,5))
          }
          if(squarepca){
              ylim <- max(abs(range(pca$x[,pcs[1]])))
              ylim <- c(-ylim, ylim)
          }else ylim <- NULL
          plotstuff <- pcaPCH(len, groups, pch, col)
          plot(pca$x[,pcs], pch = plotstuff$pch, col = plotstuff$col, bg = plotstuff$col,
               ylab= paste("PC", pcs[2], sep=""),
               xlab=paste("PC", pcs[1], sep=""),
               main = main, ylim = ylim, ...)

          if(!is.null(addtext)){
              smidge <-  (par("usr")[4] - par("usr")[3])/50
              text(pca$x[,pcs[1]], pca$x[,pcs[2]] + smidge, label = addtext,
                   cex = 0.7)
          }
          if(legend){
              if(outside){
                  if(is.null(groups))
                      unq <- unique(plotstuff)
                  else
                      unq <- unique(plotstuff[order(groups),])
                  legend(par("usr")[2], par("usr")[4], groupnames, pch = unq[,1],
                         col = unq[,2], pt.bg = unq[,2], bty = "n")
                  par(opar)
              }else{
                  pca.legend(pca, groups, groupnames, plotstuff, x.coord = x.coord,
                             y.coord = y.coord, ...)
              }
          }
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

pca.legend <- function(pca, groups, groupnames, pch.df,  x.coord = NULL, y.coord = NULL,
                       saveup = FALSE){
  ## A function to try to automagically place legend in a pca plot
  if(is.null(groups))
       unq <- unique(pch.df)
  else
      unq <- unique(pch.df[order(groups),])
  pch <- unq[,1]
  col <- unq[,2]
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

  whereto <- match(TRUE, c(upright, upleft, downleft, downright))
   if(!is.null(x.coord) & !is.null(y.coord)){
    legend(x.coord, y.coord, legend = groupnames, pch = pch, col = col, pt.bg = col, bty = "n")
  }else if(!is.na(whereto)){
    if(whereto == 1)
      legend(right, up + y.lab, legend=groupnames, pch=pch, col = col, pt.bg = col, bty = "n")
    if(whereto == 2)
      legend(left - x.lab, up + y.lab, legend=groupnames, pch=pch, col = col, pt.bg = col, bty = "n")
    if(whereto == 3)
      legend(left - x.lab, down, legend=groupnames, pch=pch, col = col, pt.bg = col, bty = "n")
    if(whereto == 4)
      legend(right, down, legend=groupnames, pch=pch, col = col, pt.bg = col, bty = "n")
  }else{
    answer <- readline("Please give the x-coordinate for a legend.")
    x.c <- as.numeric(answer)
    answer <- readline("Please give the y-coordinate for a legend.")
    y.c <- as.numeric(answer)
    legend(x.c, y.c, legend=groupnames, pch=pch, col = col, pt.bg = col, bty = "n")
  }
  if(saveup)
    return((par("usr")[4] - par("usr")[3])/50)
}

pcaPCH <- function(len, groups, pch, col){
    ## Function to minimize glyphs/colors used
    nulls <- is.null(c(pch, col)) || all(!is.null(pch), !is.null(col))
    if(nulls){
        if(!is.null(pch)){
            if(is.null(groups))
                out <- data.frame(pch, col, stringsAsFactors = FALSE)
            else
                out <- data.frame(pch[groups], col[groups], stringsAsFactors = FALSE)
        }else{
            if(is.null(groups)){
                pch <- rep(21:25, times = ceiling(len/5))[1:len]
                col <- colvec(len)[1:len]
                out <- data.frame(pch, col, stringsAsFactors = FALSE)
            }else{
                lg <- length(unique(groups))
                pch <- rep(21:25, times = ceiling(lg/5))[1:lg]
                col <- colvec(lg)[1:lg]
                out <- data.frame(pch, col, stringsAsFactors = FALSE)[groups,]
            }
        }
    }else{
        stop("You must either specify both pch and col arguments or neither\n",
             call. = FALSE)
    }
    out
}

colvec <- function(len){
    if(len > 64)
        tmp <- ceiling(sqrt(len))
    else
        tmp <- len
    mat <- matrix(1:tmp, ncol = tmp, nrow = tmp)
    tfs <- rbind(lower.tri(mat, TRUE), upper.tri(mat))
    mat <- rbind(mat,mat)
    out <- rainbow(tmp)[mat[tfs]]
    out
}

writeFit <- function(fit, annotation = NULL, eset){
    gt <- function(x) sapply(AnnotationDbi::mget(row.names(fit$t),
                                                 get(paste(annotation, x, sep = "")), ifnotfound = NA),
                             function(x) if(length(x) <= 1) x else paste(x, collapse = ";"))
    if(!is.null(annotation)){
        out <- data.frame(probeset = row.names(fit$t),
                          symbol = gt("SYMBOL"), 
                          description = gt("GENENAME"),
                          genbank = gt("ACCNUM"),
                          gene = gt("ENTREZID"),
                          unigene = gt("UNIGENE"))
    } else {
        out <- data.frame(probeset = row.names(fit$t))
    }
    nc <- ncol(fit$t)
    out2 <- lapply(1:nc, function(x) {
        z <- cbind(fit$coefficients[,x], fit$t[,x], fit$p.value[,x],
                   p.adjust(fit$p.value[,x], "BH"))
        colnames(z) <- paste(colnames(fit$t)[x], c("log fold change","t-stat","p-value","adj p-value"))
        z
    })
    out2 <- do.call("cbind", out2)
    out3 <- if(is(eset, "ExpressionSet")) exprs(eset) else eset
    return(cbind(out, out2, out3))
}


getMainProbes <- function(pdinfo){
    require(pdinfo, character.only = TRUE)
    con <- db(get(pdinfo))
    types <- dbGetQuery(con, paste("select distinct meta_fsetid, type from featureSet inner join core_mps",
                                   "on featureSet.fsetid=core_mps.fsetid;"))
    ind <- types$type %in% 1
    dbDisconnect(con)
    ind
}
  
