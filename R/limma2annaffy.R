################################################
##
##  limma2annaffy.R
##  Version 0.02
## Copyright 2005 James W. MacDonald
##
##
##  Modifications
##
##  5-27-05 Added tableFilt() and removed
##  extraneous fold change calculations
##  FC now extracted directly from topTable object
##
#################################################

limma2annaffy <- function(eset, fit, design, contrast, lib, adjust = "fdr",
                          anncols = aaf.handler()[c(1:3, 6:7, 9:12)], number = 30,
                          pfilt = NULL, fldfilt = NULL, tstat = TRUE, pval = TRUE, FC = TRUE,
                          expression = TRUE, html = TRUE, text = FALSE, save = FALSE,
                          addname = NULL, interactive = TRUE){
  ## if lib isn't a .db package, make it so
  if(length(grep("\\.db$", lib)) < 1)
      lib <- paste(lib, "db", sep = ".")

  if(!interactive){
    limma2annaffy.na(eset = eset, fit = fit, design = design, contrast = contrast,
                     lib = lib, adjust = adjust, anncols = anncols, number = number,
                     pfilt = pfilt, fldfilt = fldfilt, tstat = tstat, pval = pval,
                     FC = FC, expression = expression, html = html, text = text,
                     save = save, addname = addname)
  }else{



    tables <- vector("list", dim(contrast)[2])
    for(i in seq(along = colnames(contrast))){
      tables[[i]] <- tableFilt(fit, coef = i, number = number, pfilt = pfilt, fldfilt = fldfilt,
                               adjust = adjust)
    }

    ##Check to see if the table dimensions are ok for the end user
    ##This part needs some error handling...
    rn <- vector()
    for(i in 1:length(tables))rn[i] <- dim(tables[[i]])[1]
    cat("\nYou are going to output", length(tables),"tables,",
        "\nWith this many genes in each:\n", rn)
    cat("\n")
    doit <- getans2("Do you want to accept or change these values?", allowed=c("a","c"))
    if(doit == "c"){
      dopval <- getans2("Do you want to change the p-value?", allowed=c("y","n"))
      if(dopval == "y")
        pfilt <- getpfil()
      dofld <- getans2("Do you want to change the fold filter?", allowed = c("y", "n"))
      if(dofld == "y")
        fldfilt <- getfldfil()
      limma2annaffy(eset=eset, fit=fit, design=design, contrast=contrast,
                    lib=lib, anncols=anncols, pfilt=pfilt, fldfilt=fldfilt,
                    adjust=adjust, tstat=tstat, pval=pval, FC=FC, expression=expression,
                    html=html, text=text)
      options(show.error.messages = FALSE)
      stop()
    }
    ## Use topTables to make output tables
    ## Get filename from contrast matrix
    for(i in 1:length(tables)){
      if(dim(tables[[i]])[1]>0){
        index <- as.numeric(row.names(tables[[i]]))
        if(length(which(contrast[,i] > 0)) == 1){
          grp1 <- which(design[,which(contrast[,i] > 0)] > 0)
        }else{
          grp1 <- which(rowSums(design[,which(contrast[,i] > 0)]) > 0)
        }
        if(length(which(contrast[,i] < 0)) == 1){
          grp2 <- which(design[,which(contrast[,i] < 0)] > 0)
        }else{
          grp2 <- which(rowSums(design[,which(contrast[,i] < 0)]) > 0)
        }
        grps <- c(grp1, grp2)

        filename <- colnames(contrast)[i]
        if(!is.null(addname))
          filename <- paste(filename, addname, sep=" ")
        ## Remove illegal characters from filename
        if(length(grep("[/|\\|?|*|:|<|>|\"|\\|]", filename)) > 0)
          warning(paste("Some illegal characters have been removed from the filename",
                        filename, sep = " "), call. = FALSE)
        filename <- gsub("[/|\\|?|*|:|<|>|\"|\\|]", "", filename)

        ## Make aafTable object
        probeids <- featureNames(eset)[index]
        anntable <- aafTableAnn(probeids, lib, anncols)
        if(tstat)
          testtable <- aafTable("t-statistic" = round(tables[[i]][,"t"],2))
        if(pval){
          if(!exists("testtable")){
            testtable <- aafTable("p-value" = round(tables[[i]][,"adj.P.Val"],3))
          }else{
            testtable <- merge(testtable,
                               aafTable("p-value" = round(tables[[i]][,"adj.P.Val"],3)))
          }
        }
        if(FC){
          fld <- tables[[i]][,"logFC"]
          if(!exists("testtable")){
            testtable <- aafTable("Fold Change" = round(fld, 2))
          }else{
            testtable <- merge(testtable, aafTable("Fold Change" = round(fld, 2)))
          }
      }

        if(exists("testtable"))
          anntable <- merge(anntable, testtable)

        if(expression)
          exprtable <- aafTableInt(eset[,grps], probeids=probeids)

        if(exists("exprtable"))
          anntable <- merge(anntable, exprtable)

        if(html)
          saveHTML(anntable, paste(filename,"html", sep="."), filename)
        if(text)
          saveText(anntable, paste(filename, "txt", sep="."), header=TRUE)
      }
      if(exists("testtable")) rm(testtable)
    }
    if(save) return(tables)
  }
}

"limma2annaffy.na" <- function(eset, fit, design, contrast, lib, adjust = "fdr",
                               anncols = aaf.handler()[c(1:3, 6:7, 9:12)], number = 30,
                               pfilt = NULL, fldfilt = NULL, tstat = TRUE, pval = TRUE, FC = TRUE,
                               expression = TRUE, html = TRUE, text = FALSE, save = FALSE,
                               addname = NULL){



  tables <- vector("list", dim(contrast)[2])
  for(i in seq(along = colnames(contrast))){
    tables[[i]] <- tableFilt(fit, coef = i,  number = number, pfilt = pfilt, fldfilt = fldfilt,
                             adjust = adjust)
  }
  ## Use topTables to make output tables
  ## Get filename from contrast matrix
  for(i in 1:length(tables)){
    if(dim(tables[[i]])[1]>0){
      index <- as.numeric(row.names(tables[[i]]))
      if(length(which(contrast[,i] > 0)) == 1){
        grp1 <- which(design[,which(contrast[,i] > 0)] > 0)
      }else{
        grp1 <- which(rowSums(design[,which(contrast[,i] > 0)]) > 0)
      }
      if(length(which(contrast[,i] < 0)) == 1){
        grp2 <- which(design[,which(contrast[,i] < 0)] > 0)
      }else{
        grp2 <- which(rowSums(design[,which(contrast[,i] < 0)]) > 0)
      }
      grps <- c(grp1, grp2)

      filename <- colnames(contrast)[i]
      if(!is.null(addname))
        filename <- paste(filename, addname, sep=" ")
      ## Remove illegal characters from filename
      if(length(grep("[/|\\|?|*|:|<|>|\"|\\|]", filename)) > 0)
        warning(paste("Some illegal characters have been removed from the filename",
                      filename, sep = " "), call. = FALSE)
      filename <- gsub("[/|\\|?|*|:|<|>|\"|\\|]", "", filename)

      ## Make aafTable object
      probeids <- featureNames(eset)[index]
      anntable <- aafTableAnn(probeids, lib, anncols)
      if(tstat)
        testtable <- aafTable("t-statistic" = round(tables[[i]][,"t"],2))
      if(pval){
        if(!exists("testtable")){
            testtable <- aafTable("p-value" = round(tables[[i]][,"adj.P.Val"],3))
          }else{
            testtable <- merge(testtable,
                               aafTable("p-value" = round(tables[[i]][,"adj.P.Val"],3)))
          }
      }
      if(FC){
        fld <- tables[[i]][,"logFC"]
        if(!exists("testtable")){
          testtable <- aafTable("Fold Change" = round(fld, 2))
        }else{
          testtable <- merge(testtable, aafTable("Fold Change" = round(fld, 2)))
        }
      }

      if(exists("testtable"))
        anntable <- merge(anntable, testtable)

      if(expression)
        exprtable <- aafTableInt(eset[,grps], probeids=probeids)

      if(exists("exprtable"))
        anntable <- merge(anntable, exprtable)

      if(html)
        saveHTML(anntable, paste(filename,"html", sep="."), filename)
      if(text)
          saveText(anntable, paste(filename, "txt", sep="."), header=TRUE)
    }
    if(exists("testtable")) rm(testtable)
  }
  if(save) return(tables)
}

tableFilt <- function(fit, coef = 1,  number = 30, fldfilt = NULL, pfilt = NULL,
                      adjust = "fdr"){
  if(is.null(fldfilt) && is.null(pfilt)){
    tab <- topTable(fit, coef = coef, number = number, adjust = adjust)
  }else{
    tab <- topTable(fit, coef = coef, number = dim(fit$coefficients)[1],
                    adjust = adjust)
  }
## Filter on p-value
  if(!is.null(pfilt))
    tab <- switch(adjust,
                  none = tab[tab[,"P.Value"] < pfilt,],
                  tab[tab[,"adj.P.Val"] < pfilt,])
  ## Filter on fold change
  if(!is.null(fldfilt))
    tab <- tab[abs(tab[,"logFC"]) > fldfilt,]
  tab
}
