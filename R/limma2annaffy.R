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



#' Function to Create HTML Tables from limma Objects
#' 
#' This function is designed to take an \code{ExpressionSet} and an
#' \code{lmFit}, \code{model.matrix}, and contrast object from limma and
#' convert into HTML tables using annaffy. The alternate function
#' \code{limma2annaffy.na} is designed to be run without user intervention.
#' 
#' This function is designed to automatically output HTML or text tables, with
#' filenames taken from the column names of the contrast matrix. The number of
#' genes output can be controlled several different ways. First, if pfilt and
#' fldfilt are both \code{NULL}, the top genes will be output based on the
#' \code{number} variable. Otherwise, the genes are filtered based on p-value,
#' fold change, or both. If the genes are filtered this way, the number of
#' genes to be output will be listed and the filter(s) can then be adjusted if
#' necessary.
#' 
#' This function currently only supports Affymetrix data.
#' 
#' @aliases limma2annaffy limma2annaffy.na
#' @param eset An \code{ExpressionSet} containing affymetrix expression values.
#' @param fit An \code{lmFit} object.
#' @param design A \code{model.matrix} object.
#' @param contrast A contrasts matrix from limma.
#' @param lib An annotation package for the Affy chips used.
#' @param adjust Multiplicity adjustment. Choices are
#' "fdr","holm","hommel","bonferroni", or "none". Partial matching allowed.
#' @param anncols A vector of things to annotate, produced by a call to
#' aaf.handler().
#' @param number Number of genes to output to table. See details for more
#' information.
#' @param pfilt A p-value to filter output. See details for more information.
#' @param fldfilt A fold change to filter output. See details for more
#' information.
#' @param tstat Boolean: Output t-statistics in table? Defaults to
#' \code{FALSE}.
#' @param pval Boolean: Output (adjusted) p-values in table?  Defaults to
#' \code{FALSE}.
#' @param FC Boolean: Output fold changes in table? Defaults to \code{FALSE}.
#' @param expression Boolean: Output expression values in table? Defaults to
#' \code{TRUE}.
#' @param html Boolean: Output data in HTML tables? Defaults to \code{TRUE}.
#' @param text Boolean: Output data in text tables? Defaults to \code{TRUE}.
#' @param save Boolean: Save tables as R objects for further processing?
#' Defaults to \code{FALSE}.
#' @param addname A character vector to add to the end of the automatically
#' generated output file names. Useful for multiple calls to eliminate
#' over-writing of existing HTML or text tables.
#' @param addtitle A character vector to add to the title for the HTML table.
#' By default the title will be the same as the filename. If the addname
#' argument is not \code{NULL}, then that will be appended to the filename (and
#' will be used as the HTML title). If addtitle is not \code{NULL}, it will be
#' appended to the filename and that will then be used as the HTML table title.
#' @param interactive Boolean: Is this an interactive call, or run as part of a
#' script (e.g., in an \code{Sweave} document)? Defaults to \code{TRUE}
#' @param natFC Boolean: Add 'unlogged' fold change to output? This is intended
#' for people who don't understand logs or fractions. If the fold change is
#' positive, it is simply exponentiated (e.g., 2^x where x is the log fold
#' change). If negative, we use -2^(-x), so e.g., a log fold change of -2 will
#' result in a -4.
#' @return If \code{save} is \code{TRUE}, a list of tables from
#' \code{\link[limma:toptable]{topTable}} will be output.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @seealso \code{\link[limma:toptable]{topTable}},
#' \code{\link[annaffy]{aafTableAnn}}
#' @keywords manip
#' @export limma2annaffy
limma2annaffy <- function(eset, fit, design, contrast, lib, adjust = "fdr",
                          anncols = aaf.handler()[c(1:3, 6:7, 9:12)], number = 30,
                          pfilt = NULL, fldfilt = NULL, tstat = TRUE, pval = TRUE, FC = TRUE,
                          expression = TRUE, html = TRUE, text = FALSE, save = FALSE,
                          addname = NULL, addtitle = NULL, interactive = TRUE, natFC = FALSE){
  .Deprecated(new = "", msg = paste("limma2annaffy is being deprecated. Please see the RefactoredAffycoretools vignette",
                        "for more current ways to create annotated output"))
  ## if lib isn't a .db package, make it so
  if(length(grep("\\.db$", lib)) < 1)
      lib <- paste(lib, "db", sep = ".")

  if(!interactive){
    limma2annaffy.na(eset = eset, fit = fit, design = design, contrast = contrast,
                     lib = lib, adjust = adjust, anncols = anncols, number = number,
                     pfilt = pfilt, fldfilt = fldfilt, tstat = tstat, pval = pval,
                     FC = FC, expression = expression, html = html, text = text,
                     save = save, addname = addname, addtitle = addtitle,
                     natFC = natFC)
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
                    html=html, text=text, save=save, addname=addname, addtitle=addtitle,
                    natFC=natFC)
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
        if(!is.null(addtitle))
            title <- paste(filename, addtitle)
        else
            title <- filename
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
            testtable <- aafTable("p-value" = ifelse(tables[[i]][,"adj.P.Val"] < 1e-3,
                                  sprintf("%0.3e", tables[[i]][,"adj.P.Val"]),
                                  sprintf("%0.3f", tables[[i]][,"adj.P.Val"])))
          }else{
            testtable <- merge(testtable,
                               aafTable("p-value" = ifelse(tables[[i]][,"adj.P.Val"] < 1e-3,
                                  sprintf("%0.3e", tables[[i]][,"adj.P.Val"]),
                                  sprintf("%0.3f", tables[[i]][,"adj.P.Val"]))))
          }
        }
        if(FC){
          fld <- tables[[i]][,"logFC"]
          if(!exists("testtable")){
            testtable <- aafTable("log2 fold Change" = round(fld, 2))
          }else{
            testtable <- merge(testtable, aafTable("log2 fold Change" = round(fld, 2)))
          }
      }
      if(natFC){
        fld <- tables[[i]][,"logFC"]
        fld <- ifelse(fld > 0, round(2^fld, 2), -round(2^-fld, 2))
          if(!exists("testtable")){
            testtable <- aafTable("Fold change natural scale" = fld)
          }else{
            testtable <- merge(testtable, aafTable("Fold change natural scale" = fld))
          }
      }  
        
        if(exists("testtable"))
          anntable <- merge(anntable, testtable)

        if(expression)
          exprtable <- aafTableInt(eset[,grps], probeids=probeids)

        if(exists("exprtable"))
          anntable <- merge(anntable, exprtable)

        if(html)
          saveHTML(anntable, paste(filename,"html", sep="."), title)
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
                               addname = NULL, addtitle = NULL, natFC = FALSE){



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
      if(!is.null(addtitle))
          title <- paste(filename, addtitle)
      else
          title <- filename
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
            testtable <- aafTable("p-value" = ifelse(tables[[i]][,"adj.P.Val"] < 1e-3,
                                  sprintf("%0.3e", tables[[i]][,"adj.P.Val"]),
                                  sprintf("%0.3f", tables[[i]][,"adj.P.Val"])))
          }else{
            testtable <- merge(testtable,
                               aafTable("p-value" = ifelse(tables[[i]][,"adj.P.Val"] < 1e-3,
                                  sprintf("%0.3e", tables[[i]][,"adj.P.Val"]),
                                  sprintf("%0.3f", tables[[i]][,"adj.P.Val"]))))
          }
        }
      if(FC){
        fld <- tables[[i]][,"logFC"]
        if(!exists("testtable")){
          testtable <- aafTable("log2 fold Change" = round(fld, 2))
        }else{
          testtable <- merge(testtable, aafTable("log2 fold Change" = round(fld, 2)))
        }
      }
       if(natFC){
        fld <- tables[[i]][,"logFC"]
        fld <- ifelse(fld > 0, round(2^fld, 2), -round(2^-fld, 2))
          if(!exists("testtable")){
            testtable <- aafTable("Fold change natural scale" = fld)
          }else{
            testtable <- merge(testtable, aafTable("Fold change natural scale" = fld))
          }
      }  
        
      if(exists("testtable"))
        anntable <- merge(anntable, testtable)

      if(expression)
        exprtable <- aafTableInt(eset[,grps], probeids=probeids)

      if(exists("exprtable"))
        anntable <- merge(anntable, exprtable)

      if(html)
        saveHTML(anntable, paste(filename,"html", sep="."), title)
      if(text)
          saveText(anntable, paste(filename, "txt", sep="."), header=TRUE)
    }
    if(exists("testtable")) rm(testtable)
  }
  if(save) return(tables)
}



#' Filter a topTable object
#' 
#' This function is designed to filter genes from a \code{topTable} object
#' based on p-value and/or fold change. This is an internal function and is not
#' intended to be called by thte end user.
#' 
#' 
#' @param fit An \code{MArrayLM} object, resulting from a call to \code{eBayes}
#' @param coef The contrast to be extracted into the topTable. See ?topTable
#' for more information.
#' @param number The number of genes to output. Only used if both foldfilt and
#' pfilt are NULL.
#' @param fldfilt The absolute value of fold difference to filter on. This
#' assumes the data are log transformed.
#' @param pfilt The p-value to filter on.
#' @param adjust The multiplicity adjustment to use. Options are
#' '"bonferroni"', '"holm"', '"hochberg"', '"hommel"', '"fdr"' and '"none"'. If
#' '"none"' then the p-values are not adjusted. A 'NULL' value will result in
#' the default adjustment method, which is '"fdr"'.
#' @return Returns a \code{data.frame} containing the selected genes.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @keywords internal
tableFilt <- function(fit, coef = 1,  number = 30, fldfilt = NULL, pfilt = NULL,
                      adjust = "fdr"){
  if(is.null(fldfilt) && is.null(pfilt)){
    tab <- topTable(fit, coef = coef, number = number, adjust.method = adjust)
  }else{
    tab <- topTable(fit, coef = coef, number = dim(fit$coefficients)[1],
                    adjust.method = adjust)
  }
## Filter on p-value
  if(!is.null(pfilt))
      tab <- tab[tab[,"adj.P.Val"] < pfilt,]
  ## Filter on fold change
  if(!is.null(fldfilt))
    tab <- tab[abs(tab[,"logFC"]) > fldfilt,]
  tab
}
