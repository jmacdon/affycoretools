#########################################
##
## copyright 2006 James W. MacDonald
##
## Functions to integrate limma and biomaRt
##
## 6-16-2006 Preliminary functions. Still too much
##           hard coding on the biomaRt side.
##
##
##
## vennSelectBM - a function designed to take the output from a Venn diagram
##                and output tables for each cell of the Venn diagram. This emulates
##                the existing vennSelect() function, only using biomaRt for annotation data.
##
## limma2biomaRt - a function to emulate limma2annaffy, only using biomaRt for annotation
##                 data.
##
## probes2tableBM - a function to emulate probes2table()
##
## foldFiltBM - a function to emulate foldFilt()
##
##  TODO: Add functionality to output text tables in addition to HTML
##
##
######################################################

## A function used to output text tables using input designed for htmlpage()

makeText <- function(anntable, testtable, table.head, filename,
                     collapse = ",", matcollapse = "\t"){

  collapseList <- function(x, collapse, matcollapse){
    if(is.vector(x)) out <- lapply(x, paste, collapse = collapse)
    if(is.matrix(x)) out <- apply(x, 1, paste, collapse = matcollapse)
    out
  }  

  sep <- ""
  tmp <- c(anntable, testtable)
  rows <- vector("character", length(tmp[[1]]))
  for(i in seq(along = tmp)){
    rows <- paste(rows, collapseList(tmp[[i]], collapse, matcollapse), sep=sep)
    sep <- "\t"
  }
  rows <- gsub("&nbsp;","", rows)
  outfile <- file(paste(filename, ".txt", sep = ""), "w")
  cat(paste(table.head, collapse = "\t"),"\n",file = outfile, sep = "")
  cat(rows, file = outfile, sep = "\n")
  close(outfile)
}


## A function to either list out the sort of annotation that can be used to make hyperlinks
## if called with no arguments. Alternatively, it can be used to get the correct attribute
## names to pass to getBM(), after ensuring that these attributes exist for the given species,
## as well as getting the correct names to pass to htmlpage().

linksBM <- function(mart, annot, affyid = FALSE, ann.source = NULL){
  choices <- list("UniGene" = "unigene",
                  "RefSeq" = "refseq_dna",
                  "EntrezGene" = "entrezgene",
                  "SwissProt" = "uniprot_swissprot_accession",
                  "OMIM" = "adf_omim")
  if(missing(mart))
    return(names(choices))
  if(!class(mart) == "Mart")
    stop("The 'mart' argument must be a connection to a Biomart.")
  if(missing(annot))
    return(names(choices[unlist(choices) %in% listAttributes(mart)[,1]]))
  else{
    notchoice <- !annot %in% names(choices)
    if(any(notchoice))
      stop(paste("'",paste(annot[notchoice], collapse = " and "),
                 "' is not something I can use to make a hyperlink.\nPlease use linksBM()",
                 "to see a list of valid choices.", sep = ""), call. = FALSE)
  }

  ## check to see if the choices are available
  tmp <- unlist(choices[annot])
  tst <- tmp %in% listAttributes(mart)[,1]
  if(any(!tst))
    warning(paste("The following annotation sources are not available at this mart\n",
                  "for this species and were not used:",
                  paste(names(choices[annot][!tst]), collapse = ", ")), call. = FALSE)
  tmp <- tmp[tst]
  annot <- annot[tst]
  ## get repositories for htmlpage
  repository <- vector()
  
  for(i in tmp){
    tmp2 <-  switch(i,
                    unigene = "ug",
                    refseq_dna = "gb",
                    entrezgene = "en",
                    uniprot_swissprot_accession = "sp",
                    adf_omim = "omim")
    repository <- c(repository, tmp2)
  }

  if(affyid){
    ## kludge to account for the fact that the filter != attribute
    ## for this particular chip
    ## because of this we also need to output the corrected ann.source
    an.src <- switch(ann.source,
                     affy_hg_u133a_2="affy_hg_u133a_v2",
                     ann.source)
    tmp <- c(an.src, tmp)
    annot <- c("Affymetrix ID", annot)
    repository <- c("affy", repository)
  }
  
  out <- list(names = annot, repository = repository, links = tmp,
              ann.source = an.src)
  out
}


getSymbol <- function(mart){
  attributes <- listAttributes(mart)[,1]
  symbol <- grep("symbol", attributes, value=TRUE)
  if(length(symbol) == 0){
    symbol <- grep("hgnc", attributes, value=TRUE)
    if(length(symbol) == 0)
        symbol <- "nosymbolavailable"
  }
  symbol
}

  

## A function to either list out the sort of annotation that is available, but not
## 'hyperlinkable', or to get the correct attribute names to pass to getBM() and
## htmlpage(). As with linksBM(), we also check to ensure a given attribute is
## available for the given species.

annBM <- function(mart, annot){
  choices <- list("Symbol" = "nothingyet",
                  "Description" = "description",
                  "GO" = "go_description",
                  "GOID" = "go",
                  "Chromosome" = "chromosome_name",
                  "ChromLoc" = "chromosome_location")
  if(missing(mart))
    return(names(choices))
  if(!class(mart) == "Mart")
    stop("The 'mart' argument must be a connection to a Biomart.")
  if(missing(annot)){
    choices$Symbol <- getSymbol(mart)
    return(names(choices[unlist(choices) %in% listAttributes(mart)[,1]]))
  }
  else{
    choices$Symbol <- getSymbol(mart)
    notchoice <- !annot %in% names(choices)
    if(any(notchoice))
      stop(paste("'",paste(annot[notchoice], collapse = " and "),
                 "' is not something I can use to make a hyperlink.\nPlease use annBM()",
                 "to see a list of valid choices.", sep = ""), call. = FALSE)
  }
  repository <- vector()
  ## check to see if the choices are available
  tmp <- unlist(choices[annot])
  tst <- tmp %in% listAttributes(mart)[,1]
  if(any(!tst))
    warning(paste("The following annotation sources are not available at this mart\n",
                  "for this species and were not used:",
                  paste(names(choices[annot][!tst]), collapse = ", ")), call. = FALSE)
  out <- list(names = annot[tst], links = tmp[tst])
  out
}



vennSelectBM <- function (eset, design, x, contrast, fit, method = "same", adj.meth = "BH",
                          stat = "fstat", otherstats = c("pval", "FC"), order.by = "pval",
                          foldFilt = NULL, save = FALSE, species, links = linksBM()[1:3],
                          otherdata = annBM()[1:3], mysql = TRUE, ann.source = "entrezgene",
                          html = TRUE, text = TRUE, affyid = FALSE, ...){

  require("limma", quietly = TRUE)
  require("annotate", quietly= TRUE)
  require("biomaRt", quietly = TRUE)
  
  mart <- useMart("ensembl", dataset = paste(species, "_gene_ensembl", sep=""), mysql = mysql)

  ## check to see if the ann.source is available

  if(!ann.source %in% listFilters(mart)[,1]){
    cat(paste("Error: '", ann.source, "'is not an available annotation source for",
              "this biomaRt or this species.\nAvailable choices are listed below:\n"))
    return(listFilters(mart))
  }
   

  ## set up data to retrieve

  links <- linksBM( mart, links, affyid, ann.source)
  otherdata <- annBM(mart, otherdata)
  
  if (!is.null(foldFilt)) {
    idx <- abs(fit$coefficients) > foldFilt
    x <- as.matrix(x) * idx
  }
  ncontrasts <- ncol(x)
  if (ncontrasts < 2 || ncontrasts > 3) 
    stop("This function only works for two or three comparisons at a time.\n", 
         call. = FALSE)
  if (ncontrasts == 2) 
    name <- c(paste("Genes unique to", colnames(x)), paste("Genes in intersection of", 
                                                           intNames(x)))
  if (ncontrasts == 3) 
    name <- c(paste("Genes unique to", colnames(x)), paste("Genes in intersection of", 
                                                           intNames(x)), "Genes common to all comparisons")
  ## Remove illegal characters from filenames
  if(length(grep("[/|\\|?|*|:|<|>|\"|\\|]", name)) > 0)
    warning(paste("Some illegal characters have been removed from the filenames",
                  name, sep = " "), call. = FALSE)
  name <- gsub("[/|\\|?|*|:|<|>|\"|\\|]", "", name)
  
  indices <- makeIndices(x, method = method)
  cols <- getCols(design, contrast)
  for (i in seq(along = indices)) {
    tmp <- featureNames(eset)[indices[[i]]]
    if(stat == "fstat")
      stats <- f.stat.dat(fit, indices[[i]], contrast, adj.meth, c(stat, otherstats),
                              order.by, ncontrasts, i)
    if(stat == "tstat")
      stats <- t.stat.dat(fit, indices[[i]], contrast, adj.meth, c(stat, otherstats),
                              order.by, ncontrasts, i)
    if(is.null(stat)) stats <- NULL
    if (length(tmp) == 0){ 
      next
    }else{
      if(!affyid){
        if(!is.null(stats))
          gn <- sub("_at$", "", tmp[stats$ord])
        else
          gn <- sub("_at$", "", tmp)
      }else{
        if(!is.null(stats))
          gn <- tmp[stats$ord]
        else
          gn <- tmp
      }
      lnks <- getBM(attributes = links$links,
                    filter = ann.source, values = gn, mart = mart, output = "list",
                    na.value = "&nbsp;")
      ## If there isn't any info for a particular gene at the Biomart, an empty string
      ## is returned, so we need to put the original IDs into the output
      lnks[[match(links$ann.source, names(lnks))]] <- gn
      
      nam <- getBM(attributes = otherdata$links,
                   filter = ann.source, values = gn, mart = mart, output = "list",
                   na.value = "&nbsp;")
      table.head <- c(links$names, otherdata$names)
      if(!is.null(stats)){
        nam <- c(nam, stats$out)
        nam$expression <- round(exprs(eset[,cols[[i]]])[tmp[stats$ord], , drop = FALSE], 3)
        table.head <- c(table.head, names(stats$out))
      }else{
        nam$expression <- round(exprs(eset[,cols[[i]]])[tmp, , drop = FALSE], 3)
      }
      if(html)
        htmlpage(lnks, paste(name[i], "html", sep="."), name[i], nam,
                 c(table.head,  sampleNames(eset)[cols[[i]]]),
                 repository = as.list(links$repository))
      if(text)
        makeText(lnks, nam, c(table.head, sampleNames(eset)[cols[[i]]]), name[i])
    }
  }
  martDisconnect(mart)
  if (save) 
    sapply(indices, sum)
}


limma2biomaRt <- function (eset, fit, design, contrast, species, links = linksBM()[1:3],
                           otherdata = annBM()[1:3], ann.source = "entrezgene", adjust = "fdr",
                           number = 30, pfilt = NULL, fldfilt = NULL, tstat = TRUE,
                           pval = TRUE, FC = TRUE, expression = TRUE, html = TRUE, text = TRUE,
                           save = FALSE, addname = NULL, interactive = TRUE, affyid = FALSE,
                           mysql = TRUE){

  
  if(!interactive){
    limma2biomaRt.na(eset=eset, fit=fit, design=design, contrast=contrast, species=species,
                     links=links, otherdata=otherdata, ann.source=ann.source,
                     adjust=adjust, number=number, pfilt=pfilt, fldfilt=fldfilt,
                     tstat=tstat, pval=pval, FC=FC, expression=expression,
                     html=html, save=save, addname=addname, affyid=affyid, mysql=mysql)
  }else{
    require(biomaRt, quietly = TRUE)
    mart <- useMart("ensembl", dataset = paste(species, "_gene_ensembl", sep=""),
                    mysql = mysql)

    ## check to see if ann.source is available

    if(!ann.source %in% listFilters(mart)[,1]){
      cat(paste("Error: '", ann.source, "'is not an available annotation source for",
                "this biomaRt or this species.\nAvailable choices are listed below:\n"))
      return(listFilters(mart))
    }
    
 
    ## Set up default data to retrieve
    links <- linksBM(mart, links, affyid, ann.source)
    otherdata <- annBM(mart, otherdata)
    
    tables <- vector("list", dim(contrast)[2])
    for (i in seq(along = colnames(contrast))) {
      tables[[i]] <- tableFilt(fit, coef = i, number = number, 
                               pfilt = pfilt, fldfilt = fldfilt, adjust = adjust)
    }
    rn <- vector()
    for (i in 1:length(tables)) rn[i] <- dim(tables[[i]])[1]
    cat("\nYou are going to output", length(tables), "tables,", 
        "\nWith this many genes in each:\n", rn)
    cat("\n")
    doit <- getans2("Do you want to accept or change these values?", 
                    allowed = c("a", "c"))
    if (doit == "c") {
      dopval <- getans2("Do you want to change the p-value?", 
                        allowed = c("y", "n"))
      if (dopval == "y") 
        pfilt <- getpfil()
      dofld <- getans2("Do you want to change the fold filter?", 
                       allowed = c("y", "n"))
      if (dofld == "y") 
        fldfilt <- getfldfil()

      ## adjust links and otherdata
      links <- links[[1]]
      if(any(links == "Affymetrix ID"))
        links <- links[-which(links == "Affymetrix ID")]
      otherdata <- otherdata[[1]]
      
      limma2biomaRt(eset, fit, design, contrast, species, links,
                    otherdata, ann.source, adjust, number, pfilt, fldfilt,
                    tstat, pval, FC, expression, html, text, save, addname, interactive,
                    affyid, mysql)
      options(show.error.messages = FALSE)
      stop()
    }
   
    for (i in 1:length(tables)) {
      if (dim(tables[[i]])[1] > 0) {
        index <- as.numeric(row.names(tables[[i]]))
        if (length(which(contrast[, i] > 0)) == 1) {
          grp1 <- which(design[, which(contrast[, i] > 
                                       0)] > 0)
        }
        else {
          grp1 <- which(rowSums(design[, which(contrast[, 
                                                        i] > 0)]) > 0)
        }
        if (length(which(contrast[, i] < 0)) == 1) {
          grp2 <- which(design[, which(contrast[, i] < 
                                       0)] > 0)
        }
        else {
          grp2 <- which(rowSums(design[, which(contrast[, 
                                                        i] < 0)]) > 0)
        }
        grps <- c(grp1, grp2)
        filename <- colnames(contrast)[i]
        if (!is.null(addname)) 
          filename <- paste(filename, addname, sep = " ")
        if (length(grep("[/|\\|?|*|:|<|>|\"|\\|]", filename)) > 
            0) 
          warning(paste("Some illegal characters have been removed from the filename", 
                        filename, sep = " "), call. = FALSE)
        filename <- gsub("[/|\\|?|*|:|<|>|\"|\\|]", "", 
                         filename)
        probeids <- featureNames(eset)[index]
        if(affyid)
          gn <- probeids
        else
          gn <- sub("_at$", "", probeids)
        
        anntable <- getBM(attributes = links$links, filter = ann.source,
                          values = gn,
                            mart = mart, output = "list", na.value = "&nbsp;")
        ## Need to dump the 'values' into the anntable list - if no data at the mart,
        ## biomaRt just returns an empty string
        
        anntable[[match(links$ann.source, names(anntable))]] <-  gn
        
        testtable <- getBM(attributes = otherdata$links, filter = ann.source,
                           values = gn, mart = mart, output = "list",
                           na.value = "&nbsp;")
        
        if (tstat)
          testtable$t.statistic <-  round(tables[[i]][,"t"], 2)
        
        if (pval) 
          testtable$p.value <- round(tables[[i]][,"adj.P.Val"], 3)
        
        if (FC) {
          fld <- tables[[i]][, 2]
          testtable$fold.change <- round(fld,2)
        }
        if (expression) 
          testtable$expression <- round(exprs(eset[, grps])[probeids, , drop = FALSE], 3)
      }
      table.head <- c(links$names, otherdata$names,
                      c("t-statistic","p-value","Fold change")[c(tstat, pval, FC)],
                      sampleNames(eset)[grps])
      if(html)
        htmlpage(anntable, paste(filename, "html", sep = "."), filename, testtable,
                 table.head, repository = as.list(links$repository))
      if(text)
        makeText(anntable, testtable, table.head, filename)
      
      }
    
    martDisconnect(mart)
    if (save) 
      return(tables)
  }
}


limma2biomaRt.na <- function (eset, fit, design, contrast, species, links = linksBM()[1:3],
                              otherdata = annBM()[1:3], ann.source = "entrezgene", adjust = "fdr",
                              number = 30, pfilt = NULL, fldfilt = NULL, tstat = TRUE,
                              pval = TRUE, FC = TRUE, expression = TRUE, html = TRUE, text = TRUE,
                              save = FALSE, addname = NULL, affyid = FALSE, mysql = TRUE){


  require(biomaRt, quietly = TRUE)
  mart <- useMart("ensembl", dataset = paste(species, "_gene_ensembl", sep=""), mysql = mysql)

  ## check to see if ann.source is available

  if(!ann.source %in% listFilters(mart)[,1]){
    cat(paste("Error: '", ann.source, "'is not an available annotation source for",
              "this biomaRt or this species.\nAvailable choices are listed below:\n"))
    return(listFilters(mart))
  }
  
  ## Set up default data to retrieve
  links <- linksBM(mart, links, affyid, ann.source)
  otherdata <- annBM(mart, otherdata)
  
  tables <- vector("list", dim(contrast)[2])
  for (i in seq(along = colnames(contrast))) {
    tables[[i]] <- tableFilt(fit, coef = i, number = number, 
                             pfilt = pfilt, fldfilt = fldfilt, adjust = adjust)
  }
  rn <- vector()
  
  
  for (i in 1:length(tables)) {
    if (dim(tables[[i]])[1] > 0) {
      index <- as.numeric(row.names(tables[[i]]))
      if (length(which(contrast[, i] > 0)) == 1) {
        grp1 <- which(design[, which(contrast[, i] > 
                                     0)] > 0)
      }
      else {
        grp1 <- which(rowSums(design[, which(contrast[, 
                                                      i] > 0)]) > 0)
      }
      if (length(which(contrast[, i] < 0)) == 1) {
        grp2 <- which(design[, which(contrast[, i] < 
                                     0)] > 0)
      }
      else {
        grp2 <- which(rowSums(design[, which(contrast[, 
                                                      i] < 0)]) > 0)
      }
      grps <- c(grp1, grp2)
      filename <- colnames(contrast)[i]
      if (!is.null(addname)) 
        filename <- paste(filename, addname, sep = " ")
      if (length(grep("[/|\\|?|*|:|<|>|\"|\\|]", filename)) > 
          0) 
        warning(paste("Some illegal characters have been removed from the filename", 
                      filename, sep = " "), call. = FALSE)
      filename <- gsub("[/|\\|?|*|:|<|>|\"|\\|]", "", 
                       filename)
      probeids <- featureNames(eset)[index]

      if(affyid)
        gn <- probeids
      else
        gn <- sub("_at", "", probeids)
      anntable <- getBM(attributes = links$links, filter = ann.source,
                        values = gn, mart = mart, output = "list", na.value = "&nbsp;")
      ## Need to dump the 'values' into the anntable list - if no data at the mart,
      ## biomaRt just returns an empty string

      anntable[[match(links$ann.source, names(anntable))]] <-  gn
      
      testtable <- getBM(attributes = otherdata$links, filter = ann.source,
                         values = gn, mart = mart, output = "list",
                         na.value = "&nbsp;")
      
      if (tstat)
        testtable$t.statistic <-  round(tables[[i]][,"t"], 2)
      
      if (pval) 
        testtable$p.value <- round(tables[[i]][,"adj.P.Val"], 3)
      
      if (FC) {
        fld <- tables[[i]][, 2]
        testtable$fold.change <- round(fld,2)
      }
      if (expression) 
        testtable$expression <- round(exprs(eset[, grps])[probeids, , drop = FALSE], 3)
    }
    table.head <- c(links$names, otherdata$names,
                    c("t-statistic","p-value","Fold change")[c(tstat, pval, FC)],
                    sampleNames(eset)[grps])
    if(html)
      htmlpage(anntable, paste(filename, "html", sep = "."), filename, testtable,
               table.head, repository = as.list(links$repository))
    if(text)
      makeText(anntable, testtable, table.head, filename)
  }
  
  martDisconnect(mart)
  if (save) 
    return(tables)
}

probes2tableBM <- function(eset, probids, species, filename, otherdata = NULL,
                           links = linksBM()[1:3], otherann = annBM()[1:3],
                           ann.source = "entrezgene", express = TRUE, html = TRUE,
                           text = TRUE, affyid = FALSE, mysql = TRUE){
  require(biomaRt, quietly = TRUE)
  mart <- useMart("ensembl", dataset = paste(species, "_gene_ensembl", sep=""), mysql = mysql)
  
  ## check to see if ann.source is available
  
  if(!ann.source %in% listFilters(mart)[,1]){
    cat(paste("Error: '", ann.source, "'is not an available annotation source for",
              "this biomaRt or this species.\nAvailable choices are listed below:\n"))
    return(listFilters(mart))
  }
  
  ## Set up default data to retrieve
  links <- linksBM(mart, links, affyid, ann.source)
  otherann <- annBM(mart, otherann)

  ## get link data

  if(affyid)
    gn <- probids
  else
    gn <- sub("_at", "", probids)
  anntable <-  getBM(attributes = links$links, filter = ann.source,
                     values = gn, mart = mart, output = "list", na.value = "&nbsp;")
  ## Need to dump the 'values' into the anntable list - if no data at the mart,
  ## biomaRt just returns an empty string
  
  anntable[[match(links$ann.source, names(anntable))]] <-  gn
  
  testtable <- getBM(attributes = otherann$links, filter = ann.source,
                     values = gn, mart = mart, output = "list",
                     na.value = "&nbsp;")
  table.head <- c(links$names, otherann$names, names(otherdata),
                  sampleNames(eset))
  if(!is.null(otherdata))
    testtable <- c(testtable, otherdata)
  if(express)
    testtable <- c(testtable, as.data.frame(round(exprs(eset)[probids, , drop = FALSE], 2)))
  if(html)
    htmlpage(anntable, paste(filename, "html", sep = "."), filename, testtable,
             table.head, repository = as.list(links$repository))
  if(text)
    makeText(anntable, testtable, table.head, filename)

  martDisconnect(mart)
 }

foldFiltBM <- function(object, fold = 1, groups, comps, compnames, species,
                       links = linksBM()[1:3], otherann = annBM()[1:3], filterfun = NULL,
                       ann.source = "entrezgene", affyid = FALSE, mysql = TRUE, html = TRUE,
                       text = TRUE, save = FALSE){
 
  if(is(object, "exprSet") || is(object, "ExpressionSet"))
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
  flds <- lapply(comps, function(y) round(gps[,y[1]] - gps[,y[2]], 2))
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
      probes2tableBM(object[,idx], probes[[i]][ord], species, compnames[i],
                     list("Fold change" = FCs[[i]][ord]), links, otherann,
                     ann.source, TRUE, html, text, affyid, mysql)
    }
  }
  direct <- matrix(NA, nc = length(comps), nr = dim(gps)[1])
  for(i in seq(along = flds)){
    direct[,i] <- sign(flds[[i]] * indices[[i]])
  }
  if(save)
   return(list(sums =  sapply(indices, sum), dirs = direct))
}
