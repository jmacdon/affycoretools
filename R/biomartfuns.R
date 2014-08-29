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

## A function used to convert data.frame output from biomaRt to a list
## useful for htmlpage()

dfToList <- function(dataframe, dataToUse, gn){
    width <- dim(dataframe)[2]
    len <- dim(dataframe)[1]
    tmp <- lapply(1:width, function(x)
                  tapply(1:len, dataframe[,dataToUse],  function(y) dataframe[y,x]))
    tmp <- lapply(tmp, function(x) lapply(x, unique))
    ## biomaRt doesn't have annotation for all input values, so we have to bump
    ## this out and put &nbsp; in where we have no data. To do this, we always
    ## query biomaRt for the input values as well, and compare what we sent
    ## to what we got to figure out what is missing.
    ind <- match(gn, tmp[[dataToUse]])
    ind <- ind[!is.na(ind)]
    ind2 <- gn %in% tmp[[dataToUse]]
    out <- vector("list",width)
    for(i in seq_len(width)[-dataToUse]){
        out[[i]] <- vector("character", length = length(gn))
        out[[i]][ind2] <- tmp[[i]][ind]
        out[[dataToUse]] <- gn
    }
    out <- lapply(out, function(x) ifelse(x == "", "&nbsp;", x))
    out
}
                  

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
    ann.source <- switch(ann.source,
                     affy_hg_u133a_2="affy_hg_u133a_v2",
                     ann.source)
    tmp <- c(ann.source, tmp)
    annot <- c("Affymetrix ID", annot)
    repository <- c("affy", repository)
  }
  
  out <- list(names = annot, repository = repository, links = tmp,
              ann.source = ann.source)
  out
}


  

## A function to either list out the sort of annotation that is available, but not
## 'hyperlinkable', or to get the correct attribute names to pass to getBM() and
## htmlpage(). As with linksBM(), we also check to ensure a given attribute is
## available for the given species.



#' Select Available Annotation from a Biomart
#' 
#' These functions are designed to do two things that are useful for an end
#' user. If called with no arguments, they will output a character vector of
#' annotation sources that are typically available from a Biomart database. If
#' called with a 'mart' connection (typically created by a call to
#' \code{\link[biomaRt]{useMart}}), they will return a character vector of
#' annotation sources that exist for that particular Biomart and species. If
#' called with a 'mart' connection and a character vector of annotation
#' sources, they will return a list that is intended to be used by other
#' functions for creating HTML pages. This last function doesn't have any real
#' utility for the end user.
#' 
#' The purpose of these functions is to either give an example of typical
#' annotation sources that may be available at a particular Biomart, or to
#' output those sources that are known to exist at a Biomart.
#' 
#' \code{linksBM} is intended to list those annotation sources that may be
#' turned into hyperlinks whereas \code{annBM} is intended to list those
#' annotation sources that will not be linked.
#' 
#' These functions have only a few of the possible annotation sources, and
#' currently there is no simple way to extend these sources. Additions to the
#' list are possible, however. Please contact me if there is something in
#' particular that should be included in either list.
#' 
#' @aliases annBM linksBM
#' @param mart A 'mart' connection, typically created by a call to
#' \code{\link[biomaRt]{useMart}}.
#' @param annot A character vector of annotation sources. This is not typically
#' useful for an end user to specify.
#' @param species A species name, of the form e.g., 'hsapiens'
#' @return Normally called by an end user to output a character vector of
#' annotation sources.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @keywords manip
#' @examples
#' 
#' annBM()
#' 
#' @export annBM
annBM <- function(mart, annot, species){
  if(!missing(species)) symbol <- switch(species,
                                        "hsapiens" = "hgnc_symbol",
                                        "mmusculus" = "markersymbol",
                                        "rnorvegicus" = "markersymbol",
                                        "external_gene_id")
  else symbol <- "nothingyet"
  
  choices <- list("Symbol" = symbol,
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
    return(names(choices[unlist(choices) %in% listAttributes(mart)[,1]]))
  }
  else{
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





#' Select and Output Genelists Based on Venn Diagrams using biomaRt
#' 
#' This function is designed to output HTML and text tables based on the
#' results of a call to \code{\link[limma]{decideTests}}. This function is very
#' similar to \code{vennSelect}, except it uses the \code{biomaRt} package to
#' annotate genes, and the \code{annotate} package to create the HTML table.
#' 
#' The purpose of this function is to output HTML tables with lists of genes
#' that fulfill the criteria of a call to \code{\link[limma]{decideTests}} as
#' well as the direction of differential expression.
#' 
#' The IDs that will be used to annotate the genes depend on the source of the
#' data. If, for example, one is using an Affymetrix chip that doesn't have a
#' BioC annotation package, then the IDs will be Affymetrix IDs. To find out
#' the correct name to use for the ann.source argument, one can create a
#' connection to a Biomart database using \code{\link[biomaRt]{useMart}} and
#' then deduce the correct argument by the output from
#' \code{listFilters(mart)}. It will usually be something starting with 'affy',
#' and contain the name of the chip.
#' 
#' If one is using one of the re-mapped CDFs from MBNI at University of
#' Michigan, then the IDs to use depend on the mapping used to create the CDF.
#' At this time, only three types of CDFs can be used; EntrezGene, UniGene, and
#' RefSeq. One can determine the correct ann.source argument by creating a
#' connection to a Biomart database, and then calling \code{linksBM(mart,
#' linksBM())[[3]]}.
#' 
#' Some important things to note: First, the names of the HTML tables are
#' extracted from the \code{colnames} of the \code{TestResults} object, which
#' come from the contrasts matrix, so it is important to use something
#' descriptive. Second, the method argument is analogous to the \code{include}
#' argument from \code{\link[limma:venn]{vennCounts}} or
#' \code{\link[limma:venn]{vennDiagram}}. Choosing "both" will select genes
#' that are differentially expressed in one or more comparisons, regardless of
#' direction. Choosing "up" or "down" will select genes that are only
#' differentially expressed in one direction. Choosing "same" will select genes
#' that are differentially expressed in the same direction. Choosing "sameup"
#' or "samedown" will select genes that are differentially expressed in the
#' same direction as well as 'up' or 'down'.
#' 
#' Note that this is different than sequentially choosing "up" and then "down".
#' For instance, a gene that is upregulated in one comparison and downregulated
#' in another comparison will be listed in the intersection of those two
#' comparisons if "both" is chosen, it will be listed in only one comparison
#' for both the "up" and "down" methods, and it will be listed in the union
#' (e.g., not selected) if "same" is chosen.
#' 
#' Calling the function normally will result in the output of HTML tables.
#' 
#' Calling the function with save set to \code{TRUE} will output HTML tables as
#' well as a vector of counts for each comparison. This is useful when using
#' the function programmatically (e.g., when making reports using Sweave).
#' 
#' out <- vennSelectBM(eset, fit, design, x, <other arguments>, save = TRUE)
#' 
#' An alternative would be to use \code{vennCounts2} and
#' \code{\link[limma:venn]{vennDiagram}} to output a Venn diagram, which is
#' probably more reasonable since the tables being output are supposed to be
#' based on a Venn diagram.
#' 
#' @param eset A \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}
#' object.
#' @param design A \code{design} matrix, usually from a call to
#' \code{model.matrix}. See details for more information.
#' @param x A \code{\link[limma]{TestResults}} object, usually from a call to
#' \code{\link[limma]{decideTests}}.
#' @param contrast A contrasts matrix, produced either by hand, or by a call to
#' \code{\link[limma]{makeContrasts}}
#' @param fit An \code{\link[limma:marraylm]{MArrayLM}} object, from a call to
#' \code{\link[limma:ebayes]{eBayes}}.
#' @param method One of "same", "both", "up", "down", "sameup", or "samedown".
#' See details for more information.
#' @param adj.meth Method to use for adjusting p-values. Default is 'BH', which
#' corresponds to 'fdr'. Ideally one would set this value to be the same as was
#' used for \code{\link[limma]{decideTests}}.
#' @param stat The statistic to report in the resulting HTML tables. Choices
#' are 'fstat', 'tstat', and \code{NULL}. Ideally, the statistic chosen would
#' correspond to the method used in \code{\link[limma]{decideTests}}. In other
#' words, if one used methods such as 'separate' or 'hierarchical', which are
#' based on a t-statistic, one should choose 'tstat', however, if one used
#' 'nestedF', the logical choice would be 'fstat'.
#' @param otherstats Other statistics to be included in the HTML tables.
#' Choices include 'pval' and 'FC'.
#' @param order.by Which statistic should be used to order the probesets?
#' Choices include 'fstat', 'tstat', 'pval', and 'FC'. Note that if 'FC' is
#' chosen and there are more than one set of fold changes, the first will be
#' used.
#' @param foldFilt A fold change to use for filtering. Default is \code{NULL},
#' meaning no filtering will be done.
#' @param save Boolean - If \code{TRUE}, output a count of genes that fulfill
#' the criteria. Useful for e.g., Sweave-type reports.
#' @param species The species name. This must be in a particular format for
#' biomaRt. An example for human is "hsapiens", or for mouse "mmusculus".
#' @param links A character vector of things to annotate with hyperlinks to
#' online databases. See \code{linksBM} for possible values.
#' @param otherdata A character vector of things to annotate with text only
#' (i.e., no hyperlinks). See \code{annBM} for possible values.
#' @param ann.source The annotation source of the IDs that will be used to
#' annotate the genes. The default value is "entrezgene". See details for other
#' possibilities.
#' @param html Boolean. Output HTML tables? Defaults to \code{TRUE}
#' @param text Boolean. Output text tables? Defaulst to \code{TRUE}
#' @param affyid Boolean. Are the IDs used to annotate these data Affymetrix
#' IDs?
#' @param ... Used to pass other variables to e.g., \code{htmlpage}.
#' @return Normally called only for the side effect of producing HTML tables.
#' However, setting save to \code{TRUE} will output a vector of counts that can
#' be used for making Sweave-style reports.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @keywords manip
#' @export vennSelectBM
vennSelectBM <- function (eset, design, x, contrast, fit, method = "same", adj.meth = "BH",
                          stat = "fstat", otherstats = c("pval", "FC"), order.by = "pval",
                          foldFilt = NULL, save = FALSE, species, links = linksBM()[1:3],
                          otherdata = annBM()[1:3], ann.source = "entrezgene",
                          html = TRUE, text = TRUE, affyid = FALSE, ...){

  .Deprecated(new = "makeVenn", msg = paste("This function is being deprecated. Please see the RefactoredAffycoretools",
              "vignette for more current ways to create annotated output."))
  mart <- useMart("ensembl", dataset = paste(species, "_gene_ensembl", sep=""))

  ## check to see if the ann.source is available

  if(!ann.source %in% listFilters(mart)[,1]){
    cat(paste("Error: '", ann.source, "'is not an available annotation source for",
              "this biomaRt or this species.\nAvailable choices are listed below:\n"))
    return(listFilters(mart))
  }
   

  ## set up data to retrieve

  links <- linksBM( mart, links, affyid, ann.source)
  otherdata <- annBM(mart, otherdata, species)
  
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
      lnks <- dfToList(getBM(attributes = links$links,
                             filters = ann.source, values = gn, mart = mart),
                       which(links$links == ann.source), gn)
           
      nam <- dfToList(getBM(attributes = c(ann.source, otherdata$links),
                            filters = ann.source, values = gn, mart = mart),
                      1, gn)[-1]
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
  
  if (save) 
    sapply(indices, sum)
}




#' Function to Create HTML Tables from limma Objects using biomaRt for
#' Annotation
#' 
#' This function is designed to take an \code{ExpressionSet} and an
#' \code{lmFit}, \code{model.matrix}, and contrast object from limma and
#' convert into HTML and text tables using biomaRt. The alternate function
#' \code{limma2biomaRt.na} is designed to be run without user intervention.
#' 
#' This function is designed to automatically output HTML tables, with
#' filenames taken from the column names of the contrast matrix. The number of
#' genes output can be controlled several different ways. First, if pfilt and
#' fldfilt are both \code{NULL}, the top genes will be output based on the
#' \code{number} variable. Otherwise, the genes are filtered based on p-value,
#' fold change, or both. If the genes are filtered this way, the number of
#' genes to be output will be listed and the filter(s) can then be adjusted if
#' necessary.
#' 
#' This function currently only supports Affymetrix data. It is designed for
#' Affymetrix chips that don't have an annotation package, which includes data
#' that have been analyzed using the 're-mapped' CDFs supplied to BioC by MBNI
#' at University of Michigan.
#' 
#' The IDs that will be used to annotate the genes depend on the source of the
#' data. If, for example, one is using an Affymetrix chip that doesn't have a
#' BioC annotation package, then the IDs will be Affymetrix IDs. To find out
#' the correct name to use for the ann.source argument, one can create a
#' connection to a Biomart database using \code{\link[biomaRt]{useMart}} and
#' then get a list of available Affy arrays using
#' \code{\link[biomaRt]{listFilters}}.
#' 
#' If one is using one of the re-mapped CDFs from MBNI at University of
#' Michigan, then the IDs to use depend on the mapping used to create the CDF.
#' At this time, only three types of CDFs can be used; EntrezGene, UniGene, and
#' RefSeq. One can determine the correct ann.source argument by creating a
#' connection to a Biomart database, and then calling \code{linksBM(mart,
#' linksBM())[[3]]}.
#' 
#' @aliases limma2biomaRt limma2biomaRt.na
#' @param eset An \code{ExpressionSet} containing affymetrix expression values.
#' @param fit An \code{lmFit} object.
#' @param design A \code{model.matrix} object.
#' @param contrast A contrasts matrix from limma.
#' @param species The species name. This must be in a particular format for
#' biomaRt. An example for human is "hsapiens", or for mouse "mmusculus".
#' @param links A character vector of things to annotate with hyperlinks to
#' online databases. See \code{linksBM} for possible values.
#' @param otherdata A character vector of things to annotate with text only
#' (i.e., no hyperlinks). See \code{annBM} for possible values.
#' @param ann.source The annotation source of the IDs that will be used to
#' annotate the genes. The default value is "entrezgene". See details for other
#' possibilities.
#' @param adjust Multiplicity adjustment. Choices are
#' "fdr","holm","hommel","bonferroni", or "none". Partial matching allowed.
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
#' @param text Boolean: Output data in text tables? Defaults to \code{TRUE}
#' @param save Boolean: Save tables as R objects for further processing?
#' Defaults to \code{FALSE}.
#' @param addname A character vector to add to the end of the automatically
#' generated output file names. Useful for multiple calls to eliminate
#' over-writing of existing HTML or text tables.
#' @param interactive Boolean: Is this an interactive call, or run as part of a
#' script (e.g., in an \code{Sweave} document)? Defaults to \code{TRUE}
#' @param affyid Boolean. Are the IDs used to annotate these data Affymetrix
#' IDs?
#' @return If \code{save} is \code{TRUE}, a list of tables from
#' \code{\link[limma:toptable]{topTable}} will be output.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @seealso \code{\link[limma:toptable]{topTable}},
#' \code{\link[annaffy]{aafTableAnn}}
#' @keywords manip
#' @export limma2biomaRt
limma2biomaRt <- function (eset, fit, design, contrast, species, links = linksBM()[1:3],
                           otherdata = annBM()[1:3], ann.source = "entrezgene", adjust = "fdr",
                           number = 30, pfilt = NULL, fldfilt = NULL, tstat = TRUE,
                           pval = TRUE, FC = TRUE, expression = TRUE, html = TRUE, text = TRUE,
                           save = FALSE, addname = NULL, interactive = TRUE, affyid = FALSE){

  .Deprecated(new = "", msg = paste("This function is being deprecated. Please see the RefactoredAffycoretools",
                        "vignette for more current ways to annotate output."))
  if(!interactive){
    limma2biomaRt.na(eset=eset, fit=fit, design=design, contrast=contrast, species=species,
                     links=links, otherdata=otherdata, ann.source=ann.source,
                     adjust=adjust, number=number, pfilt=pfilt, fldfilt=fldfilt,
                     tstat=tstat, pval=pval, FC=FC, expression=expression,
                     html=html, save=save, addname=addname, affyid=affyid)
  }else{

    mart <- useMart("ensembl", dataset = paste(species, "_gene_ensembl", sep=""))

    ## check to see if ann.source is available

    if(!ann.source %in% listFilters(mart)[,1]){
      cat(paste("Error: '", ann.source, "'is not an available annotation source for",
                "this biomaRt or this species.\nAvailable choices are listed below:\n"))
      return(listFilters(mart))
    }
    
 
    ## Set up default data to retrieve
    links <- linksBM(mart, links, affyid, ann.source)
    otherdata <- annBM(mart, otherdata, species)
    
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
                    affyid)
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
        
        anntable <- dfToList(getBM(attributes = links$links, filters = ann.source,
                                   values = gn, mart = mart),
                             which(links@links == ann.source), gn)
        
        testtable <- dfToList(getBM(attributes = c(ann.source, otherdata$links),
                                    filters = ann.source,
                                    values = gn, mart = mart), 1, gn)[-1]
        
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
    
    
    if (save) 
      return(tables)
  }
}


limma2biomaRt.na <- function (eset, fit, design, contrast, species, links = linksBM()[1:3],
                              otherdata = annBM()[1:3], ann.source = "entrezgene", adjust = "fdr",
                              number = 30, pfilt = NULL, fldfilt = NULL, tstat = TRUE,
                              pval = TRUE, FC = TRUE, expression = TRUE, html = TRUE, text = TRUE,
                              save = FALSE, addname = NULL, affyid = FALSE){



  mart <- useMart("ensembl", dataset = paste(species, "_gene_ensembl", sep=""))

  ## check to see if ann.source is available

  if(!ann.source %in% listFilters(mart)[,1]){
    cat(paste("Error: '", ann.source, "'is not an available annotation source for",
              "this biomaRt or this species.\nAvailable choices are listed below:\n"))
    return(listFilters(mart))
  }
  
  ## Set up default data to retrieve
  links <- linksBM(mart, links, affyid, ann.source)
  otherdata <- annBM(mart, otherdata, species)
  
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
      anntable <- dfToList(getBM(attributes = links$links, filters = ann.source,
                                 values = gn, mart = mart),
                           which(links$links == ann.source), gn)
      
      testtable <- dfToList(getBM(attributes = c(ann.source, otherdata$links),
                                  filters = ann.source,
                                  values = gn, mart = mart), 1, gn)[-1]
      
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
  
  
  if (save) 
    return(tables)
}



#' Convert Affy Probe ids to Annotated HTML Table using biomaRt
#' 
#' A function to convert a vector of Affy ids to an annotated HTML or text
#' table. This function is very similar to \code{probes2table}, except it uses
#' the \code{biomaRt} package to annotate genes, and the \code{annotate}
#' package to create the HTML table.
#' 
#' This function is designed to output HTML tables based on a set of IDs. This
#' function currently only supports Affymetrix data. It is designed for
#' Affymetrix chips that don't have an annotation package, which includes data
#' that have been analyzed using the 're-mapped' CDFs supplied to BioC by MBNI
#' at University of Michigan.
#' 
#' The IDs that will be used to annotate the genes depend on the source of the
#' data. If, for example, one is using an Affymetrix chip that doesn't have a
#' BioC annotation package, then the IDs will be Affymetrix IDs. To find out
#' the correct name to use for the ann.source argument, one can create a
#' connection to a Biomart database using \code{\link[biomaRt]{useMart}} and
#' then get a list of available Affy arrays using
#' \code{\link[biomaRt]{listFilters}}.
#' 
#' If one is using one of the re-mapped CDFs from MBNI at University of
#' Michigan, then the IDs to use depend on the mapping used to create the CDF.
#' At this time, only three types of CDFs can be used; EntrezGene, UniGene, and
#' RefSeq. One can determine the correct ann.source argument by creating a
#' connection to a Biomart database, and then calling \code{linksBM(mart,
#' linksBM())[[3]]}.
#' 
#' @param eset An \code{ExpressionSet} containing Affy expression values.
#' @param probids A vector of probe ids.
#' @param species The species name. This must be in a particular format for
#' biomaRt. An example for human is "hsapiens" or for mouse is "mmusculus".
#' @param filename File name of the resulting HTML table.
#' @param otherdata A *named* list of additional information to include in the
#' resulting table. Examples would be t-statistics, p-values, fold change, etc.
#' Each list item should be a vector the same length as the probids vector. The
#' name associated with each list item will be used as the column name in the
#' resulting table.
#' @param links A character vector of things to annotate with hyperlinks to
#' online databases. See \code{linksBM} for possible values.
#' @param otherann A character vector of things to annotate with text only
#' (i.e., no hyperlinks). See \code{annBM} for possible values.
#' @param ann.source The annotation source of the IDs that will be used to
#' annotate the genes. The default value is "entrezgene". See details for other
#' possibilities.
#' @param express Output expression values in table?  Defaults to \code{TRUE}.
#' @param html Boolean. Output HTML table? Defaults to \code{TRUE}
#' @param text Boolean. Output text table? Defaults to \code{TRUE}
#' @param affyid Boolean. Are the IDs used to annotate these data Affymetrix
#' IDs?
#' @return This function is only used for the side effect of outputting an HTML
#' table.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @seealso \code{\link[limma:toptable]{topTable}}
#' @keywords manip
#' @export probes2tableBM
probes2tableBM <- function(eset, probids, species, filename, otherdata = NULL,
                           links = linksBM()[1:3], otherann = annBM()[1:3],
                           ann.source = "entrezgene", express = TRUE, html = TRUE,
                           text = TRUE, affyid = FALSE){
  .Deprecated(new= "", msg = paste("This function is being deprecated. Please see the RefactoredAffycoretools",
                       "vignette for more current ways to annotate output."))
  mart <- useMart("ensembl", dataset = paste(species, "_gene_ensembl", sep=""))
  
  ## check to see if ann.source is available
  
  if(!ann.source %in% listFilters(mart)[,1]){
    cat(paste("Error: '", ann.source, "'is not an available annotation source for",
              "this biomaRt or this species.\nAvailable choices are listed below:\n"))
    return(listFilters(mart))
  }
  
  ## Set up default data to retrieve
  links <- linksBM(mart, links, affyid, ann.source)
  otherann <- annBM(mart, otherann, species)

  ## get link data

  if(affyid)
    gn <- probids
  else
    gn <- sub("_at", "", probids)
  anntable <-  dfToList(getBM(attributes = links$links, filters = ann.source,
                              values = gn, mart = mart),
                        which(links$links == ann.source), gn)
  
  testtable <- dfToList(getBM(attributes = c(ann.source, otherann$links),
                              filters = ann.source,
                              values = gn, mart = mart), 1, gn)[-1]
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

  
 }



#' Output Fold Change Data using biomaRt
#' 
#' This function is designed to take an \code{ExpressionSet} and some
#' comparisons and output HTML tables. It is very similar to \code{foldFilt}
#' except it uses the \code{biomaRt} package to annotate genes and the annotate
#' package to create the HTML table(s).
#' 
#' This function is useful for outputting annotated gene lists for multiple
#' fold change comparisons. The genes will be ordered by the absolute fold
#' change.
#' 
#' This function currently only supports Affymetrix data. It is designed for
#' Affymetrix chips that don't have an annotation package, which includes data
#' that have been analyzed using the 're-mapped' CDFs supplied to BioC by MBNI
#' at University of Michigan.
#' 
#' The IDs that will be used to annotate the genes depend on the source of the
#' data. If, for example, one is using an Affymetrix chip that doesn't have a
#' BioC annotation package, then the IDs will be Affymetrix IDs. To find out
#' the correct name to use for the ann.source argument, one can create a
#' connection to a Biomart database using \code{\link[biomaRt]{useMart}} and
#' then get a list of available Affy arrays using \code{getAffyArrays}.
#' 
#' If one is using one of the re-mapped CDFs from MBNI at University of
#' Michigan, then the IDs to use depend on the mapping used to create the CDF.
#' At this time, only three types of CDFs can be used; EntrezGene, UniGene, and
#' RefSeq. One can determine the correct ann.source argument by creating a
#' connection to a Biomart database, and then calling \code{linksBM(mart,
#' linksBM())[[3]]}.
#' 
#' One can also protect against selecting probesets that have very small
#' expression values for all samples (which likely have a large fold change due
#' to noise, rather than signal) by using the filterfun argument. An example
#' would be:
#' 
#' f <- kOverA(1, 6)
#' 
#' filt <- filterfun(f)
#' 
#' Then add filterfun = filt as an argument to the call to \code{foldFilt}.
#' 
#' @param object An \code{ExpressionSet} object
#' @param fold The log fold change cutoff to use. Note that this is log base
#' two.
#' @param groups A vector of group identifiers. Probably easiest to use a
#' numeric vector
#' @param comps A list containing all the comparisons to be made. Each list
#' item should be a vector of length two. See details for more information.
#' @param compnames A character vector of the names for each of the comparisons
#' to be made. This will be the name of the resulting HTML or text file.
#' @param species The species name. This must be in a particular format for
#' biomaRt. An example for human is "hsapiens" or for mouse is "mmusculus".
#' @param links A character vector of things to annotate with hyperlinks to
#' online databases. See \code{linksBM} for possible values.
#' @param otherann A character vector of things to annotate with text only
#' (i.e., no hyperlinks). See \code{annBM} for possible values.
#' @param filterfun A filtering function created by
#' \code{\link[genefilter]{genefilter}} to filter the data using additional
#' criteria. See details for more information
#' @param ann.source The annotation source of the IDs that will be used to
#' annotate the genes. The default value is "entrezgene". See details for other
#' possibilities.
#' @param affyid Boolean. Are the IDs used to annotate these data Affymetrix
#' IDs?
#' @param html Boolean. Output HTML tables? Defaults to \code{TRUE}
#' @param text Boolean. Output text tables? Defaults to \code{TRUE}
#' @param save Boolean. If \code{TRUE}, a list will be returned. The first item
#' in the list will be a vector showing the number of 'significant' genes for
#' each comparison. The second item will be a matrix of -1's, 0's and 1's
#' indicating a significant difference, and the direction of the difference.
#' The first item is useful for creating Sweave - based reports and the second
#' is useful for making Vennn diagrams using \code{vennDiagram} from the limma
#' package.
#' @return Returns a list; see above for the elements of the list. This
#' function is mainly called for the side effect of outputting HTML or text
#' files containing annotated 'significant' gene lists.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @keywords manip
#' @export foldFiltBM
foldFiltBM <- function(object, fold = 1, groups, comps, compnames, species,
                       links = linksBM()[1:3], otherann = annBM()[1:3], filterfun = NULL,
                       ann.source = "entrezgene", affyid = FALSE, html = TRUE,
                       text = TRUE, save = FALSE){
  .Deprecated(new = "", msg = paste("This fucntion is being deprecated. Please see the RefactoredAffycoretools",
                        "vignette for more current ways to annotate output."))
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
                     ann.source, TRUE, html, text, affyid)
    }
  }
  direct <- matrix(NA, ncol = length(comps), nrow = dim(gps)[1])
  for(i in seq(along = flds)){
    direct[,i] <- sign(flds[[i]] * indices[[i]])
  }
  if(save)
   return(list(sums =  sapply(indices, sum), dirs = direct))
}
