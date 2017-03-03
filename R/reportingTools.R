## code for integrating with ReportingTools package




#' Add links to data when using ReportingTools
#' 
#' These functions are intended to add links to the Affymetrix, Entrez Gene, Nuccore,
#' and AmiGO databases when creating HTML tables using ReportingTools.
#' 
#' These functions are not actually intended to be called directly. Instead,
#' they are used as targets for the .modifyDF argument of the \code{publish}
#' function of ReportingTools. See the example below for more detail.
#' 
#' @aliases entrezLinks affyLinks goLinks nuccoreLinks ensemblLinks
#' @param df A data.frame, usually created using the \code{select} function of
#' the AnnotationDbi package. For Entrez ID data, the column name must be
#' ENTREZID. For Affy data, the column name must be PROBEID, and for GO data
#' the column name must be Term. For Nuccore the column name can be any of
#' "GI", "REFSEQ", "ACCNUM". Any other names will fail.
#' @param ... Allows one to pass arbitrary arguments to lower level functions.
#' Currently unsupported.
#' @return A data.frame is returned, with links included.
#' @author James W. MacDonald \email{jmacdon@@u.washington.edu}
#' @keywords manip
#' @examples
#' 
#' \dontrun{
#' ## say we have an ExpressionSet from HuGene 1.0 ST array
#' dat <- read.celfiles()
#' eset <- rma(dat)
#' ## annotate the ExpressionSet
#' eset <- annotateEset(eset, hugene10sttranscriptcluster.db,
#' columns = c("ENTREZID","ACCNUM","SYMBOL"))
#' ## and fit a model using limma
#' fit <- lmFit(eset, design)
#' fit2 <- eBayes(fit)
#' ## and create an HTML page with links to Affy and Entrez
#' out <- topTable(fit2, coef=2)
#' htab <- HTMLReport("The title","a_short_name")
#' publish(out, htab, .modifyDF = list(affyLinks, entrezLinks, nuccoreLinks))
#' finish(htab)}
#' 
#' @export entrezLinks affyLinks goLinks nuccoreLinks ensemblLinks
entrezLinks <- function (df, ...) {
    naind <- is.na(df$ENTREZID)
    df$ENTREZID <- hwrite(as.character(df$ENTREZID), link = paste0("http://www.ncbi.nlm.nih.gov/gene/", 
        as.character(df$ENTREZID)), table = FALSE)
    df$ENTREZID[naind] <- ""
    return(df)
}

affyLinks <- function(df, ...){
    df$PROBEID <- hwrite(as.character(df$PROBEID),
                         link = paste0("https://www.affymetrix.com/LinkServlet?probeset=", 
                         as.character(df$PROBEID)), table = FALSE)
    return(df)
}

goLinks <- function(df, ...){
    ## first column likely to be links. Extract out the GO IDs
    linkpart <- as.character(df[,1])
    if(length(grep("href", linkpart)) > 0)
        linkpart <- sapply(strsplit(linkpart, ">|<"), "[", 3)
    df$Term <- hwrite(as.character(df$Term),
                      link = paste0("http://amigo.geneontology.org/amigo/term/",
                      linkpart), table = FALSE)
    return(df)
}

nuccoreLinks <- function(df, ...){
    col <- grep("GI|REFSEQ|ACCNUM", colnames(df), ignore.case = TRUE)
    if(length(col) == 0) {
        stop(paste("There needs to be a column labeled 'GI', 'REFSEQ', or 'ACCNUM'",
                   "to add links to the nuccore database."), call. = FALSE)
    } else if(length(col) > 1) {
        stop(paste("There can only be one column labeled 'GI', 'REFSEQ', or 'ACCNUM'",
                   "to add links to the nuccore database."), call. = FALSE)
    }
    df[,col] <- hwrite(as.character(df[,col]),
                       link = paste0("http://www.ncbi.nlm.nih.gov/nuccore/",
                                     as.character(df[,col])), table = FALSE)
    return(df)
}

ensemblLinks <- function(df, ...){
    naind <- is.na(df$ENSEMBL)
    df$ENSEMBL <- hwriter::hwrite(as.character(df$ENSEMBL),
                                 link = paste0("http://www.ensembl.org/Gene/Summary?db=core;g=",
                                               as.character(df$ENSEMBL)), table = FALSE)
    df$ENSEMBL[naind] <- ""
    return(df)
}


#' Select and Output Genelists Based on Venn Diagrams
#' 
#' This function is designed to output text and/or HTML tables based on the
#' results of a call to \code{\link[limma]{decideTests}}, using the
#' ReportingTools package.
#' 
#' The purpose of this function is to output HTML and text tables with lists of
#' genes that fulfill the criteria of a call to
#' \code{\link[limma]{decideTests}} as well as the direction of differential
#' expression.
#' 
#' Some important things to note: First, the names of the HTML and text tables
#' are extracted from the \code{colnames} of the \code{TestResults} object,
#' which come from the contrasts matrix, so it is important to use something
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
#' Unlike \code{vennSelect}, this function automatically creates both HTML and
#' CSV output files.
#' 
#' @param fit An \code{\link[limma:marraylm]{MArrayLM}} object, from a call to
#' \code{\link[limma:ebayes]{eBayes}}.
#' @param contrast A contrasts matrix, produced either by hand, or by a call to
#' \code{\link[limma]{makeContrasts}}
#' @param design A design matrix.
#' @param groups This argument is used when creating a legend for the resulting
#' HTML pages. If NULL, the groups will be generated using the column names of
#' the design matrix.
#' @param cols A numeric vector indicating which columns of the fit, contrast
#' and design matrix to use. If \code{NULL}, all columns will be used.
#' @param p.value A p-value to filter the results by.
#' @param lfc A log fold change to filter the results by.
#' @param method One of "same", "both", "up", "down", "sameup", or "samedown".
#' See details for more information.
#' @param adj.meth Method to use for adjusting p-values. Default is 'BH', which
#' corresponds to 'fdr'. Ideally one would set this value to be the same as was
#' used for \code{\link[limma]{decideTests}}.
#' @param titleadd Additional text to add to the title of the HTML tables.
#' Default is NULL, in which case the title of the table will be the same as
#' the filename.
#' @param fileadd Additional text to add to the name of the HTML and CSV
#' tables. Default is NULL.
#' @param baseUrl A character string giving the location of the page in terms
#' of HTML locations. Defaults to "."
#' @param reportDirectory A character string giving the location that the
#' results will be written. Defaults to "./venns"
#' @param affy Boolean; are these Affymetrix arrays, and do you want hyperlinks
#' for each probeset to the Affy website to be generated for the resulting HTML tables?
#' @param probecol If the "affy" argument is \code{TRUE}, what is the column header
#' for the Affymetrix probeset IDs? Defaults to "PROBEID", which is the default if
#' the data are annotated using a Bioconductor annotation package.
#' @param \dots Used to pass arguments to lower level functions. 
#' @return A list with two items. First, a list of \code{HTMLReport} objects
#' from the ReportingTools package, which can be used to create an index page
#' with links to the HTML pages created by this function. See the help page for
#' HTMLReport in ReportingTools as well as the vignettes for more information.
#' The second item is a \code{vennCounts} object from limma, which can be used
#' to create a Venn diagram, e.g., in a report if this function is called
#' within a Sweave or knitR pipeline.
#' @author James W. MacDonald \email{jmacdon@@u.washington.edu}
#' @keywords manip
#' @export vennSelect2
vennSelect2 <- function(fit, contrast, design,  groups = NULL, cols = NULL, p.value = 0.05,
                        lfc = 0, method = "same", adj.meth = "BH", titleadd = NULL, fileadd = NULL, 
                        baseUrl = ".",  reportDirectory = "./venns", affy = TRUE, probecol = "PROBEID", ...){

    ## design is a design matrix from limma
    ## contrast is a conrast matrix
    ## output is a list containing the probe IDs of genes from each comparison
        
    lattice.options(default.theme = reporting.theme())  
    ## do this before we subset
    if(is.null(groups))
        groups <- factor(colnames(design)[apply(design, 1, function(x) which(x !=0))])
    else
        groups <- factor(groups)

    if(!is.null(cols)){
        contrast <- contrast[,cols]
        fit <- fit[,cols]
    }
    
    dtmat <- decideTests(fit, adjust.method = adj.meth, p.value = p.value, lfc = lfc)
    if(all(apply(dtmat, 2, sum) == 0)) stop(paste("There are no significant genes at a p.value <",
                p.value, "and a fold change >", lfc, "with",
                switch(adj.meth, none = "no multiplicity adjustment!", paste("a", adj.meth, "multiplicity adjustment!"))),
                call. = FALSE)
    colind <- getCols(design, contrast)
    
    
    ncontrasts <- ncol(dtmat)
    if(ncontrasts < 2 || ncontrasts > 4)
        stop("This function only works for 2-4 comparisons at a time.\n",
             call. = FALSE)
    if(ncontrasts == 2)
        name <- c(paste("Genes unique to", colnames(dtmat)),
                  paste("Genes in intersection of", intNames(dtmat)))
    if(ncontrasts == 3)
        name <- c(paste("Genes unique to", colnames(dtmat)),
                  paste("Genes in intersection of", intNames(dtmat)),
                  "Genes common to all comparisons")
    if(ncontrasts == 4)
        return(venn4Way(fit = fit, contrast = contrast, p.value = p.value,
                        lfc = lfc, adj.meth = adj.meth, baseUrl = baseUrl, reportDirectory = reportDirectory))
        
    
    coefind <- switch(ncol(dtmat)-1,
                      list(1,2,1:2),
                      list(1,2,3,1:2,c(1,3), 2:3, 1:3))
    
    ## Remove illegal characters from filenames
    if(length(grep("[/|\\|?|*|:|<|>|\"|\\|]", name)) > 0)
        warning(paste("Some illegal characters have been removed from the filenames",
                      name, sep = " "), call. = FALSE)
    name <- gsub("[/|\\|?|*|:|<|>|\"|\\|]", "", name)

    if(!is.null(titleadd))
        titles <- paste(name, titleadd)
    else
        titles <- name

    if(!is.null(fileadd))
        name <- paste(name, fileadd)
    
    
    indices <- makeIndices(dtmat, method = method)
    ind <- which(sapply(indices, sum) > 0)
    venncounts <- vennCounts2(dtmat, method = method)
    
    vennout <- lapply(seq(along = indices), function(x) HTMLReport(gsub(" ", "_", name[x]), titles[x],
                          baseUrl = baseUrl,reportDirectory = reportDirectory))
    ## require(annotation(eset), character.only = TRUE, quietly = TRUE)
    ps <- apply(fit$p.value, 2, p.adjust, method = adj.meth)
    colnames(ps) <- paste0(colnames(ps), ".p.value")
    coefs <- fit$coefficients
    colnames(coefs) <- paste0(colnames(coefs), ".logFC")
    csvlst <- lapply(seq(along = indices), function(x) cbind(coefs[indices[[x]], coefind[[x]], drop = FALSE],
                         ps[indices[[x]], coefind[[x]], drop = FALSE]))
    ## previously we used the ExpressionSet only to figure out the annotation source.
    ## instead of that, we now rely upon the MArrayLM object containing the annotations we want to use
    rn <- row.names(fit$coef)
    if(is.null(rn)) rn <- seq_len(nrow(fit$coef))
    if(is.null(fit$genes)){
        warning(paste("\nThere are no annotations in the MArrayLM object, using",
                      if(is.null(row.names(fit$coef))) "row indices" else "row names",
                      "from the MArrayLM object instead"))
        annotlst <- lapply(indices, function(x) if(sum(x) > 0) data.frame(ID = rn[x]))
    } else {
        annotlst <- lapply(indices, function(x) if(sum(x) > 0) fit$genes[x,])
    }
                       
    csvlst <- mapply(data.frame, annotlst, csvlst, SIMPLIFY = FALSE)
    fixed.dflst <- lapply(csvlst[ind], fixHeaderAndGo, affy = affy, probecol = probecol)
    csvlst[ind] <- lapply(fixed.dflst, function(x) x$df)
    mdf <- fixed.dflst[[1]]$mdf
    
    if(length(mdf) > 0)
        lapply(ind, function(x) publish(csvlst[[x]], vennout[[x]],  .modifyDF = mdf))
    else
        lapply(ind, function(x) publish(csvlst[[x]], vennout[[x]]))
    
    vennpaths <- gsub("html$", "txt", sapply(vennout, path))
    lapply(ind, function(x) write.table(csvlst[[x]], vennpaths[x],
                                        sep = "\t", quote = FALSE, row.names = FALSE, na = ""))
    
    lapply(vennout[ind], finish)
    return(list(vennout = vennout, venncounts = venncounts))
}

           
    
## makeLegend <- function(groups, dir, fname){
##     nam <- paste(dir, fname, sep = "/")
##     png(nam, height = 250, width = 300)
##     plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
##     cols <- lattice.options()$default.theme$superpose.symbol$col
##     legend(1,1, levels(groups), pch = 16, col = cols[sort(unique(as.numeric(groups)))],
##            xjust = 0.5, yjust=0.5, cex = 1.5, bty="n")
##     dev.off()
##     fname
## }



##' This function is designed to output CSV and HTML tables based on an analysis
##' using the limma or edgeR packages, with output generated using the ReportingTools
##' package. Please note that a DGEGLM object from edgeR is simply converted to an
##' MArrayLM object from limma and then used in the default MArrayLM method, so all
##' arguments for the MArrayLM object pertain to the DGEGLM method as well.
##' 
##' The purpose of this function is to output HTML and text tables with lists of
##' genes that fulfill the criteria of a call to
##' \code{\link[limma]{decideTests}} as well as the direction of differential
##' expression. This is a high-level function that calls \code{vennSelect2}
##' internally, and is intended to be used with \code{vennPage} to create a set
##' of Venn diagrams (on an HTML page) that have clickable links in each cell of
##' the diagram. The links will then pass the end user to individual HTML pages
##' that contain the genes that are represented by the counts in a given cell of
##' the Venn diagram.
##' 
##' In general, the only thing that is needed to create a set of Venn diagrams
##' is a list of numeric vectors that indicate the columns of the contrast
##' matrix that are to be used for a given diagram. See the example below for a
##' better explanation.
##' 
##' Some important things to note: First, the names of the HTML and text tables
##' are extracted from the \code{colnames} of the \code{TestResults} object,
##' which come from the contrasts matrix, so it is important to use something
##' descriptive. Second, the method argument is analogous to the \code{include}
##' argument from \code{\link[limma:venn]{vennCounts}} or
##' \code{\link[limma:venn]{vennDiagram}}. Choosing "both" will select genes
##' that are differentially expressed in one or more comparisons, regardless of
##' direction. Choosing "up" or "down" will select genes that are only
##' differentially expressed in one direction. Choosing "same" will select genes
##' that are differentially expressed in the same direction. Choosing "sameup"
##' or "samedown" will select genes that are differentially expressed in the
##' same direction as well as 'up' or 'down'.
##' 
##' Note that this is different than sequentially choosing "up" and then "down".
##' For instance, a gene that is upregulated in one comparison and downregulated
##' in another comparison will be listed in the intersection of those two
##' comparisons if "both" is chosen, it will be listed in only one comparison
##' for both the "up" and "down" methods, and it will be listed in the union
##' (e.g., not selected) if "same" is chosen.
##' 
##' Unlike \code{vennSelect}, this function automatically creates both HTML and
##' CSV output files.
##'
##' Also please note that this function relys on annotation information contained in
##' the "genes" slot of the "fit" object. If there are no annotation data, then
##' just statistics will be output in the resulting HTML tables.
##'
##' @title High-level function for making Venn diagrams and outputting the results from
##' the diagrams in HTML and CSV files.
##' @param object An \code{\link[limma:marraylm]{MArrayLM}} or \code{\link[edgeR:DGEGLM-class]{DGEGLM}} object.
##' @param contrast A contrasts matrix, produced either by hand, or by a call to
##' \code{\link[limma]{makeContrasts}}
##' @param design A design matrix.
##' @param groups This argument is used when creating a legend for the resulting
##' HTML pages. If NULL, the groups will be generated using the column names of
##' the design matrix. In general it is best to leave this NULL.
##' @param collist A list containing numeric vectors indicating which columns of
##' the fit, contrast and design matrix to use. If \code{NULL}, all columns will
##' be used.
##' @param p.value A p-value to filter the results by.
##' @param lfc A log fold change to filter the results by.
##' @param method One of "same", "both", "up", "down", "sameup", or "samedown".
##' See details for more information.
##' @param adj.meth Method to use for adjusting p-values. Default is 'BH', which
##' corresponds to 'fdr'. Ideally one would set this value to be the same as was
##' used for \code{\link[limma]{decideTests}}.
##' @param titleadd Additional text to add to the title of the HTML tables.
##' Default is NULL, in which case the title of the table will be the same as
##' the filename.
##' @param fileadd Additional text to add to the name of the HTML and CSV
##' tables. Default is NULL.
##' @param baseUrl A character string giving the location of the page in terms
##' of HTML locations. Defaults to "."
##' @param reportDirectory A character string giving the location that the
##' results will be written. Defaults to "./venns"
##' @param affy Boolean. Are these Affymetrix data, and should hyperlinks to the affy website
##' be generated in the HTML tables?
##' @param probecol This argument is used in concert with the preceding argument. If these are Affymetrix data
##' , then specify the column header in the \code{\link[limma:marraylm]{MArrayLM}} object that contains the Affymetrix IDs. Defaults to
##' "PROBEID", which is the expected result if the data are annotated using a BioC annotation package.
##' @param ... Used to pass other arguments to lower level functions.
##' @return A list containing the output from calling \code{vennSelect2} on the
##' columns specified by the collist argument. This is intended as input to
##' \code{vennPage}, which will use those data to create the HTML page with Venn
##' diagrams with clickable links.
##' @author James W. MacDonald \email{jmacdon@@u.washington.edu}
##' @keywords manip
##' @examples
##' 
##'   \dontrun{
##'     mat <- matrix(rnorm(1e6), ncol = 20)
##'     design <- model.matrix(~factor(1:4, each=5))
##'     colnames(design) <- LETTERS[1:4]
##'     contrast <- matrix(c(1,-1,0,0,1,0,-1,0,1,0,0,-1,0,1,-1,0,0,1,0,-1),
##'     ncol = 5)
##'     colnames(contrast) <- paste(LETTERS[c(1,1,1,2,2)],
##'     LETTERS[c(2,3,4,3,4)], sep = " vs ")
##'     fit <- lmFit(mat, design)
##'     fit2 <- contrasts.fit(fit, contrast)
##'     fit2 <- eBayes(fit2)
##'     ## two Venn diagrams - a 3-way Venn with the first three contrasts
##'     ## and a 2-way Venn with the last two contrasts
##'     collist <- list(1:3,4:5)
##'     venn <- makeVenn(fit2, contrast, design, collist = collist)
##'     vennPage(venn, "index.html", "Venn diagrams")
##'     }
##' 
##' @describeIn makeVenn Make a Venn diagram using an MArrayLM object.
##' @export makeVenn
setMethod("makeVenn", "MArrayLM",
          function(object, contrast, design, groups = NULL, collist = NULL,
                   p.value = 0.05, lfc = 0, method = "both", adj.meth = "BH",
                   titleadd = NULL, fileadd = NULL, baseUrl = ".", reportDirectory = "./venns",
                   affy = TRUE, probecol = "PROBEID", ...){
    if(is.null(collist)){
        if(ncol(object$coef) > 4) stop("You can only make Venn diagrams with four or fewer contrasts!\n\n",
                                    call. = FALSE)
        collist <- list(seq_len(ncol(object$coef)))
    }
    vennlst <- lapply(seq(along = collist), function(x) 
        vennSelect2(fit = object, contrast = contrast, design = design, 
                    groups = groups, cols = collist[[x]], p.value = p.value, lfc = lfc,
                    method = method, adj.meth = adj.meth, titleadd = titleadd,
                    fileadd = fileadd, baseUrl = baseUrl, 
                    reportDirectory = paste0(reportDirectory, "/venn", x),
                    affy = affy, probecol = probecol, ...))
    vennlst
})


##' @param comp.method Character. For DGEGLM objects, the DGEGLM object must first be processed using one of \code{\link[edgeR:glmfit]{glmLRT}},
##' \code{\link[edgeR]{glmQLFTest}}, or \code{\link[edgeR]{glmTreat}}. Choose glmLRT if you fit a model using
##' \code{\link[edgeR:glmfit]{glmFit}}, glmQLFTest if you fit a model using \code{\link[edgeR]{glmQLFit}}, or glmTreat if
##' you fit either of those models, but want to incorporate the log fold change into the comparison.
##'  
##' 
##' @describeIn makeVenn Make a Venn diagram using a DGEGLM object.
##' @export
setMethod("makeVenn", "DGEGLM",
          function(object, contrast, design, comp.method = c("glmLRT","glmQLFTest", "glmTreat"), lfc = 0, ...){
    comp.method <-  match.arg(comp.method, c("glmLRT","glmQLFTest", "glmTreat"))
    if(comp.method == "glmTreat" && lfc <= 0)
        stop(paste("When using the glmTreat method you must also choose a non-zero lfc value."), call. = FALSE)
    lst <- switch(comp.method,
                  glmLRT = lapply(seq_len(ncol(contrast)), function(x) glmLRT(object, contrast = contrast[,x])),
                  glmQLFTest = lapply(seq_len(ncol(contrast)), function(x) glmQLFTest(object, contrast = contrast[,x])),
                  glmTreat = lapply(seq_len(ncol(contrast)), function(x) glmTreat(object, contrast = contrast[,x], lfc = lfc)))
    object <- new("MArrayLM", list(p.value = do.call(cbind, lapply(lst, function(x) x$table$PValue)),
                                   coefficients = do.call(cbind, lapply(lst, function(x) x$table$logFC)),
                                   genes = object$genes))
    colnames(object$p.value) <- colnames(object$coefficients) <- colnames(contrast)
    makeVenn(object, contrast, design, lfc = lfc, ...)
})
   


##' High-level function for making Venn diagrams with clickable links to HTML
##' pages with the underlying genes.
##' 
##' This function is designed to be used in conjunction with the \code{makeVenn}
##' function, to first create a set of HTML pages containing the genes that are
##' represented by the cells of a Venn diagram, and then create an HTML page
##' with the same Venn diagrams, with clickable links that will point the end
##' user to the HTML pages.
##' 
##' This function is intended to be used as part of a pipeline, by first calling
##' \code{makeVenn} and then using the output from that function as input to
##' this function to create the HTML page with clickable links.
##' 
##' @param vennlst The output from \code{makeVenn}.
##' @param pagename Character. The file name for the resulting HTML page.
##' Something like 'venns' is reasonable. Note that the .html will automatically
##' be appended.
##' @param pagetitle Character. The heading for the HTML page.
##' @param cex.venn Numeric. Adjusts the size of the font in the Venn diagram.
##' Usually the default is OK.
##' @param shift.title Boolean. Should the right contrast name of the Venn
##' diagram be shifted down? Useful for long contrast names. If a two-way Venn
##' diagram, this will shift the right name down so they don't overlap. If a
##' three-way Venn diagram, this will shift the top right name down.
##' @param baseUrl Character. The base URL for the resulting HTML page. The
##' default of "." is usually optimal.
##' @param reportDirectory If \code{NULL}, the reportDirectory will be extracted
##' from the vennlst. This is usually what one should do.
##' @param ... To allow passing other arguments to lower level functions.
##' Currently not used.
##' @return An HTMLReport object. If used as input to the ReportingTools
##' \code{publish} function, this will create a link on an index page to the
##' Venn diagram HTML page. See e.g., the microarray analysis vignette for
##' ReportingTools for more information.
##' @author James W. MacDonald \email{jmacdon@@u.washington.edu}
##' @keywords manip
##' @examples
##' 
##'   \dontrun{
##'     mat <- matrix(rnorm(1e6), ncol = 20)
##'     design <- model.matrix(~factor(1:4, each=5))
##'     colnames(design) <- LETTERS[1:4]
##'     contrast <- matrix(c(1,-1,0,0,1,0,-1,0,1,0,0,-1,0,1,-1,0,0,1,0,-1),
##'     ncol = 5)
##'     colnames(contrast) <- paste(LETTERS[c(1,1,1,2,2)],
##'     LETTERS[c(2,3,4,3,4)], sep = " vs ")
##'     fit <- lmFit(mat, design)
##'     fit2 <- contrasts.fit(fit, contrast)
##'     fit2 <- eBayes(fit2)
##'     ## two Venn diagrams - a 3-way Venn with the first three contrasts
##'     ## and a 2-way Venn with the last two contrasts
##'     collist <- list(1:3,4:5)
##'     venn <- makeVenn(fit2, contrast, design, eset, collist = collist)
##'     vennreport <- vennPage(venn, "index.html", "Venn diagrams")
##'     indexPage <- HTMLReport("index", "My results", reportDirectory =
##'     ".", baseUrl = ".")
##'     publish(vennreport)
##'     finish(indexPage)
##'     }
##' 
##' @export vennPage

vennPage <- function(vennlst, pagename, pagetitle, cex.venn = 1, shift.title = FALSE,
                     baseUrl = ".", reportDirectory = NULL, ...){
    if(is.null(reportDirectory)){
        tmp <- strsplit(path(vennlst[[1]]$vennout[[1]]), "/")[[1]]
        reportDirectory <- paste(tmp[1:2], collapse = "/")
    }
    hpage <- openPage(paste0(pagename, ".html"), dirname = reportDirectory)
    hwrite(paste("The Venn diagrams all contain clickable links. Click on the counts",
                          "in any cell to see a table of the genes in that cell.",
                          "Also please note that the tables are sortable - simply click",
                          "on any header to sort on that column."), hpage,
                    br = TRUE)
                        
    lapply(seq(along = vennlst), function(x) drawVenn(vennlst[[x]], page = hpage, 
               dir = reportDirectory, num = x, cex = cex.venn, shift.title = shift.title, ...))
    closePage(hpage)
    paste0(reportDirectory, "/", pagename, ".html")
}

drawVenn <- function(lst, page, dir, num, cex = 1, shift.title = FALSE, ...){
    nam <- paste0(dir, "/venn", num, ".png")
    nam2 <- paste0("venn", num, ".png")
    mapname <- paste0("#venn", num)
    png(nam, height = 800, width = 800)
    if(is(lst$venncounts, "VennCounts")){
        if(shift.title) colnames(lst$venncounts)[1] <- paste0(colnames(lst$venncounts)[1], "\n\n")
        vennDiagram(lst$venncounts, cex = cex)
    } else if(is(lst$venncounts, "venn")) {
        plot(lst$venncounts)
        addLines(lst$nam)
    } else {
        stop(paste("You cannot create a Venn diagram with an object of class",
                      class(lst$venncounts)), call. = FALSE)
    }
    dev.off()
    
    hwrite(hmakeTag("img", border = 0, width = 800, height = 800,
                    src = nam2, alt = nam2, usemap = mapname), page)
    hwriteImage(paste("Venn Diagram", num), page)
    vennLinks(lst, page, mapname, dir)
    
}

vennLinks <- function(lst, page, mapname, reportDirectory){
    if(is.null(page)) page <- ""
    fun <- function(x,y) paste0('<area shape="circle" coords=',
                                    x, ' href=', y, '></area>')
    if(is(lst$venncounts, "VennCounts")){
        if(ncol(lst$venncounts) == 3){
            loclst <- list("250,400,30","550,400,30","400,400,30")
        } else {
            loclst <- list("260,310,30","540,300,30","400,550,30",
                           "400,320,30","330,420,30","470,420,30",
                           "400,390,30")
        }
    } else {
        loclst <- list("140,315,30", "315,215,30","510,215,30","685,325,30","230,270,30",
                       "230,530,30","415,625,30","415,255,30","575,530,30","595,270,30",
                       "300,345,30","485,585,30","340,585,30","530,345,30","415,470,30")
    }
    urlst <- gsub(reportDirectory, ".", sapply(lst$vennout, path))
    strng <- do.call("c", mapply(fun, loclst, urlst, SIMPLIFY = FALSE))
    strng <- c(paste0('<map name="', sub("#", "", mapname), '">'),
               strng, "</map>")
    cat(strng, file = page, sep = "\n")
}


## code for 4-way venn diagrams with colored lines under the counts

addLines <- function(comps, col = c("#66FF33","#FF9933", "#000000","#9900CC")){
    x <- rep(c(35,140,260,365,90,95,200,200,300,310,130,245,155,270,200),
             c(1,1,1,1,2,2,2,2,2,2,3,3,3,3,4))
    y <- rep(c(250,315,315,250,280,110,50,290,110,280,230,95,95,230,150),
             c(1,1,1,1,2,2,2,2,2,2,3,3,3,3,4))
    x0 <- x - 10
    x1 <- x + 10
    ind <- c(rep(1, 4), rep(1:2, times=6), rep(1:3, times=4), 1:4)
    y <- y - c(7.5, 10, 12.5, 15)[ind]
    ind2 <- c(1,2,3,4,1,2,1,3,1,4,2,3,2,4,3,4,1,2,3,1,2,4,1,3,4,2,3,4,1,2,3,4)
    segments(x0, y, x1, y, col = col[ind2], lwd = 2)
    legend("topleft", comps, lty = 1, lwd = 2, col = col, bty = "n",
           cex = 0.8)
}
##' A function to create a 4-way Venn diagram
##'
##' This function is an internal function and not really intended to be called by the end user. It is generally called by the \code{vennPage}
##' function. The goal is to create a 4-way Venn diagram in an HTML page with clickable links to tables of the genes found in a given cell.
##' In addition, the numbers in each cell are underlined with colored bars that help end users tell what contrasts are captured by that cell.
##' @title 4-way Venn Diagrams
##' @param fit An \code{MArrayLM} object, created by the limma package.
##' @param contrast A contrasts matrix, used by limma to generate the comparisons made.
##' @param p.value A p-value cutoff for significance
##' @param lfc A log fold change cutoff
##' @param adj.meth The method used to adjust for multiple comparisons.
##' @param baseUrl The base directory for the tables generated. Defaults to ".", meaning the current directory. 
##' @param reportDirectory The directory in which to put the results. Defaults to a "venns" subdirectory.
##' @param ... Allows arbitrary arguments to be passed to lower level functions
##' @return Returns a list. The first item is a (list of) HTMLReportRef objects that can be used by ReportingTools to create HTML links.
##' The second item is the output from the \code{venn} function in gtools, and the third item is the name of the contrasts used to generate
##' the Venn diagram.
##' @author James W. MacDonald \email{jmacdon@@u.washington.edu}
venn4Way <- function(fit, contrast,  p.value, lfc, adj.meth, baseUrl = ".", reportDirectory = "./venns", ...){
    ##require("gplots", character.only = TRUE) || stop("The gplots package is required for 4-way Venns.\n", call. = FALSE) 
    nam <- paste0(gsub(" ", ".", colnames(contrast)), ".")
    ## do some munging to get the right columns
    allcols <- c("PROBEID","ID","SYMBOL","GENENAME","logFC","t",switch(adj.meth, none = "P.Value", "adj.P.Value"))
    shortcols <- c("logFC","t",switch(adj.meth, none = "P.Value", "adj.P.Value"))
    outlst <- lapply(1:ncol(contrast), function(x) tableFilt(fit, x, pfilt = p.value, fldfilt = lfc, adjust = adj.meth))
    allcols <- match(allcols[allcols %in% names(outlst[[1]])], names(outlst[[1]]))
    shortcols <- match(shortcols[shortcols %in% names(outlst[[1]])], names(outlst[[1]]))
    theplot <- venn(lapply(outlst, function(x) x[,1]), show.plot = FALSE)
    for(i in seq_len(length(nam))) names(outlst[[i]]) <- paste0(rep(c("", nam[i]), c(4,3)), names(outlst[[i]]))
    vcind <- c(1,4,9,16)
    vennlst2 <- do.call("rbind", lapply(1:4, function(x) data.frame(PROBEID = outlst[[x]]$PROBEID, grp = vcind[x])))
    vennlst2 <- sapply(tapply(1:nrow(vennlst2), vennlst2$PROBEID, function(x) vennlst2[x,2]), sum)
    vennlst3 <- lapply(c(1,4,5,9,10,13,14,16,17,20,21,25,26,29,30), function(x) names(vennlst2[vennlst2 == x]))
    indmat <- rbind(diag(4),cbind(1, diag(3)), cbind(0,1,diag(2)), c(0,0,1,1))
    z <- matrix(1, 4, 4)
    diag(z) <- 0
    indmat <- rbind(indmat, z[4:1,], rep(1,4))
    indmat <- matrix(as.logical(indmat), ncol = 4)
    theind <- c(1,2,4,8,3,5,9,6,10,12,7,11,13,14,15)
    fn <- LETTERS[1:4]
    out <- lapply(seq_len(length(theind)), function(x) {
        if(sum(indmat[x,]) == 1) {
            return(outlst[[which(indmat[x,])]][outlst[[which(indmat[x,])]]$PROBEID %in% vennlst3[[theind[x]]],allcols])
        }else{
            littleind <- which(indmat[x,])
            tmp1 <- outlst[[littleind[1]]][match(vennlst3[[theind[x]]], outlst[[littleind[1]]]$PROBEID),allcols]
            tmp2 <- do.call("cbind", lapply(littleind[-1], function(y)
                                            outlst[[y]][match(vennlst3[[theind[x]]], outlst[[y]]$PROBEID),shortcols]))
            return(cbind(tmp1, tmp2))
        }
    })
    ind <- which(sapply(out, nrow) > 0)
    nam <- sapply(seq_len(length(out)), function(x) paste0(paste(fn[indmat[x,]], collapse = "_"), "_venn"))
    vennout <- lapply(seq_len(length(out)), function(x) HTMLReport(nam[x], nam[x], baseUrl = baseUrl, reportDirectory = reportDirectory))
    lapply(ind, function(x) publish(out[[x]], vennout[[x]], .modifyDF = list(affyLinks)))
    lapply(vennout[ind], finish)
    lapply(ind, function(x) write.table(out[[x]], paste0(reportDirectory, "/", nam[x], ".txt"), sep = "\t",
                                        quote = FALSE, row.names = FALSE, na = ""))
    return(list(vennout = vennout, venncounts = theplot, nam = colnames(contrast)))
}



##' A function to add dotplot glyphs and links to HTML tables
##'
##' This function is intended to create little dotplot glyphs that can be added to an HTML table of
##' results from e.g., a microarray or RNA-Seq experiment, showing graphically how much the different groups
##' are changing. The glyphs have unlabeled axes to make them small enough to fit in an HTML table, and clicking
##' on a glyph will result in a new page loading with a full sized dotplot, complete with axis labels.
##'
##' This function is very similar to the stock functions in the ReportingTools package, but the standard
##' glyphs for that package consist of a dotplot on top of a boxplot, which seems too busy to me. In addition,
##' for most microarray analyses there are not enough replicates to make a boxplot useful.
##' @title Add dotplot images
##' @param df A data.frame from calling \code{topTable}. Note that the row.names for this data.frame must
##' be consistent with the "eset" object. In other words, if "eset" is an \code{ExpressionSet}, then the row.names
##' of the data.frame must consistent with the featureNames of the \code{ExpressionSet}. 
##' @param eset A matrix, data.frame, or \code{ExpressionSet}. If using RNA-Seq data, use \code{voom} from edgeR to create
##' an \code{EList} object, and then pass in the "E" list item.
##' @param grp.factor A factor that indicates which group ALL of the samples belong to. This will be subsetted internally,
##' so do not subset yourself.
##' @param design The design matrix used by limma or edgeR to fit the model.
##' @param contrast The contrast matrix used by limma or edgeR to make comparisons.
##' @param colind Which column of the contrast matrix are we using? In other words, for which comparison are we creating a table?
##' @param boxplot Boolean. If \code{TRUE}, the output HTML table will have a boxplot showing differences between groups. If
##' \code{FALSE} (default), the table will have dotplots. 
##' @param repdir A directory in which to put the HTML tables. Defaults to a "reports" directory in the working directory.
##' @param extraname By default, the tables will go in a "reports" subdirectory, and will be named based on the column
##' name of the contrast that is specified by the colind argument (after replacing any spaces with an underscore). If this
##' will result in name collisions (e.g., a previous file will be over-written because the resulting names are the same),
##' then an extraname can be appended to ensure uniqueness.
##' @param weights Array weights, generally from \code{arrayWeights} in the limma package. These will affect the size
##' of the plotting symbols, to reflect the relative importance of each sample.
##' @param insert.after Which column should the image be inserted after? Defaults to 3.
##' @param \dots Allows arbitrary arguments to be passed down to lower level functions.
##' @return A list, two items. The first item is the input data.frame with the glyphs included, ready to be used with
##' ReportingTools to create an HTML table. The second item is a pdf of the most differentially expressed comparison. This is
##' useful for those who are using e.g., knitr or Sweave and want to be able to automatically insert an example dotplot
##' in the document to show clients what to expect.
##' @author James W. MacDonald \email{jmacdon@@u.washington.edu}
##' @export makeImages
makeImages <- function(df, eset, grp.factor, design, contrast, colind, boxplot = FALSE, repdir = "./reports", extraname = NULL, weights = NULL,
                       insert.after = 3, ...){
    ## check that this is going to work
    eclass <- class(eset)[1]
    addtrailingslash <- function(path){
        tmp <- strsplit(repdir, "")[[1]]
        if(!tmp[length(tmp)] %in% "/") path <- paste0(path, "/")
        path
    }
    repdir <- addtrailingslash(repdir)
    rn <- switch(eclass,
                 ExpressionSet = featureNames(eset),
                 matrix = row.names(eset),
                 data.frame = row.names(eset),
                 DGEList = row.names(eset$counts),
                 stop(paste0("The 'eset' argument is of class", eclass, ". This function can only use",
                            "an ExpressionSet, matrix, data.frame or DGEList object"), call. = FALSE))
    
    if(!all(row.names(df) %in% rn))
        stop(paste0("The row.names of your input data.frame do not match up with the data in your", eclass, ".",
                    " You need to fix that before proceeding!\n"), call. = FALSE)
    figure.directory <- paste0(repdir, gsub(" ", "_", colnames(contrast)[colind]))
    if(!is.null(extraname)) figure.directory <- paste0(figure.directory, extraname)
    dir.create(figure.directory, recursive = TRUE)
    ind <- apply(design[,contrast[,colind] != 0, drop = FALSE], 1, sum) > 0
    grp.factor <- factor(grp.factor[ind])
    eset <- eset[,ind]
    colnames(df)[colnames(df) == "SYMBOL"] <- "Symbol"
    makeGenePlots(df = df, expression.dat = eset, factor = grp.factor, figure.directory = figure.directory,
                  boxplot = boxplot, weights = weights, ...)
   
    figure.directory <- gsub(repdir, "", figure.directory)
    mini.image <- file.path(figure.directory, paste("mini", 
                                                   rownames(df), "png", sep = "."))
    pdf.image <- file.path(figure.directory, paste("boxplot", 
                                                  rownames(df), "pdf", sep = "."))
    df <- data.frame(df[,1:insert.after, drop = FALSE], Image = hwriteImage(mini.image, link = pdf.image, table = FALSE),
                     df[,(insert.after + 1):ncol(df), drop = FALSE])
    return(list(df = df, top.pdf = pdf.image[1]))
}

## this I just stole from ReportingTools because I don't like their stupid dotplots on top of boxplots

makeGenePlots <- function (df, expression.dat, factor, figure.directory, boxplot, ylab.type = "Expression Value", 
    scales = list(), par.settings = list(), xlab = NULL, weights = NULL, ...) {
    scales <- c(scales, list(x = list(rot = 45)))
    eclass <- class(expression.dat)[1]
    expression.dat <- switch(eclass,
                             ExpressionSet = exprs(expression.dat),
                             matrix = expression.dat,
                             data.frame = as.matrix(expression.dat),
                             DGEList = cpm(expression.dat, log = TRUE),
                             stop(paste0("The 'expression.dat' argument is of class", eclass, ". This function can only use",
                            "an ExpressionSet, matrix, data.frame or DGEList object"), call. = FALSE))
        
    if (any(!rownames(df) %in% rownames(expression.dat))) {
        stop(paste("Can't find expression data for some features\n"))
    }
    for (probe in rownames(df)) {
        if ("Symbol" %in% colnames(df)) {
            ylab <- paste(df[probe, "Symbol"], ylab.type)
        }
        else {
            ylab <- paste(probe, ylab.type)
        }
        if(boxplot) {
            bigplot <- bwplot(expression.dat[probe, ] ~ factor, 
                              groups = factor, ylab = ylab, 
                              scales = scales, par.settings = par.settings, xlab = xlab)
        } else {
            if(is.null(weights)) {
                bigplot <- dotplot(expression.dat[probe, ] ~ factor, 
                                   groups = factor, ylab = ylab, 
                                   scales = scales, par.settings = par.settings, xlab = xlab)
            } else {
                bigplot <- dotplot(expression.dat[probe, ] ~ factor, 
                                   ylab = ylab, scales = scales, 
                                   par.settings = par.settings, xlab = xlab,
                                   col = trellis.par.get()$superpose.symbol$col[-4], col.var = factor,
                                   cex = weights,
                                   panel = function(x, y, cex, col, col.var, subscripts, ...) {
                                       panel.dotplot(x, y, col = col[col.var[subscripts]],
                                                     cex = cex[subscripts], ...)})
            }
        }
        miniplot <- remove.axis.and.padding(bigplot)
        minipng.filename <- paste("mini", probe, "png", sep = ".")
        minipng.file <- file.path(figure.directory, minipng.filename)
        png(minipng.file, height = 40, width = 200)
        grid.newpage()
        pushViewport(viewport(angle = 270, height = unit(220, 
            "points"), width = unit(44, "points"), name = "VP"))
        print(miniplot, newpage = FALSE)
        upViewport()
        dev.off()
        pdf.filename <- paste("boxplot", probe, "pdf", sep = ".")
        pdf.file <- file.path(figure.directory, pdf.filename)
        pdf(pdf.file, height = 4.5, width = 4.5)
        print(bigplot)
        dev.off()
    }
}

## this stolen from ReportingTools as well
remove.axis.and.padding <- function(plot) {
    plot$par.settings$layout.heights <- list(top.padding = 0, 
        main.key.padding = 0, key.axis.padding = 0, axis.xlab.padding = 0, 
        xlab.key.padding = 0, key.sub.padding = 0, bottom.padding = 0)
    plot$par.settings$layout.widths <- list(left.padding = 0, 
        key.ylab.padding = 0, ylab.axis.padding = 0, axis.key.padding = 0, 
        right.padding = 0)
    plot$x.scales$draw <- FALSE
    plot$y.scales$draw <- FALSE
    plot$xlab <- NULL
    plot$ylab <- NULL
    return(plot)
}
##' A function to create an HTML table showing genes that gave rise to a significant GO term
##'
##' This is an internal function, not intended to be called by the end user. Documentation here for clarity.
##' After running a GO analysis, it is advantageous to output a table listing those
##' genes that gave rise to a significant GO term. This function creates the table, along with links to Netaffx
##' (if the data are Affymetrix) and to the NCBI Gene database (if there are Entrez Gene IDs).
##' @title Make Gene table from GO analysis results
##' @param fit.table The output from \link[limma]{topTable}
##' @param probe.sum.table The output from running \link[GOstats]{probeSetSummary} on a \link[GOstats:GOHyperGResult-class]{GOHyperGResults} object.
##' @param go.id The GO ID of interest
##' @param cont.name The contrast name.
##' @param base.dir Character. Where should the HTML tables be generated? Defaults to NULL. 
##' @param extraname Character. An extra name that can be used if the contrast name isn't descriptive enough.
##' @param probecol The column name in the topTable object that contains probe IDs. Defaults to PROBEID.
##' @param affy Boolean. Are the arrays from Affymetrix?
##' @return Returns an \link[ReportingTools:HTMLReportRef-class]{HTMLReportRef} object.
##' @author Jim MacDonald
makeGoGeneTable <- function(fit.table, probe.sum.table, go.id, cont.name, base.dir = NULL, extraname = NULL, 
                        probecol = "PROBEID", affy = TRUE){
    ## Windows doesn't like ':' in a file name and will silently replace with '_' anyway
    go.id <- gsub(":", "_", go.id)
    prbs <- probe.sum.table$ProbeSetID[as.logical(probe.sum.table$selected)]
    out <- fit.table[fit.table[,probecol] %in% prbs, , drop = FALSE]
    rep.dir <- gsub(" ", "_", cont.name)
    if(!is.null(extraname)) rep.dir <- paste0(rep.dir, "/", gsub(" ", "_", extraname))
    fixed <- fixHeaderAndGo(out, affy = affy, probecol = probecol)
    htab <- HTMLReport(go.id, paste0("Genes responsible for significance of ", go.id, " in contrast ", cont.name,
                                         if(!is.null(extraname)) paste(",", extraname)),
                       reportDirectory = rep.dir, basePath = base.dir)
    publish(fixed$df, htab, if(length(fixed$mdf) > 0) .modifyDF = fixed$mdf)
    finish(htab)
    htab
}
##' Internal function used to automatically test for columns that can be converted to links
##'
##' This is an internal function designed to test for the presence of Affymetrix Probeset IDs or
##' Entrez Gene IDs, and if found, generate a list that can be passed to the ReportingTools publish
##' function in order to generate hyperlinks. The underlying assumption is that the data will have been
##' annotated using a Bioconductor annotation package, and thus Affy probeset IDs will have a column header
##' "PROBEID", and Entrez Gene IDs will have a header "ENTREZID" (or any combination of upper and lowercase letters).
##' @title Fix data.frame header for use with ReportingTools
##' @param df A data.frame
##' @param affy Boolean; does the data.frame contain Affymetrix probeset IDs?
##' @param probecol Character. The column header containing Affymetrix probeset IDs. Defaults to "PROBEID".
##' @return Returns a list of length two (with names mdf and df). The mdf object can be passed to the
##' \code{\link[ReportingTools:publish-methods]{publish}} using the .modifyDF argument, and the df object is
##' input dat.frame with column names corrected to conform to \code{affyLinks} and \code{entrezLinks}, so links
##' will be generated correctly.
##' @author Jim MacDonald
fixHeaderAndGo <- function(df, affy = TRUE, probecol = "PROBEID"){
     any.entrez <- grep("entrezid", names(df), ignore.case = TRUE)
     if(length(any.entrez) > 0){
         names(df)[any.entrez] <- "ENTREZID"
         any.entrez <- TRUE
     }
     if(affy){
         any.affy <- grep(probecol, names(df), ignore.case = TRUE)
         if(length(any.affy) > 0) 
             names(df)[any.affy] <- "PROBEID"
         else
             stop(paste("\n\nFailed trying to generate links to the Affymetrix",
                        "website. If these are not Affy data, please set affy = FALSE",
                        "otherwise set the probecol argument to match the column",
                        "containing the Affymetrix Probeset IDs.\n"), call. = FALSE)
     }
     mdf <- list(affyLinks, entrezLinks)[c(affy, any.entrez)]
     return(list(mdf = mdf, df = df))
 }

##' This function is used to create HTML tables to present the results from a Gene Ontology (GO) analysis.
##'
##' After running a GO analysis, it is often useful to first present a table showing the set of significant GO terms
##' for a given comparison, and then have links to a sub-table for each GO term that shows the genes that were responsible
##' for the significance of that term. The first table can be generated using the \link[GOstats:GOHyperGResult-class]{summary} function, but
##' it will not contain the links to the sub-table. The ReportingTools package has functionality to make these tables and sub-tables
##' automatically, but the default is to include extra glyphs in the main table that are not that useful.
##' 
##'This function is intended to generate a more useful version of the table that one normally gets from ReportingTools.
##' @title Create HTML tables for Gene Ontology (GO) analyses
##' @param fit.table The output from \link[limma]{topTable}
##' @param go.summary The output from running \code{summary} on a \link[GOstats:GOHyperGResult-class]{GOHyperGResults} object.
##' @param probe.summary The output from running \link[GOstats]{probeSetSummary} on a \link[GOstats:GOHyperGResult-class]{GOHyperGResults} object.
##' @param cont.name The contrast name.
##' @param base.dir Character. Where should the HTML tables be generated? Defaults to GO_results.
##' @param extraname Character. An extra name that can be used if the contrast name isn't descriptive enough.
##' @param probecol  The column name in the topTable object that contains probe IDs. Defaults to PROBEID.
##' @param affy Boolean. Are the arrays from Affymetrix?
##' @return Returns an \link[ReportingTools:HTMLReportRef-class]{HTMLReportRef} object, which can be used when creating an index page to link
##' to the results.
##' @export
##' @author Jim MacDonald

makeGoTable <- function(fit.table, go.summary, probe.summary, cont.name, base.dir = "GO_results",
                        extraname = NULL, probecol = "PROBEID", affy = TRUE){
    htablst <- lapply(seq_len(nrow(go.summary)), function(x) makeGoGeneTable(fit.table, probe.summary[[x]],
                                                                             as.character(go.summary[x,1]),
                                                                             cont.name, base.dir, extraname,
                                                                             probecol, affy))
    go.summary[,1] <- hwrite(as.character(go.summary[,1]), 
                             link = sapply(htablst, function(x) gsub(paste0(base.dir, "/"), "", path(x))), table = FALSE)
    short.name <- paste0(gsub(" ", "_", cont.name), if(!is.null(extraname)) paste0("_", gsub(" ", "_", extraname)))
    long.name <- paste0(cont.name, if(!is.null(extraname)) paste(",", extraname))
    htab <- HTMLReport(short.name, long.name, reportDirectory = base.dir, baseDir = ".")
    publish(go.summary, htab, .modifyDF = goLinks)
    finish(htab)
    htab
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
