## code to create HTML tables for gene set analyses, specifically output from romer() in limma
## we start with lower level code and work our way up to top-level



#' A function to create a simple heatmap and key.
#' 
#' This is an internal function called by \code{runRomer} and is not intended
#' to be used directly. It is documented here only because arguments may be
#' passed down via the \code{dots} argument.
#' 
#' As noted above, this is only intended to be called indirectly by
#' \code{runRomer}. However, certain arguments such as scale.row, or col, etc,
#' can be passed down to this function via the \code{dots} argument, allowing
#' the end user to have more control over the finished product.
#' 
#' @param eset An \code{\link{ExpressionSet}} containing normalized, summarized
#' gene expression data.
#' @param ind Numeric vector indicating which rows of the
#' \code{\link{ExpressionSet}} to use.
#' @param filename The filename for the heatmap and associated key.
#' @param columns Numeric vector indicating which columns of the
#' \code{\link{ExpressionSet}} to use. If \code{NULL}, all columns will be
#' used.
#' @param colnames Character. Substitute column names for the heatmap. If
#' \code{NULL}, the sampleNames will be used.
#' @param col A vector of colors to use for the heatmap. If \code{NULL}, the
#' \code{\link{bluered}} function will be used.
#' @param annot Character. The name of the array annotation package to use to
#' convert probeset IDs to gene symbols. If \code{NULL}, either the probeset
#' IDs will be used, or the input from the 'nam' argument below will be used.
#' @param nam Character. A vector of row.names to use for the heatmap. If
#' \code{NULL}, and the above 'annot' argument is \code{NULL} as well, probeset
#' IDs will be used.
#' @param scale.row Boolean. Should the data be scaled by row? Defaults to
#' \code{FALSE}.
#' @param key Boolean. Should a key be produced that shows the numeric range
#' for the colors of the heatmap? Defaults to \code{TRUE}.
#' @param bline A numeric vector, usually extracted from a contrast matrix,
#' used to \sQuote{sweep} the mean baseline sample means from the heatmap data.
#' The end result will be a heatmap in which the colors correspond to log fold
#' changes from the baseline samples.
#' @return Nothing is returned. Called only for the side effect of creating
#' heatmaps in 'png' format.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
gsHeatmap <- function(eset, ind, filename, columns = NULL, colnames = NULL, col = NULL,
                      annot = NULL, nam = NULL, scale.row = FALSE, key = TRUE, bline = NULL){
    if(!is.null(annot) && !is.null(nam))
        stop(paste("You can specify an annotation package (annot argument) to use for row names\n",
                   "or directly specify the row names (nam argument), but not both!\n\n", sep = ""),
             call. = FALSE)
    if(!is.null(bline) && scale.row)
        stop(paste("The scale.row and bline arguments are mutually exclusive - you can either\n",
                   "scale each row of the heatmap, or you can convert to fold changes, but not both.\n\n"),
             call. = FALSE)
    if(!is.null(annot))
        require(annot, quietly = TRUE, character.only = TRUE)
    
    if(is.null(columns)) columns <- 1:dim(eset)[2]
    if(is(eset, "ExpressionSet"))
        mat <- exprs(eset)[ind,columns, drop = FALSE]
    else
        mat <- eset[ind,columns, drop = FALSE]
    if(nrow(mat) < 2) return(NULL)
    if(!is.null(colnames)) colnames(mat) <- colnames
    if(!is.null(annot)){
        rn <- sapply(AnnotationDbi::mget(row.names(mat), get(paste(sub("\\.db", "", annot),
                                                                   "SYMBOL", sep = ""))), "[", 1)
        rn <- ifelse(is.na(rn), row.names(mat), rn)
        row.names(mat) <- rn
    }
    ## fragile - there is no assurance that the order of the matrix matches the names supplied!
    if(!is.null(nam)) row.names(mat) <- nam
    if(scale.row){
        mat <- sweep(mat, 1, rowMeans(mat))
        mat <- sweep(mat, 1, apply(mat, 1, sd), "/")
    }
     if(!is.null(bline)){
        mn <- rowMeans(mat[,columns %in% bline])
        mat <- sweep(mat, 1, mn)
    }
   
    if(is.null(col))
        col <- bluered(25)
    if(key){
        png(sub("\\.png", "_key.png", filename))
        z <- seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length = length(col))
        image(z = matrix(z, ncol = 1), col = col, xaxt = "n", yaxt = "n")
        axis(1, at = c(0,0.5, 1), labels = round(quantile(z, c(0, 0.5, 1), include.lowest = TRUE),1),
             cex.axis = 3)
        dev.off()
    }
    png(filename, pointsize = 18, width = 480 + ncol(mat)*10, height = 480 + nrow(mat)*10)
    par(mar = c(8.1,3.1,3.1,6.1))
    image(1:ncol(mat), 1:nrow(mat), t(mat[nrow(mat):1,]), col = col, axes = FALSE,
          xlab = "", ylab = "", main = "Gene set\nHeatmap",
          ylim = c(0, nrow(mat)) + 0.5, xlim = c(0, ncol(mat)) + 0.5)
    axis(1, at = 1:ncol(mat), colnames(mat), las = 2, tick = 0)
    axis(4, at = nrow(mat):1, row.names(mat), tick = 0,
         cex.axis = 0.2 + 1/log10(nrow(mat)), line = -0.5, las = 2)
    dev.off()
}



#' A function to create an annotated HTML table for all genes in a significant
#' gene set as well as a heatmap of these data.
#' 
#' This is intended to be an internal function to \code{runRomer}. It is
#' documented here only because it may be necessary to pass alternative
#' arguments to this function from \code{runRomer}.
#' 
#' This function creates an annotation table using \code{probes2table} if an
#' annotation file is used, otherwise data will be output in a simple HTML
#' table. A heatmap showing the expression values for all the genes in the gene
#' set is then placed below this table, along with a key that indicates the
#' range of the expression values.
#' 
#' @param eset An \code{\link{ExpressionSet}} containing normalized, summarized
#' gene expression data.
#' @param fit An \code{MArrayLM} object.
#' @param ind Numeric vector indicating which rows of the
#' \code{\link{ExpressionSet}} to use.
#' @param columns Numeric vector indicating which columns of the
#' \code{\link{ExpressionSet}} to use. If \code{NULL}, all columns will be
#' used.
#' @param annot Character. The name of the annotation package for the array
#' used. If \code{NULL}, data will be output without any additional annotation.
#' @param fname The filename of the resulting output, without the 'html' file
#' extension.
#' @param heatmap Character. The filename of the heatmap to append to the
#' bottom of the HTML page.
#' @param title Title to be placed at the top of the resulting HTML page.
#' @param key Character. The filename of the heatmap key to append to the
#' bottom of the HTML page.
#' @param fitind Numeric. Which column of the \code{MArrayLM} object to use for
#' output in the HTML table.
#' @author James W. MacDonald
dataAndHeatmapPage <- function(eset, fit, ind, columns = NULL, annot = NULL, fname, heatmap, title,
                               key = TRUE, fitind = NULL){
    fnhtml <- paste(fname, "_heatmap", sep = "")
    prbs <- featureNames(eset)[ind]
    if(is.null(columns)) columns <- 1:dim(eset)[2]
    if(is.null(fitind)) fitind <- 1:ncol(fit$t)
    otherdata <- list(t.statistic = fit$t[ind,fitind],
                      p.value = fit$p.value[ind,fitind],
                      fold = fit$coefficients[ind,fitind])
    ind2 <- order(otherdata$p.value)
    otherdata <- lapply(otherdata, function(x) x[ind2])
    prbs <- prbs[ind2]
    if(!is.null(annot)){
        probes2table(eset = eset[,columns], probids = prbs, lib = annot,
                     otherdata = otherdata, filename = fname, title = title,
                     anncols = annaffy::aaf.handler()[c(1:3,6,7,9)])
    }else{
        d.f <- data.frame(Probesets = prbs,
                          do.call("cbind", otherdata),
                          if(is(eset, "ExpressionSet")) exprs(eset)[ind,columns] else eset[ind,columns])
        print(xtable(d.f, caption = title), file = fnhtml, type = "html", include.rownames = FALSE)
    }
    target <- HTMLInitFile(".", fnhtml, Title = paste(title, "heatmap"))
    cat(paste("\n <h>", paste(title, "heatmap"), "</h>\n<br><br><br>", sep = ""), file = target, append = TRUE)
    if(key) HTMLInsertGraph(sub("\\.png", "_key.png", heatmap), caption = "Heatmap key",
                            Align = "left", WidthHTML = 300, HeightHTML = 150, file = target)
    HTMLInsertGraph(heatmap, file = target, WidthHTML = NULL)
    return(ind[ind2])
}



#' A function to create an HTML page for each gene set, as well as the HTML
#' pages for each significant gene set.
#' 
#' This is intended to be an internal function to \code{runRomer}, and is not
#' intended to be called by end users. However, the \dots{} argument to
#' \code{runRomer} allows one to pass arguments to lower level functions, so
#' the arguments are described here.
#' 
#' This function creates a \sQuote{midlevel} HTML table that contains each gene
#' set that was significant, with a link to an HTML table that shows data for
#' each gene in that gene set (with annotation), as well as a heatmap showing
#' the expression levels. Normally this is not run by end users, but is called
#' as part of the \code{runRomer} function.
#' 
#' @param rslts The results from running \code{\link{romer}} on one gene set.
#' @param genesets Character. A vector of gene symbols for one gene set.
#' @param eset An \code{\link{ExpressionSet}} containing normalized, summarized
#' gene expression data.
#' @param fit An \code{MArrayLM} object, containing the fitted data.
#' @param file Filename for the resulting HTML page.
#' @param cutoff Numeric. The cutoff for significance for a given gene set.
#' Defaults to 0.05.
#' @param dir The directory to write the results. Defaults to the working
#' directory.
#' @param subdir The subdirectory to write the individual gene set results.
#' Defaults to the working directory.
#' @param columns Numeric. The columns of the \code{ExpressionSet} to use for
#' the individual gene set output pages. See \code{dataAndHeatmapPage} for more
#' information.
#' @param colnames Character. Alternative column names for the resulting
#' heatmap. See \code{dataAndHeatmapPage} for more information.
#' @param col A vector of colors for the heatmap. Defaults to
#' \code{\link{bluered}}.
#' @param caption Caption to put at the top of the HTML page.
#' @param annot Character. The name of the annotation package for the arrray
#' used. If \code{NULL}, the row names of the heatmap will be the probeset IDs
#' unless the \sQuote{nam} argument below is used.
#' @param nam Character. The row names for the resulting heatmap. If
#' \code{NULL} and no \sQuote{annot} argument is used, the probeset IDs will be
#' used for the row names.
#' @param fitind Numeric. The columns of the \code{MArrayLM} object to use for
#' the individual HTML tables.
#' @param bline Defaults to \code{NULL}. Otherwise, a numeric vector indicating
#' which columns of the data are the baseline samples. The data used for the
#' heatmap will be centered by subtracting the mean of these columns from all
#' data.
#' @param \dots Allows arguments to be passed to lower-level functions. See
#' \code{dataAndHeatmapPage} and \code{gsHeatmap} for available arguments.
#' @return Nothing is returned. Called only for the side effect of creating
#' HTML tables.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
geneSetPage <- function(rslts, genesets, eset, fit, file, cutoff = 0.05, dir = ".", subdir = ".",
                        columns = NULL, colnames = NULL, col = NULL, caption = NULL,
                        annot = NULL, nam = NULL, fitind = NULL, bline = NULL, ...){
    ## first go through results and select genesets with significant results
    ind <- apply(rslts[,c("Up","Down","Mixed")], 1, function(x) any(x < cutoff))
    rslts <- rslts[ind,]
    genesets <- genesets[ind]
    ind2 <- sapply(genesets, length) > 1
    rslts <- rslts[ind2,]
    colnames(rslts) <- paste(colnames(rslts), c("", rep(" (p-value)", 3)))
    genesets <- genesets[ind2]
    fun <- function(x, ...){
        fn <- if(subdir == ".") paste(dir, x, sep = "/") else paste(dir, subdir, x, sep = "/")
        ind2 <- dataAndHeatmapPage(eset = eset, fit = fit, ind = genesets[[x]], columns = columns,
                                   annot = annot, fname = fn,
                                   heatmap = paste(x, "png", sep = "."),
                                   title = paste(x, "geneset,", if(subdir != ".") paste(subdir, "contrast")),
                                   fitind = fitind)
        gsHeatmap(eset = eset, ind = ind2, filename = paste(fn, ".png", sep = ""),
                  columns = columns, colnames = colnames, col = col,
                  annot = annot, nam = nam, bline = bline, ...)
        
        return(if(subdir == ".") paste(x, "html", sep = ".") else paste(subdir, "/", x, ".html", sep = ""))
    }
    fnames <- lapply(names(genesets), fun)
    links <- paste("<a href=\"", fnames, "\">Data table</a>", sep = "")
    hlinks <- paste("<a href=\"", gsub("\\.html", "_heatmap.html", fnames), "\">Heatmap</a>", sep = "")
    cat("<script src=\"../sorttable.js\"></script>\n\n", file = paste(dir, file, sep = "/"))
    d.f <- cbind(names(genesets), rslts, links, hlinks)
    colnames(d.f) <- c("Genesets", colnames(rslts), "Data","Heatmaps")
    print(xtable(d.f, caption = caption, digits = rep(c(0,3,0), c(3,3,2))),
          caption.placement = "top", type = "html", file = paste(dir, file, sep = "/"),
          sanitize.text.function = function(x) x,
          include.rownames = FALSE, append = TRUE, html.table.attributes = c("border=1 class=\"sortable\""))
}



#' A function to create HTML output from the results of running romer on a set
#' of contrasts.
#' 
#' This function is actually intended to be a sub-function of \code{runRomer},
#' but can hypothetically run by itself if the \code{\link{romer}} step has
#' already been done.
#' 
#' This function is intended to be an internal function for \code{runRomer}.
#' However, it is possible that \code{runRomer} errored out after saving the
#' results from running \code{\link{romer}} on a set of contrasts, and all that
#' remains is to create the output HTML.
#' 
#' Please note that the first two arguments to this function have certain
#' expectations. The rsltlst should be the output from running
#' \code{\link{romer}}. If using the saved output from \code{runRomer}, one
#' should first \code{load} the 'romer.Rdata' file, which will introduce a list
#' object with the name 'romerlst' into the working directory, so the first
#' argument should be rsltlst = romerlst.
#' 
#' Second, see the code for runRomer, specifically the line that creates the
#' 'sets' object, which will show how to create the correct genesetlst object.
#' 
#' @param rsltlst A list of results, generated by the \code{\link{romer}}
#' function. See discussion for more information.
#' @param genesetlst A list of genesets, usually created by loading in the
#' RData files that can be downloaded from
#' http://bioinf.wehi.edu.au/software/MSigDB/. See details for more
#' information.
#' @param eset An \code{\link{ExpressionSet}} containing normalized, summarized
#' gene expression data.
#' @param fit An \code{MArrayLM} object, containing the fitted data.
#' @param design A design matrix describing the model.
#' @param contrast A contrast matrix describing the contrasts that were fit.
#' This matrix should have colnames, which will be used to name subdirectories
#' containing results.
#' @param changenames Boolean. When creating heatmaps of the gene sets, should
#' the columns be appended with the colnames from the design matrix? If
#' \code{FALSE}, the sampleNames will be used.
#' @param dir Character. The subdirectory to use for the output data. Defaults
#' to 'genesets'.
#' @param explanation If \code{NULL}, a generic paragraph will be placed at the
#' top of the indexRomer.html page, giving a brief explanation of the analysis.
#' Alternatively, this can be replaced with other text. Please note that this
#' text should conform to HTML standards (e.g., will be pasted into the HTML
#' document as-is, so should contain any required HTML markup).
#' @param annot Character. The name of the annotation package for the array
#' being used. This is used to add annotation data to the data tables for each
#' gene set. If \code{NULL} (the default), the tables will only contain
#' probeset IDs and statistics.
#' @param baseline.hmap Boolean. If \code{TRUE}, then the resulting heatmaps
#' will be centered by subtracting the mean of the baseline sample. As an
#' example, in a contrast of treatment A - treatment B, the mean of the
#' treatment B samples will be subtracted. The heatmap colors then represent
#' the fold change between the A and B samples.
#' @param file Character. The filename to output. Defaults to
#' 'indexRomer.html'.
#' @param \dots Arguments to be passed to lower-level functions. See
#' \code{geneSetPage}, \code{dataAndHeatmapPage} and \code{gsHeatmap} for
#' available arguments.
#' @return Nothing is returned. The function is run only for the side effect of
#' creating HTML tables with output for each significant gene set.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @export outputRomer
outputRomer <- function(rsltlst, genesetlst, eset, fit, design = NULL, contrast = NULL, changenames = TRUE,
                        dir = "genesets", explanation = NULL, annot = NULL, baseline.hmap = TRUE,
                        file = "indexRomer.html", ...){
    if(!is.null(design) && !is.null(contrast)){
        ## here I am assuming that if design and contrast are used, it's because this is
        ## being run by runRomer, so I don't have to check that things line up. Otherwise caveat emptor
        cols <- lapply(1:ncol(contrast), function(x) c(which(as.logical(design[,contrast[,x] > 0])),
                                                       which(as.logical(design[,contrast[,x] < 0]))))
        if(changenames)
            nams <- lapply(1:ncol(contrast), function(x) c(rep(colnames(design)[contrast[,x] > 0],
                                                               sum(design[,contrast[,x] > 0])),
                                                           rep(colnames(design)[contrast[,x] < 0],
                                                               sum(design[,contrast[,x] < 0]))))
        if(baseline.hmap)
            bline <- lapply(1:ncol(contrast), function(x) which(as.logical(design[,contrast[,x] < 0])))
        else
            bline <- NULL
    }
    
    crushit <- function(x) sapply(strsplit(x, " "), paste, collapse = "_")
    ## fixit <- function(x) sapply(strsplit(x, "\\/|_|\\."), grep, pattern = "c1|c2|c3|c4|c5", value = TRUE)
    topnam <- crushit(names(rsltlst))
    nextnam <- names(rsltlst[[1]])
    ## mapper <- c("C1-positional","C2-Curated","C3-Motif","C4-Computational","C5-GeneOntology","C6-Oncogenic")
    ## names(mapper) <- paste("c", 1:6, sep = "")
    ## nextnam2 <- mapper[nextnam]
    for(i in seq(along = topnam)){
        for(j in seq(along = nextnam)){
            if(!file.exists(dir))
                dir.create(dir)
            if(!file.exists(paste(dir, topnam[[i]], sep = "/")))
                dir.create(paste(dir, topnam[[i]], sep = "/"))
            geneSetPage(rslts = rsltlst[[i]][[j]], genesets = genesetlst[[j]], eset = eset,
                        fit = fit, file = paste(paste(topnam[i], nextnam[j], sep = "_"), "html", sep = "."),
                        dir = dir, subdir = topnam[[i]], fitind = i,
                        columns = cols[[i]], colnames = nams[[i]],
                        caption = paste("Results for", topnam[i], "contrast, using", nextnam[j], "geneset"),
                        annot = annot, bline = bline[[i]], ...)
        }
    }

    contlinks <- paste(dir, "/", paste(rep(topnam, each = length(nextnam)),
                                       rep(nextnam, length(topnam)), sep = "_"),
                       ".html", sep = "")
    contrastlnks <- paste("<a href=\"", contlinks, "\">", rep(topnam, each = length(nextnam)),
                          "</a>", sep = "")
    if(is.null(explanation)){
        target <- HTMLInitFile(".", sub("\\.html", "", file),"html", BackGroundColor = "#FFFFFF",
                               Title = "Gene Set Enrichment Analysis based on Broad Gene Sets")
        cat('\n<script type="text/javascript" src="sorttable.js"></script>',file=target, append = TRUE)
        cat('\n<style type="text/css">',file=target, append = TRUE)
        cat("\n<!--",file=target, append = TRUE)
        cat("\nbody,td,th {",file=target, append = TRUE)
        cat("\n	font-size: small;",file=target, append = TRUE)
        cat("\n	font-family: Verdana, Geneva, sans-serif;",file=target, append = TRUE)
        cat("\n	color: #000",file=target, append = TRUE)
            cat("\n  width:100px;", file = target, append = TRUE)
        cat("\n}",file=target, append = TRUE)
        cat("\n-->",file=target, append = TRUE)
        cat("\n</style>",file=target, append = TRUE)
        HTML(paste("<br>This analysis is based on several gene set collections from the",
                   "<a href=\"http://www.broadinstitute.org/gsea/msigdb/index.jsp\">Broad Institute MSigDB:</a>.",
                   "The basic idea is to look for sets of genes that appear to be either higher, lower, or",
                   "mixed (some higher, some lower) in a ranked list of genes than would be expected by chance.",
                   "<br> In other words, if we were to randomly select a group of genes, we would expect them",
                   "to be uniformly distributed in a ranked set of genes (where the genes are ranked based on",
                   "a particular comparison of samples). If the set of genes in aggregate is found to be closer",
                   "to the top of the gene list, then we assume this didn't occur by chance, and assume that",
                   "that set of genes is in aggregate up-regulated, and possibly perturbed in the comparison of",
                   "interest.<br>The Broad Institute genesets are arranged in six collections, based on either",
                   "genetic location (C1-Positional), being part of a pathway (C2-Curated), being in a motif",
                   "based on conserved cis-regulatory motifs (C3-Motif), computational gene sets based on",
                   "expression of cancer associated genes (C4-Computational), Gene Ontology terms",
                   "(C5-GeneOntology), or oncogenic signatures from cancer microarray experiments (C6-Oncogenic).",
                   "<br>The analysis is based on the romer function in the Bioconductor limma package.",
                   "Note that the table below is sortable - simply click on any header to re-sort."),
             file= target)
        
    }else{
        target <- HTMLInitFile(".", sub("\\.html", "", file),"html", BackGroundColor = "#FFFFFF",
                               Title = "Gene Set Enrichment Analysis based on Broad Gene Sets")
        HTML(explanation, file = target)
    }
    d.f <- data.frame(Genesets = rep(nextnam, length(topnam)),
                      Contrasts = contrastlnks)
    print(xtable(d.f), include.rownames = FALSE, type = "html", file = file,
          sanitize.text.function = function(x) x, append = TRUE,
          html.table.attributes = "border=1, class=sortable")
}

## runRomer <- function(setloc, annot = NULL, eset, design = NULL, contrast = NULL, fit, wts = NULL, save = TRUE,
##                      baseline.hmap = TRUE, ...){
##     require(annot, quietly = TRUE, character.only = TRUE)
##     sets <- sapply(paste(setloc, dir(setloc, "[Rr]data"), sep = "/"), load, envir = .GlobalEnv)
##     symb <- select(get(annot), featureNames(eset), "SYMBOL")
##     symb <- symb[!is.na(symb$SYMBOL),]
##     symb <- symb[!duplicated(symb$SYMBOL),]
##     eset <- eset[as.character(symb$PROBEID),]
##     fit <- fit[as.character(symb$PROBEID),]
##     setlst <- lapply(sets, function(x) symbols2indices(get(x), as.character(symb$SYMBOL)))
##     romerlst <- lapply(seq_len(ncol(contrast)), function(x)
##                        lapply(setlst, romer, y = exprs(eset), design = design, contrast = contrast[,x],
##                               if(!is.null(wts)) array.weights = wts))
##     names(romerlst) <- colnames(contrast)
##     if(save) save(list = "romerlst", file = "romer.Rdata")
##     outputRomer(rsltlst = romerlst, genesetlst = setlst, eset = eset, fit = fit,
##                 design = design, contrast = contrast, annot = annot, baseline.hmap = baseline.hmap, ...)
## }



#' A function to run the romer function on a set of contrasts.
#' 
#' This function automates both running \code{\link{romer}} on a set of
#' contrasts as well as the creation of output HTML tables that can be used to
#' explore the results.  The basic idea here is that one might have used limma
#' to fit a model and compute some contrasts, and then want to do a GSEA using
#' \code{\link{romer}}.
#' 
#' The \code{\link{romer}} expects as input a list or lists of gene symbols
#' that represent individual gene sets. One example is the various gene sets
#' from the Broad Institute that are available at
#' http://bioinf.wehi.edu.au/software/MSigDB/, which are distributed as RData
#' files. The default assumption for this function is that the end user will
#' have downloaded these files, and the setloc argument simply tells
#' \code{runRomer} where to find them.
#' 
#' Alternatively, user-based gene sets could be created (these should consist
#' of lists of character vectors of gene symbols - see one of the Broad gene
#' sets for an example).
#' 
#' This function will run \code{\link{romer}} using all the gene sets in the
#' referenced directory, on all the contrasts supplied, and then output the
#' results in a (default) 'genesets' subdirectory. There will be an HTML file
#' in the working directory with a (default) filename of 'indexRomer.html' that
#' will point to individual HTML files in the genesets subdirectory, which will
#' point to individual files in subdirectories within the genesets subdirectory
#' (named after the colnames of the contrast matrix).
#' 
#' @param setloc A character vector giving the path for gene set RData files
#' (see description for more information).
#' @param annot Character. The name of the array annotation package.
#' @param eset An \code{\link{ExpressionSet}} containing normalized expression
#' data.
#' @param design A design matrix describing the model fit to the data.
#' @param contrast A contrast matrix describing the contrasts that were
#' computed from the data. This contrast should have colnames, which will be
#' used to create parts of the resulting directory structure.
#' @param fit An \code{MArrayLM} object, containing the fitted model data.
#' @param wts Optional weights vector - if array weights were used to fit the
#' model, they should be supplied here as well.
#' @param save Boolean. If true, after running the \code{\link{romer}} step,
#' the results will be saved in a file 'romer.Rdata', which can be used as
#' input for \code{outputRomer} to create HTML tables. Since
#' \code{\link{romer}} can take a long time to run, it is advantageous to keep
#' the default.
#' @param baseline.hmap Boolean. If \code{TRUE}, then the resulting heatmaps
#' will be centered by subtracting the mean of the baseline sample. As an
#' example, in a contrast of treatment A - treatment B, the mean of the
#' treatment B samples will be subtracted. The heatmap colors then represent
#' the fold change between the A and B samples.
#' @param \dots Used to pass arguments to lower-level functions. See
#' \code{outputRomer} \code{geneSetPage}, \code{dataAndHeatmapPage} and
#' \code{gsHeatmap} for available arguments.
#' @return Nothing is returned. This function is called only for the
#' side-effects of creating output HTML files in the working and
#' sub-directories.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @export runRomer
runRomer <- function(setloc, annot = NULL, eset, design = NULL, contrast = NULL, fit, wts = NULL, save = TRUE,
                     baseline.hmap = TRUE, ...){
    require(annot, quietly = TRUE, character.only = TRUE)
    thefile <- paste(setloc, dir(setloc, "Rdata$"), sep = "/")
    if(!file.exists(thefile)) stop(paste("Ensure that the Broad dataset exists in", thefile), call. = FALSE)
    broad <- get(load(thefile))
    symb <- select(get(annot), featureNames(eset), "SYMBOL")
    symb <- symb[!is.na(symb$SYMBOL),]
    symb <- symb[!duplicated(symb$SYMBOL),]
    eset <- eset[as.character(symb$PROBEID),]
    fit <- fit[as.character(symb$PROBEID),]
    setlst <- lapply(broad, function(x) symbols2indices(x, as.character(symb$SYMBOL)))
    names(setlst) <- names(broad)
    romerlst <- lapply(seq_len(ncol(contrast)), function(x)
                       lapply(setlst, romer, y = exprs(eset), design = design, contrast = contrast[,x],
                              if(!is.null(wts)) array.weights = wts))
    names(romerlst) <- colnames(contrast)
    if(save) save(list = "romerlst", file = "romer.Rdata")
    outputRomer(rsltlst = romerlst, genesetlst = setlst, eset = eset, fit = fit,
                design = design, contrast = contrast, annot = annot, baseline.hmap = baseline.hmap, ...)
}
