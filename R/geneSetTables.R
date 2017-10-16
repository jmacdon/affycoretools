## code to create HTML tables for gene set analyses, specifically output from romer() in limma
## we start with lower level code and work our way up to top-level



##' A function to create a simple heatmap and key.
##' 
##' This is an internal function called by \code{runRomer} and is not intended
##' to be used directly. It is documented here only because arguments may be
##' passed down via the \code{dots} argument.
##' 
##' As noted above, this is only intended to be called indirectly by
##' \code{runRomer}. However, certain arguments such as scale.row, or col, etc,
##' can be passed down to this function via the \code{dots} argument, allowing
##' the end user to have more control over the finished product.
##' 
##' @param eset An \code{\link{ExpressionSet}} containing normalized, summarized
##' gene expression data.
##' @param ind Numeric vector indicating which rows of the
##' \code{\link{ExpressionSet}} to use.
##' @param filename The filename for the heatmap and associated key.
##' @param columns Numeric vector indicating which columns of the
##' \code{\link{ExpressionSet}} to use. If \code{NULL}, all columns will be
##' used.
##' @param colnames Character. Substitute column names for the heatmap. If
##' \code{NULL}, the sampleNames will be used.
##' @param col A vector of colors to use for the heatmap. If \code{NULL}, the
##' \code{\link{bluered}} function will be used.
##' @param annot A matrix or data.frame containing gene symbols to annotate the heatmap.
##' This will normally be extracted automatically from the 'fit' object passed to geneSetPage.
##' If there is no annotation in the fit object, then the probe IDs will be used instead.
##' @param scale.row Boolean. Should the data be scaled by row? Defaults to
##' \code{FALSE}.
##' @param key Boolean. Should a key be produced that shows the numeric range
##' for the colors of the heatmap? Defaults to \code{TRUE}.
##' @param bline A numeric vector, usually extracted from a contrast matrix,
##' used to \sQuote{sweep} the mean baseline sample means from the heatmap data.
##' The end result will be a heatmap in which the colors correspond to log fold
##' changes from the baseline samples.
##' @return Nothing is returned. Called only for the side effect of creating
##' heatmaps in 'png' format.
##' @author James W. MacDonald <jmacdon@@u.washington.edu>
gsHeatmap <- function(eset, ind, filename, columns = NULL, colnames = NULL, col = NULL,
                      annot = NULL, scale.row = FALSE, key = TRUE, bline = NULL){
    if(!is.null(bline) && scale.row)
        stop(paste("The scale.row and bline arguments are mutually exclusive - you can either\n",
                   "scale each row of the heatmap, or you can convert to fold changes, but not both.\n\n"),
             call. = FALSE)
    if(is.null(columns)) columns <- 1:dim(eset)[2]
    if(is(eset, "ExpressionSet"))
        mat <- exprs(eset)[ind,columns, drop = FALSE]
    else
        mat <- eset[ind,columns, drop = FALSE]
    if(nrow(mat) < 2) return(NULL)
    if(!is.null(colnames)) colnames(mat) <- colnames
    if(!is.null(annot)){
        if("SYMBOL" %in% toupper(colnames(annot))){
            rn <- annot[,grep("symbol", colnames(annot), ignore.case = TRUE)]
            rn <- ifelse(is.na(rn), row.names(mat), rn)
            row.names(mat) <- rn
        }
    }
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



##' A function to create an annotated HTML table for all genes in a significant
##' gene set as well as a heatmap of these data.
##' 
##' This is intended to be an internal function to \code{runRomer}. It is
##' documented here only because it may be necessary to pass alternative
##' arguments to this function from \code{runRomer}.
##' 
##' This function creates an annotation table using \code{probes2table} if an
##' annotation file is used, otherwise data will be output in a simple HTML
##' table. A heatmap showing the expression values for all the genes in the gene
##' set is then placed below this table, along with a key that indicates the
##' range of the expression values.
##' 
##' @param object An ExpressionSet, containing normalized, summarized
##' gene expression data.
##' @param fit An MArrayLM object containing the fitted data.
##' @param ind Numeric vector indicating which rows of the
##' data object to use.
##' @param columns Numeric vector indicating which columns of the
##' data object to use. If \code{NULL}, all columns will be
##' used.
##' @param fname The filename of the resulting output, without the 'html' file
##' extension.
##' @param heatmap Character. The filename of the heatmap to append to the
##' bottom of the HTML page.
##' @param title Title to be placed at the top of the resulting HTML page.
##' @param key Character. The filename of the heatmap key to append to the
##' bottom of the HTML page.
##' @param fitind Numeric. Which column of the \code{MArrayLM} object to use for
##' output in the HTML table.
##' @param affy Boolean. Are these Affymetrix arrays? If \code{TRUE}, then links
##' will be generated to netaffx for the probeset IDs.
##' @param \dots Included to allow arbitrary commands to be passed to lower level functions.
##' @author James W. MacDonald <jmacdon@u.washington.edu>

dataAndHeatmapPage <- function(object, fit, ind, columns = NULL, fname, heatmap, title,
                               key = TRUE, fitind = NULL, affy = TRUE, ...){
    fnhtml <- paste(fname, "_heatmap.html", sep = "")
    prbs <- featureNames(object)[ind]
    if(is.null(columns)) columns <- 1:dim(object)[2]
    if(is.null(fitind)) fitind <- 1:ncol(fit$coefficients)
    
    otherdata <- switch(grep("t|LR|F", names(fit), value = TRUE),
                        t = data.frame(t.statistic = fit$t[ind,fitind],
                                       p.value = fit$p.value[ind,fitind],
                                       fold = fit$coefficients[ind,fitind]),
                        LR = data.frame(LR.statistic = fit$LR[ind,fitind],
                                       p.value = fit$p.value[ind,fitind],
                                       fold = fit$coefficients[ind,fitind]),
                        F =  data.frame(F.statistic = fit$F[ind,fitind],
                                       p.value = fit$p.value[ind,fitind],
                                       fold = fit$coefficients[ind,fitind]))
    ind2 <- order(otherdata$p.value)
    otherdata <- otherdata[ind2,]
    prbs <- prbs[ind2]
    ## Here we just rely on any annotations that exist in the MArrayLM object
    rn <- row.names(fit$coef)
    if(is.null(rn)) rn <- seq_len(nrow(fit$coef))
    if(is.null(fit$genes)){
        warning(paste("\nThere are no annotations in the MArrayLM object, using",
                      if(is.null(row.names(fit$coef))) "row indices" else "row names",
                      "from the MArrayLM object instead"))
        annot <- data.frame(ID = rn[ind])
    } else {
        annot <- fit$genes[ind,]
    }
    d.f <- cbind(annot, otherdata)
    d.f <- fixHeaderAndGo(d.f, affy = affy, ...)
    tab <- HTMLReport(fname, title)
    publish(d.f$df, tab, if(length(d.f$mdf) > 0) .modifyDF = d.f$mdf)
    finish(tab)
    page <- openPage(fnhtml, title = paste(title, "heatmap"))
    if (key) hwriteImage(sub("\\.png", "_key.png", heatmap), center = FALSE,
                         width = 300, height = 150, page = page)
    hwriteImage(heatmap, page)
    closePage(page)
    return(list(ind = ind[ind2], annot = annot))
}


    

##' A function to create an HTML page for each gene set, as well as the HTML
##' pages for each significant gene set.
##' 
##' This is intended to be an internal function to \code{runRomer}, and is not
##' intended to be called by end users. However, the \dots{} argument to
##' \code{runRomer} allows one to pass arguments to lower level functions, so
##' the arguments are described here.
##' 
##' This function creates a \sQuote{midlevel} HTML table that contains each gene
##' set that was significant, with a link to an HTML table that shows data for
##' each gene in that gene set (with annotation), as well as a heatmap showing
##' the expression levels. Normally this is not run by end users, but is called
##' as part of the \code{runRomer} function.
##' 
##' @param rslts The results from running \code{\link{romer}} on one gene set.
##' @param genesets Character. A vector of gene symbols for one gene set.
##' @param object An ExpressionSet, DGEList or EList containing normalized, summarized
##' gene expression data.
##' @param fit An MArrayLM or DGEGLM object, containing the fitted data.
##' @param file Filename for the resulting HTML page.
##' @param cutoff Numeric. The cutoff for significance for a given gene set.
##' Defaults to 0.05.
##' @param dir The directory to write the results. Defaults to the working
##' directory.
##' @param subdir The subdirectory to write the individual gene set results.
##' Defaults to the working directory.
##' @param columns Numeric. The columns of the \code{ExpressionSet} to use for
##' the individual gene set output pages. See \code{dataAndHeatmapPage} for more
##' information.
##' @param colnames Character. Alternative column names for the resulting
##' heatmap. See \code{dataAndHeatmapPage} for more information.
##' @param col A vector of colors for the heatmap. Defaults to
##' \code{\link{bluered}}.
##' @param caption Caption to put at the top of the HTML page.
##' @param fitind Numeric. The columns of the \code{MArrayLM} object to use for
##' the individual HTML tables.
##' @param bline Defaults to \code{NULL}. Otherwise, a numeric vector indicating
##' which columns of the data are the baseline samples. The data used for the
##' heatmap will be centered by subtracting the mean of these columns from all
##' data.
##' @param affy Boolean; are these Affymetrix arrays? If \code{TRUE}, the Affymetrix
##' probeset IDs will contain links to the netaffx site.
##' @param \dots Allows arguments to be passed to lower-level functions. See
##' \code{dataAndHeatmapPage} and \code{gsHeatmap} for available arguments.
##' @return Nothing is returned. Called only for the side effect of creating
##' HTML tables.
##' @author James W. MacDonald <jmacdon@@u.washington.edu>
geneSetPage <- function(rslts, genesets, object, fit, file, cutoff = 0.05, dir = ".", subdir = ".",
                        columns = NULL, colnames = NULL, col = NULL, caption = NULL,
                        fitind = NULL, bline = NULL, affy = TRUE,  ...){
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
        out <- dataAndHeatmapPage(eset = eset, fit = fit, ind = genesets[[x]], columns = columns,
                                   fname = fn, heatmap = paste(x, "png", sep = "."),
                                   title = paste(x, "geneset,", if(subdir != ".") paste(subdir, "contrast")),
                                   fitind = fitind, affy = affy, ...)
        gsHeatmap(object = object, ind = out$ind, filename = paste(fn, ".png", sep = ""),
                  columns = columns, colnames = colnames, col = col,
                  annot = out$annot, bline = bline, ...)
        
        return(if(subdir == ".") paste(x, "html", sep = ".") else paste(subdir, "/", x, ".html", sep = ""))
    }
    fnames <- lapply(names(genesets), fun)
    links <- hwrite(rep("Data table", length(fnames)), link = fnames, table = FALSE)
    hlinks <- hwrite(rep("Heatmap", length(fnames)), link = gsub("\\.html","_heatmap.html", fnames), table = FALSE)
    d.f <- data.frame(Genesets = names(genesets), rslts, Data = links, Heatmaps = hlinks)
    tab <- HTMLReport(sub("\\.html$", "", file), caption, dir)
    publish(d.f, tab)
    finish(tab)
    tab
}




##' @param rsltlst A list of results, generated by the \code{\link{romer}}
##' function. See discussion for more information.
##' @param genesetlst A list of genesets, usually created by loading in the
##' RData files that can be downloaded from
##' http://bioinf.wehi.edu.au/software/MSigDB/. See details for more
##' information.
##' @param design A design matrix describing the model.
##' @param contrast A contrast matrix describing the contrasts that were fit.
##' This matrix should have colnames, which will be used to name subdirectories
##' containing results.
##' @param changenames Boolean. When creating heatmaps of the gene sets, should
##' the columns be appended with the colnames from the design matrix? If
##' \code{FALSE}, the sampleNames will be used.
##' @param dir Character. The subdirectory to use for the output data. Defaults
##' to 'genesets'.
##' @param explanation If \code{NULL}, a generic paragraph will be placed at the
##' top of the indexRomer.html page, giving a brief explanation of the analysis.
##' Alternatively, this can be replaced with other text. Please note that this
##' text should conform to HTML standards (e.g., will be pasted into the HTML
##' document as-is, so should contain any required HTML markup).
##' @param baseline.hmap Boolean. If \code{TRUE}, then the resulting heatmaps
##' will be centered by subtracting the mean of the baseline sample. As an
##' example, in a contrast of treatment A - treatment B, the mean of the
##' treatment B samples will be subtracted. The heatmap colors then represent
##' the fold change between the A and B samples.
##' @param file Character. The filename to output. Defaults to
##' indexRomer.html.
##' @param affy Boolean. Are these Affymetrix arrays? if \code{TRUE}, then
##' thre will be links generated in the HTML table to the netaffx site.
##' @param \dots Arguments to be passed to lower-level functions. See
##' \code{geneSetPage}, \code{dataAndHeatmapPage} and \code{gsHeatmap} for
##' available arguments.
##' @return Nothing is returned. The function is run only for the side effect of
##' creating HTML tables with output for each significant gene set.
##' @describeIn outputRomer Output romer results using microarray data
setMethod("outputRomer", c("ExpressionSet","MArrayLM"),
          function(object, fit, rsltlst, genesetlst, design = NULL, contrast = NULL, changenames = TRUE,
                   dir = "genesets", explanation = NULL, baseline.hmap = TRUE, 
                   file = "indexRomer.html", affy = TRUE,  ...){
    ## Ensure we are using a cell means design
    if(any(apply(design, 1, sum) > 1)) stop(paste("The design matrix needs to be for a cell-means model (e.g., one '1' per row)",
                                                  "for this function to work correctly. If you have an intercept or batch variables",
                                                  "please construct a simplified design matrix that just specifies the individual sample",
                                                  "types.\n"), call. = FALSE)
    
    if(!is.null(design) && !is.null(contrast)){
        ## here I am assuming that if design and contrast are used, it's because this is
        ## being run by runRomer, so I don't have to check that things line up. Otherwise caveat emptor
        cols <- lapply(1:ncol(contrast), function(x) which(apply(design[,contrast[,x] != 0, drop = FALSE], 1, function(y) any(y != 0))))
        
        if(changenames)
            nams <- lapply(cols, function(x) colnames(design)[apply(design[x,],1,function(y) which(y != 0))])
        
        if(baseline.hmap)
            bline <- lapply(1:ncol(contrast), function(x) which(apply(design[,contrast[,x] < 0, drop = FALSE], 1, function(y) any(y != 0))))
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
                        bline = bline[[i]], affy = affy,...)
        }
    }

    contlinks <- paste(dir, "/", paste(rep(topnam, each = length(nextnam)),
                                       rep(nextnam, length(topnam)), sep = "_"),
                       ".html", sep = "")
    contrastlnks <- paste("<a href=\"", contlinks, "\">", rep(topnam, each = length(nextnam)),
                          "</a>", sep = "")
    indx <- HTMLReport(sub("\\.html", "", file), "Gene Set Enrichment Analysis based on Broad Gene Sets", reportDirectory = dir)
    d.f <- data.frame(Genesets = rep(nextnam, length(topnam)),
                      Contrasts = contrastlnks)
    if(is.null(explanation))
        explanation <- paste("This analysis is based on several gene set collections from the",
                             hwrite("Broad Institute MSigDB;", link = "http://www.broadinstitute.org/gsea/msigdb/index.jsp"),
                             "The basic idea is to look for sets of genes that appear to be either higher, lower, or",
                             "mixed (some higher, some lower) in a ranked list of genes than would be expected by chance.",
                             "In other words, if we were to randomly select a group of genes, we would expect them",
                             "to be uniformly distributed in a ranked set of genes (where the genes are ranked based on",
                             "a particular comparison of samples). If the set of genes in aggregate is found to be closer",
                             "to the top of the gene list, then we assume this didn't occur by chance, and assume that",
                             "that set of genes is in aggregate up-regulated, and possibly perturbed in the comparison of",
                             "interest.<br>The Broad Institute genesets are arranged in six collections, based on either",
                             "genetic location (C1-Positional), being part of a pathway (C2-Curated), being in a motif",
                             "based on conserved cis-regulatory motifs (C3-Motif), computational gene sets based on",
                             "expression of cancer associated genes (C4-Computational), Gene Ontology terms",
                             "(C5-GeneOntology), or oncogenic signatures from cancer microarray experiments (C6-Oncogenic).",
                             "The analysis is based on the romer function in the Bioconductor limma package.",
                             "Note that the table below is sortable - simply click on any header to re-sort.")

        
    publish(hwrite(explanation, br = TRUE), indx)
    publish(d.f, indx)
    finish(indx)
})

##' @describeIn outputRomer Output romer results using RNA-Seq data processed using edgeR
setMethod("outputRomer", c("DGEList","DGEGLM"),
          function(object, fit, rsltlst, genesetlst, design = NULL, contrast = NULL, changenames = TRUE,
                   dir = "genesets", explanation = NULL, baseline.hmap = TRUE, 
                   file = "indexRomer.html", ...){
    ## Convert to an ExpressionSet and MArrayLM and run through base method
    object <- ExpressionSet(cpm(object, log = TRUE), featureData = AnnotatedDataFrame(object$genes))
    if(is.null(fit$var.prior)){
        fitlst <- lapply(1:ncol(contrast)) function(x) topTags(glmLRT(fit, contrast = contrast[,x]), Inf, sort.by = "none")$table
        fitnew <- list(coefficient = do.call(cbind(lapply(fitlst, function(x) x$logFC))),
                       p.value = do.call(cbind(lapply(fitlst, function(x) x$PValue))),
                       LR = do.call(cbind(lapply(fitlst, function(x) x$LR))))
    } else {
        fitlst <- lapply(1:ncol(contrast)) function(x) topTags(glmQLFTest(fit, contrast = contrast[,x]), Inf, sort.by = "none")$table
        fitnew <- list(coefficient = do.call(cbind(lapply(fitlst, function(x) x$logFC))),
                       p.value = do.call(cbind(lapply(fitlst, function(x) x$PValue))),
                       F = do.call(cbind(lapply(fitlst, function(x) x$F))))
    }
    outputRomer(object, new("MArrayLM", fitnew), ...)
})

##` @describeIn outputRomer Output romer results using RNA-Seq data processed using voom.
setMethod("outputRomer", c("EList","MArrayLM"),
          function(object, fit, rsltlst, genesetlst, design = NULL, contrast = NULL, changenames = TRUE,
                   dir = "genesets", explanation = NULL, baseline.hmap = TRUE, 
                   file = "indexRomer.html", ...){
    ## Convert to an ExpressionSet and run through base method
    object <- ExpressionSet(object$E, featureData = AnnotatedDataFrame(object$genes))
    outputRomer(object, fit,...)
})


##' @param fit A fitted model from either limma (e.g., MArrayLM) or edgeR (e.g., DGEGLM)
##' @param setloc A character vector giving the path for gene set RData files
##' (see description for more information), or a named list (or list of lists),
##' where the top-level names consist of gene set grouping names (like KeGG or GO),
##' the next level names consist of gene set names (like NAKAMURA_CANCER_MICROENVIRONMENT_UP),
##' and the list items themselves are gene symbols, matching the expected capitalization for
##' the species being used (e.g., for human, they are ALL CAPS. For most other species only the
##' First Letter Is Capitalized).
##' @param annot Character. The name of the array annotation package.
##' @param eset An \code{\link{ExpressionSet}} containing normalized expression
##' data.
##' @param design A design matrix describing the model fit to the data.
##' @param contrast A contrast matrix describing the contrasts that were
##' computed from the data. This contrast should have colnames, which will be
##' used to create parts of the resulting directory structure.
##' @param fit An \code{MArrayLM} object, containing the fitted model data.
##' @param wts Optional weights vector - if array weights were used to fit the
##' model, they should be supplied here as well.
##' @param save Boolean. If true, after running the \code{\link{romer}} step,
##' the results will be saved in a file 'romer.Rdata', which can be used as
##' input for \code{outputRomer} to create HTML tables. Since
##' \code{\link{romer}} can take a long time to run, it is advantageous to keep
##' the default.
##' @param baseline.hmap Boolean. If \code{TRUE}, then the resulting heatmaps
##' will be centered by subtracting the mean of the baseline sample. As an
##' example, in a contrast of treatment A - treatment B, the mean of the
##' treatment B samples will be subtracted. The heatmap colors then represent
##' the fold change between the A and B samples.
##' @param affy Boolean; are these Affymetrix arrays? If \code{TRUE}, the output tables
##' will contain links to the netaffx site.
##' @param \dots Used to pass arguments to lower-level functions. See
##' \code{outputRomer} \code{geneSetPage}, \code{dataAndHeatmapPage} and
##' \code{gsHeatmap} for available arguments.
##' @return Nothing is returned. This function is called only for the
##' side-effects of creating output HTML files in the working and
##' sub-directories.
##' @describeIn runRomer Perform gene set analysis using microarray data.
##' @author James W. MacDonald <jmacdon@@u.washington.edu>
##' @export 
setMethod("runRomer", c("ExpressionSet"),
          function(object, fit, setloc, annot = NULL,  design = NULL, contrast = NULL,  wts = NULL, save = TRUE,
                     baseline.hmap = TRUE, affy = TRUE, ...){
    ## Ensure we are using a cell means design
    if(any(apply(design, 1, sum) > 1)) stop(paste("The design matrix needs to be for a cell-means model (e.g., one '1' per row)",
                                                  "for this function to work correctly. If you have an intercept or batch variables",
                                                  "please construct a simplified design matrix that just specifies the individual sample",
                                                  "types.\n"), call. = FALSE)
   
    broad <- switch(class(setloc),
                    character = get(load(paste(setloc, dir(setloc, "Rdata$"), sep = "/"))),
                    list = setloc,
                    stop("The setloc argument should be either a file location or a named list.", call. = FALSE))
    if(!is.null(annot)){
        require(annot, character.only = TRUE, quietly = TRUE)
        object <- annotateEset(object, annot)
    }
    symbind <- grep("symbol", names(fData(object)), ignore.case = TRUE)
    if(length(symbind) == 1){
        symb <- fData(object)[,symbind]
        dontuse <- is.na(symb) | duplicated(symb)
        symb <- symb[-dontuse]
        object <- object[-dontuse,]
        fit <- fit[-dontuse,]
        } else {
            stop("Please either specify an annotation package, or supply gene symbols in the fData slot of your ExpressionSet.", call. = FALSE)
        }
    }
        
    setlst <- lapply(broad, function(x) ids2indices(x, as.character(symb)))
    names(setlst) <- names(broad)
    romerlst <- lapply(seq_len(ncol(contrast)), function(x)
                       lapply(setlst, romer, y = exprs(object), design = design, contrast = contrast[,x],
                              if(!is.null(wts)) array.weights = wts))
    names(romerlst) <- colnames(contrast)
    if(save) save(list = "romerlst", file = "romer.Rdata")
    outputRomer(object = object, fit = fit, rsltlst = romerlst, genesetlst = setlst, 
                design = design, contrast = contrast, baseline.hmap = baseline.hmap, affy = affy, ...)
})



##' @describeIn runRomer Perform gene set analysis using RNA-Seq data processed using edgeR.
##' @export
setMethod("runRomer", c("DGEList"),
          function(object, fit, setloc, design = NULL, contrast = NULL,
                   save = TRUE, baseline.hmap = TRUE, ...){
    if (any(apply(design, 1, sum) > 1))
        stop(paste("The design matrix needs to be for a cell-means model (e.g., one '1' per row)", 
            "for this function to work correctly. If you have an intercept or batch variables", 
            "please construct a simplified design matrix that just specifies the individual sample", 
            "types.\n"), call. = FALSE)
    broad <- switch(class(setloc),
                    character = get(load(paste(setloc, dir(setloc, "Rdata$"), sep = "/"))),
                    list = setloc,
                    stop("The setloc argument should be either a file location or a named list.", call. = FALSE))    
    symbind <- grep("SYMBOL", colnames(object$genes), ignore.case = TRUE)
    if(length(symbind != 1))
        stop("Please include gene symbols in the annotation for your DGEList!", call. = FALSE)
    
    symb <- object$genes[,symbind]
    ind <- 1:nrow(object)
    ind  <- ind[is.na(symb)|duplicated(symb)]
    object <- object[-ind,]
    fit <- fit[-ind, ]
    setlst <- lapply(broad, function(x) ids2indices(x, as.character(object$genes[,symbind])))
    names(setlst) <- names(broad)
    romerlst <- lapply(seq_len(ncol(contrast)), function(x) lapply(setlst, 
                       romer, y = object, design = design, contrast = contrast[, x], ...))
    names(romerlst) <- colnames(contrast)
    if (save) 
        save(list = "romerlst", file = "romer.Rdata")
    outputRomer(object = object, fit = fit, rsltlst = romerlst, genesetlst = setlst,
                design = design, contrast = contrast, 
                baseline.hmap = baseline.hmap, ...)
})


##' @describeIn runRomer Perform gene set analysis using RNA-Seq data processed using voom.
##' @export
setMethod("runRomer", c("EList"),
          function(object, fit, setloc, design = NULL, contrast = NULL,
                   save = TRUE, baseline.hmap = TRUE, ...){
     callNextMethod()
})


    
