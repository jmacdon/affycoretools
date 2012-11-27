## code to create HTML tables for gene set analyses, specifically output from romer() in limma
## we start with lower level code and work our way up to top-level

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

geneSetPage <- function(rslts, genesets, eset, fit, file, cutoff = 0.05, dir = ".", subdir = ".",
                        columns = NULL, colnames = NULL, col = NULL, caption = NULL,
                        annot = NULL, nam = NULL, fitind = NULL, bline = NULL, ...){
    ## first go through results and select genesets with significant results
    ind <- apply(rslts[,c("Up","Down","Mixed")], 1, function(x) any(x < cutoff))
    rslts <- rslts[ind,]
    genesets <- genesets[ind]
    ind2 <- sapply(genesets, length) > 1
    rslts <- rslts[ind2,]
    names(rslts) <- paste(names(rslts), c("", rep(" (p-value)", 3), sep = "")) 
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
    cat("<script src=\"../sorttable.js\"></script>\n\n", paste(dir, file, sep = "/"))
    d.f <- data.frame(Genesets = names(genesets), rslts, Data = links, Heatmaps = hlinks)
    print(xtable(d.f, caption = caption, digits = rep(c(0,3,0), c(3,3,2))),
          caption.placement = "top", type = "html", file = paste(dir, file, sep = "/"),
          sanitize.text.function = function(x) x,
          include.rownames = FALSE, append = TRUE, html.table.attributes = c("border=1 class=\"sortable\""))
}

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
