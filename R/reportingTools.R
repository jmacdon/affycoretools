## code for integrating with ReportingTools package


entrezLinks <- function(df, ...){
    df$ENTREZID <- hwriter::hwrite(as.character(df$ENTREZID),
                                   link = paste0("http://www.ncbi.nlm.nih.gov/gene/",
                                   as.character(df$ENTREZID)), table = FALSE)
    return(df)
}

affyLinks <- function(df, ...){
    df$PROBEID <- hwriter::hwrite(as.character(df$PROBEID),
                                  link = paste0("https://www.affymetrix.com/LinkServlet?probeset=", 
                                  as.character(df$PROBEID)), table = FALSE)
    return(df)
}

goLinks <- function(df, ...){
    df$Term <- hwriter::hwrite(as.character(df$Term),
                               link = paste0("http://www.godatabase.org/cgi-bin/amigo/go.cgi?view=details&search_constraint=terms&depth=0&query=",
                               as.character(df$Term)), table = FALSE)
    return(df)
}



vennSelect2 <- function(fit, contrast, design, eset, groups = NULL, cols = NULL, p.value = 0.05,
                        lfc = 0, method = "same", adj.meth = "BH", titleadd = NULL, fileadd = NULL, 
                        baseUrl = ".",  reportDirectory = "./venns", ...){

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
    require(annotation(eset), character.only = TRUE, quietly = TRUE)
    ps <- apply(fit$p.value, 2, p.adjust, method = adj.meth)
    colnames(ps) <- paste0(colnames(ps), ".p.value")
    coefs <- fit$coefficients
    colnames(coefs) <- paste0(colnames(coefs), ".logFC")
    csvlst <- lapply(seq(along = indices), function(x) cbind(coefs[indices[[x]], coefind[[x]], drop = FALSE],
                         ps[indices[[x]], coefind[[x]], drop = FALSE]))
    
    annotlst <- lapply(indices, function(x) if(sum(x) > 0) AnnotationDbi::select(get(annotation(eset)), as.character(featureNames(eset)[x]),
                                                                  c("ENTREZID","SYMBOL","GENENAME")) else NULL)
    csvlst <- mapply(data.frame, annotlst, csvlst, SIMPLIFY = FALSE)
    
    
    lapply(ind, function(x) publish(csvlst[[x]], vennout[[x]], .modifyDF = list(entrezLinks, affyLinks)))
    vennpaths <- gsub("html$", "txt", sapply(vennout, path))
    lapply(ind, function(x) write.table(csvlst[[x]], vennpaths[x],
                                        sep = "\t", quote = FALSE, row.names = FALSE, na = ""))
    
    lapply(vennout[ind], finish)
    return(list(vennout = vennout, venncounts = venncounts))
}

           
    
makeLegend <- function(groups, dir, fname){
    nam <- paste(dir, fname, sep = "/")
    png(nam, height = 250, width = 300)
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
    cols <- lattice.options()$default.theme$superpose.symbol$col
    legend(1,1, levels(groups), pch = 16, col = cols[sort(unique(as.numeric(groups)))],
           xjust = 0.5, yjust=0.5, cex = 1.5, bty="n")
    dev.off()
    fname
}

makeVenn <- function(fit, contrast, design, eset, groups = NULL, collist = NULL,
                     p.value = 0.05, lfc = 0, method = "both", adj.meth = "BH",
                     titleadd = NULL, fileadd = NULL, baseUrl = ".", reportDirectory = "./venns",
                     ...){
    if(is.null(collist))
        vennlst <- vennSelect2(fit = fit, contrast = contrast, design = design, eset = eset,
                               groups = groups, cols = seq_len(ncol(fit$coefficients)), p.value = p.value, lfc = lfc,
                               method = method, adj.meth = adj.meth, titleadd = titleadd,
                               fileadd = fileadd, baseUrl = baseUrl, 
                               reportDirectory = paste0(reportDirectory, "/venn", x), ...)
    else
        vennlst <- lapply(seq(along = collist), function(x) 
                          vennSelect2(fit = fit, contrast = contrast, design = design, eset = eset,
                                      groups = groups, cols = collist[[x]], p.value = p.value, lfc = lfc,
                                      method = method, adj.meth = adj.meth, titleadd = titleadd,
                                      fileadd = fileadd, baseUrl = baseUrl, 
                                      reportDirectory = paste0(reportDirectory, "/venn", x), ...))
    vennlst
}

vennPage <- function(vennlst, pagename, pagetitle, cex.venn = 1, shift.title = FALSE,
                     baseUrl = ".", reportDirectory = NULL, ...){
    if(is.null(reportDirectory)){
        tmp <- strsplit(path(vennlst[[1]]$vennout[[1]]), "/")[[1]]
        reportDirectory <- paste(tmp[1:2], collapse = "/")
    }
    hpage <- hwriter::openPage(paste0(pagename, ".html"), dirname = reportDirectory)
    hwrite(paste("The Venn diagrams all contain clickable links. Click on the counts",
                          "in any cell to see a table of the genes in that cell.",
                          "Also please note that the tables are sortable - simply click",
                          "on any header to sort on that column."), hpage,
                    br = TRUE)
                        
    lapply(seq(along = vennlst), function(x) drawVenn(vennlst[[x]], page = hpage, 
               dir = reportDirectory, num = x, cex = cex.venn, shift.title = shift.title))
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

