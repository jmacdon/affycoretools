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
    if(ncontrasts < 2 || ncontrasts > 3)
        stop("This function only works for two or three comparisons at a time.\n",
             call. = FALSE)
    if(ncontrasts == 2)
        name <- c(paste("Genes unique to", colnames(dtmat)),
                  paste("Genes in intersection of", intNames(dtmat)))
    if(ncontrasts == 3)
        name <- c(paste("Genes unique to", colnames(dtmat)),
                  paste("Genes in intersection of", intNames(dtmat)),
                  "Genes common to all comparisons")
    
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
    ## add legend to HTML pages here
   ## browser()
    ## lapply(vennout, function(x) publish(hwrite("Legend for boxplot colors", 
    ##                                             heading = 2), x))
    ## leglst <- lapply(seq(along = indices), function(x) 
    ##                  makeLegend(groups[colind[[x]], drop = TRUE], 
    ##                             reportDirectory,
    ##                             paste0("legend", x, ".png")))
    ## leglst <- lapply(leglst, function(x) hwriteImage(x, links = x))
    ## lapply(seq(along = vennout), function(x) publish(hwrite(leglst[[x]], br = TRUE), vennout[[x]]))
    
    
    ## venncsv <- lapply(seq(along = indices), function(x) CSVFile(gsub(" ", "_", name[x]), 
    ##                       reportDirectory = reportDirectory))
    require(annotation(eset), character.only = TRUE, quietly = TRUE)
    ps <- apply(fit$p.value, 2, p.adjust, method = adj.meth)
    colnames(ps) <- paste0(colnames(ps), ".p.value")
    coefs <- fit$coefficients
    colnames(coefs) <- paste0(colnames(coefs), ".logFC")
    csvlst <- lapply(seq(along = indices), function(x) cbind(coefs[indices[[x]], coefind[[x]], drop = FALSE],
                         ps[indices[[x]], coefind[[x]], drop = FALSE]))
    
    annotlst <- lapply(indices, function(x) select(get(annotation(eset)), as.character(featureNames(eset)[x]),
                                                   c("ENTREZID","SYMBOL","GENENAME")))
    csvlst <- mapply(data.frame, annotlst, csvlst, SIMPLIFY = FALSE)
    
    
    ## lapply(ind, function(x) publish(fit[indices[[x]],], vennout[[x]], eSet = eset[indices[[x]],colind[[x]]],
    ##                                 factor = groups[colind[[x]], drop = TRUE],
    ##                                 coef = coefind[[x]], adjust.method = adj.meth, pvalueCutoff = 1,
    ##                                 n = Inf))
    ## lapply(ind, function(x) publish(csvlst[[x]], venncsv[[x]]))

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
    vennlst <- lapply(seq(along = collist), function(x) 
                      vennSelect2(fit = fit, contrast = contrast, design = design, eset = eset,
                                 cols = collist[[x]], p.value = p.value, lfc = lfc,
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

drawVenn <- function(lst, page, dir, num, cex = 1, shift.title = FALSE){
    nam <- paste0(dir, "/venn", num, ".png")
    nam2 <- paste0("venn", num, ".png")
    mapname <- paste0("#venn", num)
    png(nam, height = 800, width = 800)
    if(shift.title) colnames(lst$venncounts)[1] <- paste0(colnames(lst$venncounts)[1], "\n\n")
    vennDiagram(lst$venncounts, cex = cex)
    dev.off()
    
    hwrite(hmakeTag("img", border = 0, width = 800, height = 800,
                    src = nam2, alt = nam2, usemap = mapname), page)
    hwriteImage(paste("Venn Diagram", num), page)
    vennLinks(lst, page, mapname, dir)
    
}

vennLinks <- function(lst, page, mapname, reportDirectory){
    if(!ncol(lst$venncounts) %in% 3:4)
        stop(paste("You can currently only create Venn diagrams for two and three",
                   "comparisons.\n"), call. = FALSE)
    fun <- function(x,y) paste0('<area shape="circle" coords=',
                                    x, ' href=', y, '></area>')
    if(ncol(lst$venncounts) == 3){
        loclst <- list("250,400,30","550,400,30","400,400,30")
    } else {
        loclst <- list("260,310,30","540,300,30","400,550,30",
                       "400,320,30","330,420,30","470,420,30",
                       "400,390,30")
    }
    urlst <- gsub(reportDirectory, ".", sapply(lst$vennout, path))
    strng <- do.call("c", mapply(fun, loclst, urlst, SIMPLIFY = FALSE))
    strng <- c(paste0('<map name="', sub("#", "", mapname), '">'),
               strng, "</map>")
    cat(strng, file = page, sep = "\n")
}
   
