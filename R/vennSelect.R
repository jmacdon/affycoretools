################################################
##
##  Copyright 2005 James W. MacDonald
##
##  Functions to both create Venn diagrams
##  and to output the subsets from the Venn
##  diagrams in HTML and text tables
##
################################################


vennCounts2 <- function(x, method = "same", fit = NULL,
                        foldFilt = NULL){

  ## x is a TestResults object from a call to
  ## decideTests()
  ## alternatively it can be the 'dirs' matrix from
  ## a call to foldFilt()
  ## fit is an MArrayLM object and foldFilt is a fold change to
  ## filter on. Note that fit is only required if foldFilt != NULL

  if(!is.null(foldFilt) && is.null(fit))
    stop("If you want to filter by fold change, you must also pass an MArrayLM object\n",
         "that was produced from a call to lmFit() and then eBayes().", call. = FALSE)
  if(!is.null(foldFilt)){
        idx <- abs(fit$coefficients) > foldFilt
        x <- as.matrix(x) * idx
    }
  tmp <- vennCounts(x)
  newcnts <- sapply(makeIndices(x, method), sum)
  all <- sum(tmp[,dim(tmp)[2]])
  if(dim(x)[2] == 2)
    ord <- c(2, 1, 3)
  if(dim(x)[2] == 3)
    ord <- c(3, 2, 6, 1, 5, 4, 7)
  tmp[,dim(tmp)[2]] <- c(all - sum(newcnts), newcnts[ord])
  tmp
}

intNames <- function(x){
  tmp <- colnames(x)
  if(is.null(dim(x)) || dim(x)[2] > 3)
    stop("This function only works for two or three comparisons at a time\n",
         call. = FALSE)
  if(dim(x)[2] == 2)
    return(paste(tmp[1], "and", tmp[2], sep = " "))
  if(dim(x)[2] == 3)
    return(paste(c(tmp[1], tmp[1], tmp[2]), "and", c(tmp[2], tmp[3], tmp[3])))
}

makeIndices <- function(x, method = "same"){
  method <- match.arg(method, c("same", "up", "down", "both", "sameup", "samedown"))
  indices <- list()
  x <- sign(switch(method, both = abs(x), up = x > 0, down = x < 0, same = x,
                   sameup = x, samedown = x))
  if(ncol(x) == 2){
    if(method != "sameup" && method != "samedown"){
      indices[[1]] <- abs(x[,1]) == 1 & x[,2] != x[,1]
      indices[[2]] <- x[,1] != x[,2] & abs(x[,2]) == 1
      indices[[3]] <- abs(rowSums(x)) == 2
    }
    if(method == "sameup"){
      indices[[1]] <- x[,1] == 1 & x[,2] != x[,1]
      indices[[2]] <- x[,1] != x[,2] & x[,2] == 1
      indices[[3]] <- rowSums(x) == 2
    }
    if(method == "samedown"){
      indices[[1]] <- x[,1] == -1 & x[,2] != x[,1]
      indices[[2]] <- x[,1] != x[,2] & x[,2] == -1
      indices[[3]] <- rowSums(x) == -2
    }
  }
  if(ncol(x) == 3){
    if(method != "sameup" && method != "samedown"){
      indices[[1]] <- abs(x[,1]) == 1 & x[,2] != x[,1] & x[,3] != x[,1]
      indices[[2]] <- x[,1] != x[,2] & abs(x[,2]) == 1 & x[,3] != x[,2]
      indices[[3]] <- x[,1] != x[,3] & x[,2] != x[,3] & abs(x[,3]) == 1
      indices[[4]] <- abs(rowSums(x[,1:2])) == 2 & x[,3] != x[,1]
      indices[[5]] <- abs(rowSums(x[,c(1,3)])) == 2 & x[,2] != x[,1]
      indices[[6]] <- abs(rowSums(x[,2:3])) == 2 & x[,1] != x[,2]
      indices[[7]] <- abs(rowSums(x)) == 3
    }
    if(method == "sameup"){
      indices[[1]] <- x[,1] == 1 & x[,2] != x[,1] & x[,3] != x[,1]
      indices[[2]] <- x[,1] != x[,2] & x[,2] == 1 & x[,3] != x[,2]
      indices[[3]] <- x[,1] != x[,3] & x[,2] != x[,3] & x[,3] == 1
      indices[[4]] <- rowSums(x[,1:2]) == 2 & x[,3] != x[,1]
      indices[[5]] <- rowSums(x[,c(1,3)]) == 2 & x[,2] != x[,1]
      indices[[6]] <- rowSums(x[,2:3]) == 2 & x[,1] != x[,2]
      indices[[7]] <- rowSums(x) == 3
    }
    if(method == "samedown"){
      indices[[1]] <- x[,1] == -1 & x[,2] != x[,1] & x[,3] != x[,1]
      indices[[2]] <- x[,1] != x[,2] & x[,2] == -1 & x[,3] != x[,2]
      indices[[3]] <- x[,1] != x[,3] & x[,2] != x[,3] & x[,3] == -1
      indices[[4]] <- rowSums(x[,1:2]) == -2 & x[,3] != x[,1]
      indices[[5]] <- rowSums(x[,c(1,3)]) == -2 & x[,2] != x[,1]
      indices[[6]] <- rowSums(x[,2:3]) == -2 & x[,1] != x[,2]
      indices[[7]] <- rowSums(x) == -3
    }
  }
  indices
}
      

getCols <- function(design, contrast){
  ## A function to get the correct columns of the ExpressionSet
  ## to output based on the design and contrasts matrices
  ## This is sort of a kludge and could probably be improved
  ncontrasts <- ncol(contrast)
  ids <- 1:dim(design)[1]
  tmp <- function(design, contrast, x){
    a <- design[,contrast[,x] > 0]
    if(is.matrix(a))
      a <- rowSums(abs(a)) > 0
    b <- design[,contrast[,x] < 0]
    if(is.matrix(b))
      b <- rowSums(abs(b)) > 0
    c(ids[as.logical(a)], ids[as.logical(b)])
  }
  cols <- list()
  if(ncontrasts == 2){
    for(i in 1:2)
      cols[[i]] <- tmp(design, contrast, i)
    cols[[3]] <- unique(unlist(cols[1:2]))
  }else{
    for(i in 1:3)
      cols[[i]] <- tmp(design, contrast, i)
    comb <- list(c(1, 2), c(1, 3), c(2, 3), 1:3)
    cols[4:7] <- lapply(comb, function(x) unique(unlist(cols[x])))
  }
  cols
}

getStat <- function(stat, ncontrasts, num, fit, contrast, index, adj.meth){

  ## function to get the correct statistics based on the cell of the Venn diagram

  if(ncontrasts == 2)
    tmp <- list(1, 2, 1:2)
  else
    tmp <- list(1, 2, 3, 1:2, c(1, 3), 2:3, 1:3)
  if(stat == "tstat")
    if(length(tmp[[num]]) == 1){
      out <- list("t-statistic" = fit$t[index,tmp[[num]]])
    }else{
      out <- fit$t[index,tmp[[num]]]
      if(is.matrix(out)){
          out <- as.list(as.data.frame(out))
      }else{
          out <- as.list(out)
      }
      names(out) <- paste("t-statistic for ", colnames(contrast)[tmp[[num]]], sep = "")
    }
  if(stat == "pval")
    if(length(tmp[[num]]) == 1){
      out <- list("p-value" = p.adjust(fit$p.value[,tmp[[num]]], adj.meth)[index])
    }else{
      out <- fit$p.value[,tmp[[num]]]
      if(is.matrix(out)){
          out <- as.list(as.data.frame(out))
      }else{
          out <- as.list(out)
      }
      out <- lapply(out, p.adjust, method = adj.meth)
      out <- lapply(out, function(x) x[index])
      names(out) <- paste("p-value for ", colnames(contrast)[tmp[[num]]], sep = "")
    }
  if(stat == "FC")
    if(length(tmp[[num]]) == 1){
      out <- list("log2 fold change" = fit$coefficients[index,tmp[[num]]])
    }else{
      out <- fit$coefficients[index,tmp[[num]], drop = FALSE]
      if(is.matrix(out)){
          out <- as.list(as.data.frame(out))
      }else{
          out <- as.list(out)
      }
      names(out) <- paste("log2 fold change for ", colnames(contrast)[tmp[[num]]], sep = "")
    }
  
  out
}


f.stat.dat <- function(fit, index, contrast, adj.meth = "BH", stats,
                       order.by = "pval", ncontrasts, num){

  ## A function to get the f-statistics, fold change, and p-value (or some subset)
  ## for each cell of the Venn diagram
 
  if(!order.by %in% stats)
    stop(paste("If you choose to order by ", order.by, " you must also include that as ",
               "one of the statistics to output.",sep = ""), call. = FALSE)
  
  if("fstat" %in% stats)
    f.stat <- list("F-statistic" = fit$F[index])
  if("pval" %in% stats){
    p.val <- list("p-value" = p.adjust(fit$F.p.value, method = adj.meth)[index])
    if(exists("f.stat")) out <- c(f.stat, p.val) else out <- p.val
  }
  if("FC" %in% stats){
    fc <- getStat("FC", ncontrasts, num, fit, contrast, index)
    if(exists("out")) out <- c(out, fc) else out <- fc
  }
  
  ord <- switch(order.by, pval = order(p.val[[1]]),
                fstat = order(f.stat[[1]], decreasing = TRUE),
                fc = order(fc[[1]], decreasing = TRUE))
  out <- lapply(out, function(x) x[ord])
  out <- lapply(out, function(x)
                ifelse(abs(x) < 1e-3, sprintf("%0.3e", x), sprintf("%0.3f", x)))
  return(list(out = out, ord = ord))
}

t.stat.dat <- function(fit, index, contrast, adj.meth = "BH", stats,
                       order.by = "pval", ncontrasts, num){
  ## A function to get the t-statistics, fold change, and p-value (or some subset)
  ## for each cell of the Venn diagram
 
  if(!order.by %in% stats)
    stop(paste("If you choose to order by ", order.by, " you must also include that as ",
               "one of the statistics to output.",sep = ""), call. = FALSE)
  
  if("tstat" %in% stats)
    t.stat <- getStat("tstat", ncontrasts, num, fit, contrast, index)
  if("pval" %in% stats){
    p.val <- getStat("pval", ncontrasts, num, fit, contrast, index, adj.meth)
    if(exists("t.stat")) out <- c(t.stat, p.val) else out <- p.val
  }
  if("FC" %in% stats){
    fc <- getStat("FC", ncontrasts, num, fit, contrast, index)
    if(exists("out")) out <- c(out, fc) else out <- fc
  }
  
  ord <- switch(order.by, pval = order(p.val[[1]]),
                tstat = order(t.stat[[1]], decreasing = TRUE),
                fc = order(fc[[1]], decreasing = TRUE))
  out <- lapply(out, function(x) x[ord])
  out <- lapply(out, function(x)
                ifelse(abs(x) < 1e-3, sprintf("%0.3e", x), sprintf("%0.3f", x)))
  return(list(out = out, ord = ord))
}


vennSelect <- function(eset, design, x, contrast, fit, method = "same", adj.meth = "BH",
                       stat = "fstat", otherstats = c("pval", "FC"), order.by = "pval",
                       foldFilt = NULL, save = FALSE, titleadd = NULL, ...){
  ## eset is an ExpressionSet containing data used for comparisons
  ## design is a design matrix from limma
  ## x is a TestResults object from a call to decideTests()
  ## output is a list containing the probe IDs of genes from each comparison

  if(!stat %in% c("fstat", "tstat")) stop(paste(stat, " is not an acceptable argument for the ",
                                               "'stat' argument.\nPlease use either ",
                                               "'fstat' or 'tstat'.", sep = ""), call. = FALSE)

  if(!is.null(foldFilt)){
    idx <- abs(fit$coefficients) > foldFilt
    x <- as.matrix(x) * idx
  }
  ncontrasts <- ncol(x)
  if(ncontrasts < 2 || ncontrasts > 3)
    stop("This function only works for two or three comparisons at a time.\n",
         call. = FALSE)
  if(ncontrasts == 2)
    name <- c(paste("Genes unique to", colnames(x)),
              paste("Genes in intersection of", intNames(x)))
  if(ncontrasts == 3)
    name <- c(paste("Genes unique to", colnames(x)),
              paste("Genes in intersection of", intNames(x)),
              "Genes common to all comparisons")
  
   ## Remove illegal characters from filenames
  if(length(grep("[/|\\|?|*|:|<|>|\"|\\|]", name)) > 0)
    warning(paste("Some illegal characters have been removed from the filenames",
                  name, sep = " "), call. = FALSE)
  name <- gsub("[/|\\|?|*|:|<|>|\"|\\|]", "", name)
  
  indices <- makeIndices(x, method = method)
  cols <- getCols(design, contrast)
  for(i in seq(along = indices)){
    tmp <- featureNames(eset)[indices[[i]]]
    if(stat == "fstat")
      otherdata <- f.stat.dat(fit, indices[[i]], contrast, adj.meth, c(stat, otherstats),
                              order.by, ncontrasts, i)
    if(stat == "tstat")
      otherdata <- t.stat.dat(fit, indices[[i]], contrast, adj.meth, c(stat, otherstats),
                              order.by, ncontrasts, i)
    if(is.null(stat)) otherdata <- NULL
    if(!is.null(titleadd))
        title <- paste(name[i], titleadd)
    else
        title <- name[i]
    if(length(tmp) == 0) next
    else
      probes2table(eset[,cols[[i]]], tmp[otherdata$ord], annotation(eset), text = TRUE,
                   filename = name[i], otherdata = otherdata$out, title = title, ...)
  }
  if(save)
    sapply(indices, sum)
}


vennSelect2 <- function(fit, contrast, design, eset, groups = NULL, cols = NULL, p.value = 0.05,
                        lfc = 0, method = "same", adj.meth = "BH", titleadd = NULL, fileadd = NULL, 
                        baseUrl = ".",  reportDirectory = "./venns", ...){

    ## design is a design matrix from limma
    ## contrast is a conrast matrix
    ## output is a list containing the probe IDs of genes from each comparison
    require("limma", character.only = TRUE, quietly = TRUE)
    require("lattice", character.only = TRUE, quietly = TRUE)
    require("ReportingTools", character.only = TRUE, quietly = TRUE)
    
    
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
    lapply(vennout, function(x) hwriter::hwrite("Legend for boxplot colors", p = page(x), 
                                                heading = 2))
    leglst <- lapply(seq(along = indices), function(x) 
                     makeLegend(groups[colind[[x]], drop = TRUE], 
                                reportDirectory(vennout[[x]]),
                                paste0("legend", x, ".png")))
    leglst <- lapply(leglst, function(x) hwriter::hwriteImage(x, links = x))
    lapply(seq(along = vennout), function(x) hwriter::hwrite(leglst[[x]], page(vennout[[x]]), br = TRUE))
    
    
    venncsv <- lapply(seq(along = indices), function(x) CSVFile(gsub(" ", "_", name[x]), 
                          reportDirectory = reportDirectory))
    require(annotation(eset), character.only = TRUE, quietly = TRUE)
    ps <- apply(fit$p.value, 2, p.adjust, method = adj.meth)
    colnames(ps) <- paste0(colnames(ps), ".p.value")
    coefs <- fit$coefficients
    colnames(coefs) <- paste0(colnames(coefs), ".logFC")
    csvlst <- lapply(seq(along = indices), function(x) cbind(coefs[indices[[x]], coefind[[x]], drop = FALSE],
                         ps[indices[[x]], coefind[[x]], drop = FALSE]))
    
    annotlst <- lapply(indices, function(x) select(get(annotation(eset)), as.character(featureNames(eset)[x]),
                                                   c("ENTREZID","SYMBOL")))
    csvlst <- mapply(data.frame, annotlst, csvlst, SIMPLIFY = FALSE)
    
    
    
    
    lapply(ind, function(x) publish(fit[indices[[x]],], vennout[[x]], eSet = eset[indices[[x]],colind[[x]]],
                                    factor = groups[colind[[x]], drop = TRUE],
                                    coef = coefind[[x]], adjust.method = adj.meth, pvalueCutoff = 1,
                                    n = Inf))
           
           
    lapply(ind, function(x) publish(csvlst[[x]], venncsv[[x]]))
           
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
                     baseUrl = ".", reportDirectory = NULL){
    require("ReportingTools", character.only = TRUE, quietly = TRUE)
    if(is.null(reportDirectory)){
        tmp <- strsplit(reportDirectory(vennlst[[1]]$vennout[[1]]), "/")[[1]]
        reportDirectory <- paste(tmp[-length(tmp)], collapse = "/")
    }
    hpage <- HTMLReport(pagename, pagetitle, baseUrl = baseUrl, reportDirectory = reportDirectory)
    hwriter::hwrite(paste("The Venn diagrams all contain clickable links. Click on the counts",
                          "in any cell to see a table of the genes in that cell.",
                          "Also please note that the tables are sortable - simply click",
                          "on any header to sort on that column."), page(hpage),
                    br = TRUE)
                        
    lapply(seq(along = vennlst), function(x) drawVenn(vennlst[[x]], page = page(hpage), 
               dir = reportDirectory, num = x, cex = cex.venn, shift.title = shift.title))
    finish(hpage)
    hpage
}

drawVenn <- function(lst, page, dir, num, cex = 1, shift.title = FALSE){
    require("ReportingTools", character.only = TRUE, quietly = TRUE)
    require("limma", character.only = TRUE, quietly = TRUE)
    require("hwriter", character.only = TRUE, quietly = TRUE)
    nam <- paste0(dir, "/venn", num, ".png")
    nam2 <- paste0("venn", num, ".png")
    mapname <- paste0("#venn", num)
    png(nam, height = 800, width = 800)
    if(shift.title) colnames(lst$venncounts)[1] <- paste0(colnames(lst$venncounts)[1], "\n\n")
    vennDiagram(lst$venncounts, cex = cex)
    dev.off()
    hwrite(paste("Venn Diagram", num), page, heading = 2)
    tag <- hmakeTag("img", border = 0, width = 800, height = 800,
                    src = nam2, alt = nam2, usemap = mapname)
    hwrite(tag, page)                
    #hwrite(tmp, page, br = TRUE)
    vennLinks(lst, page, mapname, paste0("./venn", num, "/"))
}

vennLinks <- function(lst, page, mapname, loc){
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
    urlst <- paste0(loc, sapply(lst$vennout, filename))
    strng <- do.call("c", mapply(fun, loclst, urlst, SIMPLIFY = FALSE))
    strng <- c(paste0('<map name="', sub("#", "", mapname), '">'),
               strng, "</map>")
    cat(strng, file = page, sep = "\n")
}
   
