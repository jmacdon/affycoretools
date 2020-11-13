##' This function is designed to parse a design and contrasts matrix in order to
##' generate Glimma's interactive MA-plots for all contrasts. The heading for the
##' resulting plot is based on the colnames of the contrasts matrix, so it is important
##' to include the colnames and make them explanatory.
##'
##' In addition, if there are multiple contrasts with the same name, say if the same
##' comparisons are being made for different tissue types, the extraname argument will cause
##' the output to be placed in glimma-plots/<extraname>, to eliminate over-writing of existing
##' files.
##' @title A function to generate MA-plots from Glimma, for all contrasts.
##' @param tablst An MArrayLM object or list of DGEExact or DGELRT objects
##' @param datobj A DGEList, ExpressionSet, EList or matrix
##' @param dsgn A design matrix
##' @param cont A contrast matrix
##' @param grpvec A vector of groups, usually what was used to make the design matrix 
##' @param padj Method for multiplicity adjustment. BH by default
##' @param sigfilt Significance cutoff for selecting genes
##' @param extraname Used to add a sub-directory to the glimma-plots directory, mainly used to
##' disambiguate contrasts with the same name (see below).
##' @param ID Character, default is NULL. The column name of the genes list item (from the MArrayLM,
##' DGEExact or DGELRT objects) that will be used for labeling the individual gene expression plot
##' at the top right of the resulting plot. The default of NULL will search for a column labeled
##' SYMBOL to use. If the ID value doesn't match any column, and there is no SYMBOL column, then
##' the default for glMDPlot will be used.
##' @param sample.cols A vector of numbers or hex strings, the same length as the grpvec. This
##' will cause the samples in the upper right plot to have different colors. Useful primarily to
##' check for batch differences. The default is the default from glMDPlot.
##' @param ... Allows end users to pass other arguments to the Glimma glMDPlot function
##' @return A character vector of the files generated, useful for using as links to the output.
##' @author James W. MacDonald \email{jmacdon@@u.washington.edu}
##' @keywords manip
##' @examples
##'   \dontrun{
##'     mat <- matrix(rnorm(1e6), ncol = 20)
##'     grps <- factor(1:4, each=5)
##'     design <- model.matrix(~0 + grps)
##'     colnames(design) <- LETTERS[1:4]
##'     contrast <- matrix(c(1,-1,0,0,1,0,-1,0,1,0,0,-1,0,1,-1,0,0,1,0,-1),
##'     ncol = 5)
##'     colnames(contrast) <- paste(LETTERS[c(1,1,1,2,2)],
##'     LETTERS[c(2,3,4,3,4)], sep = " vs ")
##'     fit <- lmFit(mat, design)
##'     fit2 <- contrasts.fit(fit, contrast)
##'     fit2 <- eBayes(fit2)
##'     htmllinks <- doGlimma(fit2, mat, design, contrast, grps)
##'    }
##' @export doGlimma
doGlimma <- function(tablst, datobj, dsgn, cont, grpvec, padj = "BH", sigfilt = 0.05,
                     extraname = NULL, ID = NULL, sample.cols = rep("#1f77b4", ncol(datobj)), ...){
    getSymb <- function(x, ID){
        if(!is.null(ID)) {
            symb <- grep(ID, colnames(x$genes), ignore.case = TRUE, value = TRUE)
        } else {
            symb <- grep("symbol", colnames(x$genes), ignore.case = TRUE, value = TRUE)
        }
        if(length(symb) == 1) return(symb) else return(NULL)
    }

    counts <- switch(class(datobj)[1],
                     DGEList = datobj$counts,
                     ExpressionSet = exprs(datobj),
                     EList = datobj$E,
                     matrix = datobj,
                     stop(paste("Please use either a DGEList, ExpressionSet, matrix",
                                "or EList object for the datobj argument."), call. = FALSE))
    for(i in seq_len(ncol(cont))){
        folder <- if(is.null(extraname)) "glimma-plots" else paste0("glimma-plots/", extraname)
        if(!file.exists(folder)) dir.create(folder, recursive = TRUE)
        html <- gsub(" ", "_", colnames(cont))
        ind <- as.logical(dsgn %*% cont[,i])
        if(is(tablst[[i]], "DGELRT") | is(tablst[[i]],"DGEExact")){
            symb <- getSymb(tablst[[i]], ID)
            status <- decideTests(tablst[[i]], p.value = sigfilt, adjust.method = padj)
            glMDPlot(tablst[[i]], counts = counts[,ind], groups = factor(grpvec[ind]), status = status,
                     transform = TRUE, folder = folder, side.main = symb,
                     html = html[i], launch = FALSE, main = colnames(cont)[i], p.adj.method = padj,
                     sample.cols = sample.cols[ind,], ...)
        } else if(is(tablst, "MArrayLM")) {
            symb <- getSymb(tablst, ID)
            status <- decideTests(tablst, p.value = sigfilt, adjust.method = padj, coefficients = i)
            glMDPlot(tablst, counts = counts[,ind], groups = factor(grpvec[ind]), status = status, coef = i,
                     transform = FALSE, folder = folder, side.main = symb,
                     html = html[i], launch = FALSE, main = colnames(cont)[i], p.adj.method = padj,
                     sample.cols = sample.cols[ind], ...)
        } else {
            stop("Please provide either a DGELRT, DGExact or MArrayLM object!", call. = FALSE)
        }
    }
    return(paste0(paste(folder, html, sep = "/"), ".html"))
}
