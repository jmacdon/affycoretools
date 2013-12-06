###################################################
##
##  Copyright 2005 James W. MacDonald
##
##  MA plotting functions
##
####################################################




#' A Function to make MA plots from all arrays.
#' 
#' This function creates an MA plot for all arrays in either an ExpressionSet
#' or a matrix. A 'baseline' array is created using the median expression for
#' each gene, and each array is then compared to the baseline array.
#' 
#' 
#' @param object An ExpressionSet or matrix containing log-transformed array
#' data.
#' @param layout A numeric vector, length two. Best results will be obtained if
#' both values are the same, and between 2 and 5 (e.g., c(3,3))
#' @param ... Other arguments that will be passed down to the \code{xyplot}
#' function from the lattice package.
#' @return No output. Used only for the side effect of creating MA plots.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @keywords hplot
#' @export maplot
maplot <- function(object, layout = NULL, ...){
    if(is(object, "ExpressionSet")){
        mat <- exprs(object)
    }else{
        mat <- as.matrix(object)
    }
    if(!is.null(layout) && length(layout) !=2)
        stop("The layout argument should be a numeric vector of length two!\n",
             call. = FALSE)
    if(!is.null(layout) && length(unique(layout)) != 1)
        warning(paste("You will get better results if the layout argument is a vector\n",
                      "of two equal numbers, usually between 2 and 5.\n"),
                call. = FALSE, immediate. = TRUE)
    if(is.null(layout)){
        if(ncol(mat) < 10){
            layout <- c(3,3)
        }else{
            if(ncol(mat) < 17){
                layout <- c(4,4)
            }else{
                layout <- c(5,5)
            }
        }
    }

    med <- apply(mat, 1, median, na.rm = TRUE)
    M <- mat - med
    A <- (mat + med)/2
    df <- data.frame(M = as.vector(M), A = as.vector(A),
                     Id = colnames(mat)[col(mat)])
    print(xyplot(M~A|Id, df, panel = function(x, y, ...){panel.xyplot(x,y,...)
                                                         panel.abline(h = 0, lty = 2)},
                 layout = layout, pch = ".", ...))
}
