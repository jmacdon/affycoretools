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
#' @return No output. Used only for the side effect of creating MA plots.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @keywords hplot
#' @export maplot
maplot <- function(object){
    if(is(object, "ExpressionSet")){
        mat <- exprs(object)
    }else{
        mat <- as.matrix(object)
    }
    med <- apply(mat, 1, median, na.rm = TRUE)
    M <- mat - med
    A <- (mat + med)/2
    df <- data.frame(M = as.vector(M), A = as.vector(A),
                     Id = colnames(mat)[col(mat)])
    g <- ggplot(df, aes(A,M)) + geom_point(size = 0.05) + facet_wrap(~Id)
    g
}
