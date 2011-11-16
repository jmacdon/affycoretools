###################################################
##
##  Copyright 2005 James W. MacDonald
##
##  MA plotting functions
##
####################################################


maplot <- function(object, layout = NULL){
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
                 layout = layout, pch = "."))
}
