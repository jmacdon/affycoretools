###################################################
##
##  Copyright 2005 James W. MacDonald
##
##  MA plotting functions
##
####################################################


maplot <- function(object){
    if(is(object, "ExpressionSet")){
        mat <- exprs(object)
    }else{
        mat <- as.matrix(object)
    }
    if(ncol(mat) < 10){
        layout <- c(3,3)
    }else{
        if(ncol(mat) < 17){
            layout <- c(4,4)
        }else{
            layout <- c(5,5)
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
