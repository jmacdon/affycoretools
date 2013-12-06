#################################################
##
##  Copyright 2005 James W. MacDonald
##
##  Utility functions for limma2annaffy
##
#################################################
getfldfil <- function(){
  fldfilt.ans <- readline("Please type in the new fold change value. ")
  fldfilt <- try(as.numeric(fldfilt.ans))
  if(is.na(fldfilt)){
    cat(paste(fldfilt.ans, "is not a numeric value!\n"))
    cat("Please enter a number.\n")
    getfldfil()
  }
  return(fldfilt)
}



#' Interactive Functions for limma2anaffy
#' 
#' These are internal functions used to filter the number of genes for
#' \code{limma2annaffy} as well as to automagically get required annotation
#' packages for \code{annaffy}.
#' 
#' These functions were added primarily for the error checking they afford. Now
#' the end user can filter the number of genes based on p-value and fold
#' change, and incorrect values entered will not result in unexpected
#' behaviour.
#' 
#' @aliases getans2 getfldfil getpfil
#' @param msg A character string to pass to the end user.
#' @param allowed A character vector of acceptable responses.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @keywords internal
getans2 <- function(msg, allowed = c("a","c")){
  repeat{
    out <- paste("[", paste(allowed, collapse = "/"),
                 "]")
    outmsg <- paste(msg, out)
    cat(outmsg)
    answer <- readLines(n = 1)
    if(answer %in% allowed)
      break
    else cat(paste(answer, "is not a valid response, try again.\n"))
  }
  answer
}

#getans <- function(cdfname){
#  answer <- readline("Do you want the release or devel version (r/d)?\n")
#  if(substr(answer,1,1) == "r" | substr(answer,1,1) == "R"){
#    do.call("install.packages2",list(cdfname, develOK=FALSE))
#  }
#  if(substr(answer,1,1) == "d" | substr(answer,1,1) == "D"){
#    do.call("install.packages2", list(cdfname, develOK=TRUE))
#  }
#  if(substr(answer,1,1) %in% c("d","D","r","R") == FALSE){
#    cat(paste(substr(answer,1,1), "is not one of the available choices!\n"))
#    cat("Please choose either (r)elease or (d)evel.\n")
#    getans()
#  }
#}

getpfil <- function(){
  pfilt.ans <- readline("Please type in the new p-value. ")
  pfilt <- as.numeric(pfilt.ans)
  if(is.na(pfilt)){
    cat(paste(pfilt.ans, "is not a numeric value!\n"))
    cat("Please enter a number between 0 and 1.\n")
    getpfil()
  }

  if(pfilt > 1 | pfilt < 0){
    cat(paste(pfilt, "is not between 0 and 1!\n"))
    cat("Please enter a number between 0 and 1.\n")
    getpfil()
  }
  return(pfilt)
}
