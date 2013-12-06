#####################################################3
##
##  Copyright 2005 James W. MacDonald
##
##  probes2table - convert a vector of affy probe IDs
##                 to HTML or text tables
##
##  11-08-2005 Removed FC, tstat, and pval in lieu of
##      a list (otherdata) containing arbitrary information to add to table
##
####################################################



#' Convert Affy Probe ids to Annotated HTML Table
#' 
#' A function to convert a vector of Affy ids to an annotated HTML table.
#' 
#' 
#' @param eset An \code{ExpressionSet} containing Affy expression values.
#' @param probids A vector of probe ids.
#' @param lib An annotation package for the Affy chips used.
#' @param otherdata A *named* list of additional information to include in the
#' resulting table. Examples would be t-statistics, p-values, fold change, etc.
#' Each list item should be a vector the same length as the probids vector. The
#' name associated with each list item will be used as the column name in the
#' resulting table.
#' @param anncols A vector of things to annotate, produced by a call to
#' aaf.handler().
#' @param html Output data in HTML tables? Defaults to \code{TRUE}.
#' @param text Output data in text tables? Defaults to \code{TRUE}.
#' @param express Output expression values in table?  Defaults to \code{TRUE}.
#' @param save Should tables be saved as R objects for further processing?
#' Defaults to \code{FALSE}.
#' @param filename Filename of the resulting HTML table.
#' @param title Title for HTML table. If \code{NULL}, the filename will be
#' used.
#' @return If \code{save} is \code{TRUE}, a \code{data.frame} is saved
#' containing the data.
#' @author James W. MacDonald <jmacdon@@u.washington.edu>
#' @seealso \code{\link[limma:toptable]{topTable}},
#' \code{\link[annaffy]{aafTableAnn}}
#' @keywords manip
#' @export probes2table
probes2table <- function(eset, probids, lib, otherdata = NULL,
                         anncols=aaf.handler()[c(1:3, 6:7, 9:12)], html=TRUE,
                         text=FALSE, express = TRUE, save=FALSE,  filename,
                         title = NULL){


  ## test that lib has a .db extension
  if(length(grep("\\.db$", lib)) < 1)
      lib <- paste(lib, "db", sep = ".")
  anntable <- aafTableAnn(probids, lib, anncols)

  if(!is.null(otherdata)){
    if(!is(otherdata, "list")) stop("otherdata must be a list!")
    if(is.null(names(otherdata))) stop("otherdata must be a *named* list!")
    testtable <- aafTable(items = otherdata)
  }
  if(exists("testtable"))
    anntable <- merge(anntable, testtable)
  
  if(express){
    exprtable <- aafTableInt(eset, probeids=probids)
    anntable <- merge(anntable, exprtable)
  }
  if(is.null(title))
      title <- filename
  
  if(html)
    saveHTML(anntable, paste(filename, "html", sep="."), title=title)
  if(text)
    saveText(anntable, paste(filename, "txt", sep="."), header=TRUE)
  if(save)
    return(anntable)
}
