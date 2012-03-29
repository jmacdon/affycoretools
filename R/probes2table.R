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
