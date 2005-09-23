#####################################################3
##
##  Copyright 2005 James W. MacDonald
##
##  probes2table - convert a vector of affy probe IDs
##                 to HTML or text tables
##
####################################################

probes2table <- function(eset, probids, lib, FC=NULL, tstat=NULL, pval=NULL,
                         anncols=aaf.handler()[c(1:3, 7:8, 10:13)], html=TRUE,
                         text=FALSE, express = TRUE, save=FALSE,  filename){

  require(annaffy, quietly=TRUE)

 

  anntable <- aafTableAnn(probids, lib, anncols)
  if(!is.null(tstat))
    testtable <- aafTable("t-statistic" = round(tstat, 2))
  if(!is.null(pval)){
    if(exists("testtable")){
      testtable <- merge(testtable, aafTable("p-value" = round(pval,3)))
    }else{
      testtable <- aafTable("p-value" = round(pval, 3))
    }
  }
  
  if(!is.null(FC)){
    if(exists("testtable")){
      testtable <- merge(testtable, aafTable("Fold change" = round(FC,2)))
    }else{
      testtable <- aafTable("Fold change" = round(FC, 2))
    }
  }
  
  if(exists("testtable"))
    anntable <- merge(anntable, testtable)
  
  if(express){
    exprtable <- aafTableInt(eset, probeids=probids)
    anntable <- merge(anntable, exprtable)
  }
  
  if(html)
    saveHTML(anntable, paste(filename, "html", sep="."), title=filename)
  if(text)
    saveText(anntable, paste(filename, "txt", sep="."), header=TRUE)
  if(save)
    return(anntable)
}
