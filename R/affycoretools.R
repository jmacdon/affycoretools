#' Make repetitive analyses of microarray and RNA-Seq data simpler with affycoretools.
#' 
#' The affycoretools package is primarily intended to make analyses of Affymetrix GeneChip
#' data simpler and more straightforward. There are any number of packages designed for
#' preprocessing or analyzing Affy data, but there are not so many that help streamline
#' the analysis to help create useful output that can be given to collaborators.
#'
#' The affycoretools package is primarily intended to be used as a way to do reproducible
#' research, where the analysis and documentation are all held in a single file, that is then
#' processed by R to create the output data, as well as a nicely formatted pdf that documents
#' the analysis. The affycoretools package can be used with either Sweave or knitr documents,
#' although these days knitr is really the way to go.
#'
#' In addition, affycoretools can be used with either annaffy or ReportingTools to create useful
#' output in HTML or text format to share with your collaborators. However, ReportingTools is being
#' actively developed and maintained, whereas annaffy is not, so the intention is to slowly
#' convert all the functions to primarily use ReportingTools.
#'
#' @importMethodsFrom affy hist mas5calls
#' @importMethodsFrom AnnotationDbi Term get dbconn select mapIds 
#' @importMethodsFrom Biobase annotation exprs featureNames sampleNames write.exprs pData featureData<- AnnotatedDataFrame
#' @importFrom S4Vectors DataFrame
#' @importFrom affy AffyRNAdeg mas5 plotAffyRNAdeg plotDensity ReadAffy rma
#' @importFrom Biobase addVigs2WinMenu assayDataElement
#' @importFrom graphics layout legend par text
#' @importFrom grDevices dev.off rainbow recordPlot replayPlot
#' @importFrom stats median p.adjust prcomp
#' @importFrom utils write.table
#' @importFrom limma topTable vennCounts decideTests romer ids2indices vennDiagram
#' @importFrom edgeR cpm
#' @importFrom GOstats probeSetSummary
#' @importFrom gcrma gcrma
#' @importFrom splines ns
#' @importFrom xtable xtable
#' @importFrom lattice xyplot panel.xyplot panel.abline lattice.options dotplot bwplot trellis.par.get panel.dotplot
#' @importFrom grid grid.newpage pushViewport viewport unit upViewport
#' @importFrom gplots heatmap.2 bluered venn plot.venn
#' @importFrom oligoClasses db
#' @importFrom ReportingTools HTMLReport CSVFile publish finish reporting.theme filename path
#' @importFrom hwriter hwrite hwriteImage hmakeTag closePage openPage
#' @importFrom ggplot2 ggplot geom_point facet_wrap aes
#' @importFrom grDevices pdf png
#' @importFrom graphics axis image plot segments
#' @importFrom stats cor quantile sd
#' @importFrom utils read.table
#' @importFrom RSQLite dbGetQuery dbDisconnect dbListTables
#' @importFrom BiocGenerics plotPCA
#' @import methods
#' @docType package
#' @name affycoretools
#' @author
#' James W. MacDonald \email{jmacdon@@u.washington.edu}
NULL
