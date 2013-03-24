\name{vennPage}
\alias{vennPage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{High-level function for making Venn diagrams with clickable links
  to HTML pages with the underlying genes.}
\description{
  This function is designed to be used in conjunction with the
  \code{makeVenn} function, to first create a set of HTML pages
  containing the genes that are represented by the cells of a Venn
  diagram, and then create an HTML page with the same Venn diagrams,
  with clickable links that will point the end user to the HTML pages.
}
\usage{
vennPage(vennlst, pagename, pagetitle, cex.venn = 1, shift.title = FALSE,
                     baseUrl = ".", reportDirectory = NULL)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vennlst}{The output from \code{makeVenn}.}
  \item{pagename}{Character. The file name for the resulting HTML
    page. Something like 'venns' is reasonable. Note that the .html will
  automatically be appended.}
  \item{pagetitle}{Character. The heading for the HTML page.}
  \item{cex.venn}{Numeric. Adjusts the size of the font in the Venn
    diagram. Usually the default is OK.}
  \item{shift.title}{Boolean. Should the right contrast name of the Venn diagram
    be shifted down? Useful for long contrast names. If a two-way Venn
    diagram, this will shift the right name down so they don't
    overlap. If a three-way Venn diagram, this will shift the top right
    name down.}
  \item{baseUrl}{Character. The base URL for the resulting HTML
    page. The default of "." is usually optimal.}
  \item{reportDirectory}{If \code{NULL}, the reportDirectory will be
    extracted from the vennlst. This is usually what one should do.}
}

\details{
  This function is intended to be used as part of a pipeline, by first
  calling \code{makeVenn} and then using the output from that function
  as input to this function to create the HTML page with clickable links.
}
\value{
  An HTMLReport object. If used as input to the ReportingTools
  \code{publish} function, this will create a link on an index page to
  the Venn diagram HTML page. See e.g., the microarray analysis vignette
  for ReportingTools for more information.
}
\author{James W. MacDonald <jmacdon@u.washington.edu>}
\examples{
  \dontrun{
    mat <- matrix(rnorm(1e6), ncol = 20)
    design <- model.matrix(~factor(1:4, each=5))
    colnames(design) <- LETTERS[1:4]
    contrast <- matrix(c(1,-1,0,0,1,0,-1,0,1,0,0,-1,0,1,-1,0,0,1,0,-1),
    ncol = 5)
    colnames(contrast) <- paste(LETTERS[c(1,1,1,2,2)],
    LETTERS[c(2,3,4,3,4)], sep = " vs ")
    fit <- lmFit(mat, design)
    fit2 <- contrasts.fit(fit, contrast)
    fit2 <- eBayes(fit2)
    ## two Venn diagrams - a 3-way Venn with the first three contrasts
    ## and a 2-way Venn with the last two contrasts
    collist <- list(1:3,4:5)
    venn <- makeVenn(fit2, contrast, design, eset, collist = collist)
    vennreport <- vennPage(venn, "index.html", "Venn diagrams")
    indexPage <- HTMLReport("index", "My results", reportDirectory =
    ".", baseUrl = ".")
    publish(vennreport)
    finish(indexPage)
    }
}

\keyword{ manip}% at least one, from doc/KEYWORDS

  
