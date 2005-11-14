\name{tableFilt}
\alias{tableFilt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Filter a topTable object}
\description{
  This function is designed to filter genes from a \code{topTable}
  object based on p-value and/or fold change. This is an internal
  function and is not intended to be called by thte end user.
}
\usage{
tableFilt(fit, coef = 1,  number = 30, fldfilt = NULL, pfilt = NULL,
adjust = "fdr")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{ An \code{MArrayLM} object, resulting from a call to \code{eBayes}}
  \item{coef}{The contrast to be extracted into the topTable. See
    ?topTable for more information. }
  \item{number}{ The number of genes to output. Only used if both
    foldfilt and pfilt are NULL. }
  \item{fldfilt}{The absolute value of fold difference to filter
    on. This assumes the data are log transformed.}
  \item{pfilt}{The p-value to filter on.}
  \item{adjust}{The multiplicity adjustment to use. Options are
    '"bonferroni"', '"holm"', '"hochberg"', '"hommel"', '"fdr"' and
    '"none"'. If '"none"' then the p-values are not adjusted. A 'NULL'
    value will result in the default adjustment method, which is
    '"fdr"'.}
  }

\value{
  Returns a \code{data.frame} containing the selected genes.
}

\author{James W. MacDonald <jmacdon@med.umich.edu> }

\keyword{internal }% __ONLY ONE__ keyword per line
