\name{getMainProbes}
\alias{getMainProbes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Remove control probesets from ST arrays}
\description{
  This function is designed to remove all but the 'main' type of
  probesets from the Gene ST array types.
}
\usage{
getMainProbes(pdinfo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{pdinfo}{Character. The name of a pdinfo package, used by oligo to summarize
  probesets. This will generally be of the form e.g.,
  pd.hugene.1.0.st.v1}
 }

\value{
  Returns a boolean vector that can be used to subset e.g., an
  ExpressionSet of Gene ST data, or an MArrayLM object created from Gene
  ST data.
}

\author{James W. MacDonald <jmacdon@u.washington.edu>}
\keyword{manip}% __ONLY ONE__ keyword per line
