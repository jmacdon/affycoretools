\name{pca.legend}
\alias{pca.legend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A Function to Automagically Place a Legend in a PCA Plot }
\description{
 This function places a legend in a PCA plot depending on which corner
 is available. If there are no available corners, user intervention will
 be required.
}
\usage{
pca.legend(pca, groupnames, groups, x.coord = NULL, y.coord = NULL,
saveup = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pca}{A 'pca' object }
  \item{groupnames}{A vector of names describing the different groups }
  \item{groups}{A numeric vector delineating the group memberships }
  \item{x.coord}{X-coordinate for legend. Used if automatic placement
    will fail}
  \item{y.coord}{Y-coordinate for legend. Used if automatic placement
    will fail}
  \item{saveup}{Boolean: Save a small value for plotting additional text?}
}

\value{
 This function returns nothing. It is used only for the side effect of
 placing a legend in a plot.
}
\author{James W. MacDonald}

\keyword{internal }% at least one, from doc/KEYWORDS
