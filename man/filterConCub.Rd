\name{filterConCub}
\alias{filterConCub}

\title{
Filter results from two- or three-way tests
}
\description{
Performs filtering on results from two- or three-way tests
}
\usage{
filterConCub(obj, filter, p.adjust.method = "none")
}
\arguments{
  \item{obj}{
  object of type \code{concub}
}
  \item{filter}{
  object of type \code{concubfilter}
}
  \item{p.adjust.method}{
set adjustment-method for p-values. Must match any of \code{\link{p.adjust.methods}}.
}
}


\details{
You have to execute \code{\link{runConCub}} before filtering.
}


\value{
an (extended) object of type \code{concub} with filtered test results
}

\examples{
# a character vector listing possible
# adjustment approaches
p.adjust.methods
}