\name{pval2star}
\alias{pval2star}

\title{
Transform P-values to stars}
\description{
Transform P-values to stars}

\usage{
pval2star(x)
}
\arguments{
  \item{x}{ A matrix of P-values  }
}

\details{
Use stars as simplification of P-values
}

\value{
A character matrix of same dimension and names as \code{x} with stars instead of P-values.
}


\examples{
x <- matrix( runif(25), nrow=5, dimnames=list(LETTERS[1:5], letters[1:5]) )
pval2star(x)
}
