\name{sortAscii}
\alias{sortAscii}
%\alias{intersectPresort}
%\alias{setdiffPresort}

\title{
Optimized operations of sets of character-vectors
}
\description{
Sort and use pre-sorted character vectors in set-operations
}
\usage{
sortAscii(x)
%intersectPresort(pop, x)
%setdiffPresort(pop, x)
}
\arguments{
  \item{x}{
	an unsorted vectors of strings
  }
 % \item{pop}{a sorted vector of strings }
}

\details{
By default, sorting is done lexicographically in R. The routine \code{sortAscii} does sorting accoring to the ASCII-order as done in C/C++.

For routines \code{intersectPresort} and \code{setdiffPresort} the first argument has to be sorted according to ASCII-order. This first argument is expected to be large compared to the second argument. Both functions are wrappers for optimized C++-functions performing the set-operation.
}

\value{
An character-vector. In case of \code{intersectPresort} and \code{setdiffPresort}, these vectors are unnamed.
}


\examples{
AA <- matrix( sample( c(LETTERS, letters), 10*30000, rep=TRUE ), ncol=10 )
A <- unique(apply(AA, 1, paste, collapse=""))
B <- sample(AA, 100, replace=FALSE); B <- c(B, "1234")

%res <- intersectPresort( sortAscii(A), B )
}
