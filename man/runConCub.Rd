\name{runConCub}
\alias{runConCub}

\title{
Enrichment analysis on two- or three-way contingency tables.
}
\description{
Perform the enrichment analysis on two- or three-way contingency tables.
}
\usage{
runConCub(obj, filter, nthreads = 2, rng = NULL,
	verbose=list(output.step=0, show.cat1=FALSE,
	show.cat2=FALSE, show.cat3=FALSE))
}
\arguments{
  \item{obj}{an object with class concub}
  \item{filter}{an object with class concubfilter}
  \item{nthreads}{number of threads to use in \code{hypergeom.test}}
  \item{rng}{a named list. Restrict enrichment analysis to these category variables}
  \item{verbose}{
	A list to control verbosity:
		\describe{
			\item{output.step:}{ after how many variables passed of category 2 a control output should be printed }
			\item{show.cat1:}{ show current level of category 1 }
			\item{show.cat2:}{ show current level of category 2 }
			\item{show.cat3:}{ show current level of category 3 }
		}
	}
}

\details{
	This function applies a test for association for all combinations of all variables of all categories to be tested. Depending on the settings in the concubfilter-object, a one-sided or two-sided test is made, using the exact hypergeometric test as implemented in the hypergea-package if the smallest expected value is smaller than 5, or using the chi-squared test as implemented in the loglm-function implemented in the MASS-package. The minimum expected value can be changed in the concub-object by the user (\code{approx}-parameter).
	In this function only filter-settings those filter settings are used, which skip the tests.
}
\examples{
##
## a completely artificial example run
## through the routines of the package
##
R <- 500
#generate R random gene-ids
ID <- sapply(1:R, function(r){paste( sample(LETTERS, 10), collapse="" ) } )
ID <- unique(ID)

#assign artificial differentially expressed genes randomly
category1 <- list( deg.smallFC=sample(ID, 100, rep=FALSE),
	deg.hughFC=sample(ID, 100, rep=FALSE) )
#assign artificial GO terms of genes randomly
category2 <- list( go1=sample(ID, 50, replace=FALSE),
	go2=sample(ID, 166, replace=FALSE),
	go3=sample(ID, 74, replace=FALSE),
	go4=sample(ID, 68, replace=FALSE) )
#assign artificial sequence length of genes randomly
LEN <- setNames(sample(seq(100, 1000, 100), length(ID), replace=TRUE), ID)
category3 <- split( ID, f=factor(LEN, levels=seq(100, 1000, 100)) )
CatList <- list(deg=category1, go=category2, len=category3)

ConCubFilter.obj <- new("concubfilter", names=names(CatList))
ConCub.obj <- new("concub", fact=CatList)
ConCub.obj.2 <- runConCub( obj=ConCub.obj, filter=ConCubFilter.obj, nthreads=1 )
ConCub.obj.2
}
