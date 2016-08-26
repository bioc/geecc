\name{plotConCub}
\alias{plotConCub}
\title{
Generate a heatmap showing \eqn{log_2}{log2} odds ratios and \eqn{P}{P}-values.
}
\description{
The function generates a heatmap by calling the \code{heatmap.2}-function from the \code{gplots}-package. Each cell shows the \eqn{log_2}{log2} odds ratio of the test for the corresponding variables. In addition, stars indicate the \eqn{P}{P}-value for this test.
}
\usage{
plotConCub(obj, filter, fix.cat = 1, show=list(), dontshow=list(),
	args_heatmap.2 = list(), col = list(range = NULL),
	alt.names = list(), t = FALSE)
}
\arguments{
  \item{obj}{ An object with class concub }
  \item{filter}{ An object with class concubfilter }
  \item{fix.cat}{ The heatmap can only visualize a two-dimensional table. In case of three-dimensions, one dimension (category) must be fixed. }
  \item{show}{
A named list. The names are the names of the categories. Each item is a character vector of variables that should be shown in the plot.
}
  \item{dontshow}{
A named list. The names are the names of the categories. Each item is a character vector of variables that should not be shown in the plot.
}
  \item{args_heatmap.2}{ Arguments passed to \sQuote{\code{\link{heatmap.2}}}. Can be used to change size of fonts etc. }
  \item{col}{ A vector of colors, for instance from \code{\link{heat.colors}} }
  \item{alt.names}{ Substitute variables by alternative terms. For instance, if variables are artificial ids, they can be substituted by descriptive text for the heatmap. }
 \item{t}{ logical; transpose matrix for heatmap. Default \code{FALSE}. }
}

\value{An object with class concub.}

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
ConCub.obj <- new("concub", categories=CatList)
ConCub.obj.2 <- runConCub( obj=ConCub.obj, filter=ConCubFilter.obj, nthreads=1 )
ConCub.obj.3 <- filterConCub( obj=ConCub.obj.2, filter=ConCubFilter.obj )
plotConCub( obj=ConCub.obj.3, filter=ConCubFilter.obj )

}
