\name{geecc-package}
\alias{geecc-package}
\alias{geecc}
\docType{package}
\title{
Gene set enrichment for two or three categories
}
\description{
This package performs gene set enrichment analyses considering two or three categories. Categories might be regulated genes, sequence length, GC content, GO terms, KEGG pathways and so on.
}

\author{
Markus Boenn
Maintainer: Markus Boenn <markus.boenn@ufz.de>
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
ConCub.obj.3 <- filterConCub( obj=ConCub.obj.2, filter=ConCubFilter.obj )
plotConCub( obj=ConCub.obj.3, filter=ConCubFilter.obj )
x <- getTable(ConCub.obj.3)
}

\keyword{ package }
