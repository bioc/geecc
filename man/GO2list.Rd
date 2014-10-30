\name{GO2list}
\alias{GO2list}
\alias{KEGG2list}
\alias{GO2level}
\alias{GO2offspring}

\title{
Filter GO and KEGG database
}
\description{
Filter GO and KEGG database and transform database to list
}
\usage{
GO2list(dbase, go.cat = NULL, rm = NULL, keep = NULL)
KEGG2list(dbase, rm = NULL, keep = NULL)
GO2offspring(x)
GO2level(x, go.level=-1)
}
\arguments{
  \item{dbase}{
	A datastructure storing identifieres of GO/KEGG terms and assigned genes. Can be one of
  	\describe{
	\item{database}{ usually of class \sQuote{ProbeGo3AnnDbBimap} (as defined in package \dQuote{AnnotationDbi}) }
	\item{named list}{ with keys being the identifiers and values being genes }
	\item{dataframe}{ with first column being the identifiers and second column being genes. Additional columns are ignored. }
	}
  }
  \item{x}{a list with keys being the identifiers and values being genes (e.g. output of \code{GO2list})}
  \item{go.cat}{GO category ("MF", "BP", "CC") that should be returned and filtered}
  \item{go.level}{Level in the DAG of GO terms. Defaults to \dQuote{-1} for pass through without modification. Otherwise: a positive integer giving the level at which GO terms should be grouped together.}
  \item{rm}{remove these terms}
  \item{keep}{keep only these terms}
}

\details{
The settings for \dQuote{rm} and \dQuote{keep} can be combined, allowing for efficient reduction of the number of GO terms and KEGG pathways, respectively.

Providing a named list instead of a database can be useful for non-model organisms, where only a draft Blast2GO-annotation is available. In this case, the names of the list are the GO terms (or KEGG pathways) and the content of each list item is a character vector with tag-ids.

The function \code{GO2offspring} does the same as the databaseGO2ALLPROBES function does (e.g. hgu133plus2GO2ALLPROBES). I.e. instead of representing only features (probe sets, genes, ...) assigned to the GO terms directly, it also contains all features assigned to all children (offsprings).

The function \code{GO2level} groups GO terms together at a more general level to simplify data interpretation and speed up runtime. This function works according to the level option provided by DAVID, but the number of levels is not restricted.
}

\value{
A named list with each slot containing the ids for the term or pathway.
}


\examples{
library(hgu133plus2.db)
x <- GO2list(dbase=hgu133plus2GO2PROBE, go.cat="CC",
	rm=c("GO:0000139", "GO:0000790", "GO:0005730", "GO:0005739"))
}
