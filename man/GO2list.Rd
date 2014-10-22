\name{GO2list}
\alias{GO2list}
\alias{KEGG2list}

\title{
Filter GO and KEGG database
}
\description{
Filter GO and KEGG database and transforms database to list
}
\usage{
GO2list(dbase, go.cat = NULL, rm = NULL, keep = NULL)
KEGG2list(dbase, rm = NULL, keep = NULL)
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
  \item{go.cat}{GO category ("MF", "BP", "CC") that should be returned and filtered}
  \item{rm}{remove these terms}
  \item{keep}{keep only these terms}
}

\details{
The settings for \dQuote{rm} and \dQuote{keep} can be combined, allowing for efficient reduction of the number of GO terms and KEGG pathways, respectively.

Providing a named list instead of a database can be useful for non-model organisms, where only a draft Blast2GO-annotation is available. In this case, the names of the list are the GO terms (or KEGG pathways) and the content of each list item is a character vector with tag-ids.
}

\value{
A named list with each slot containing the ids for the term or pathway.
}


\examples{
library(hgu133plus2.db)
x <- GO2list(dbase=hgu133plus2GO2PROBE, go.cat="CC",
	rm=c("GO:0000139", "GO:0000790", "GO:0005730", "GO:0005739"))
}
