\name{concub-class}
\Rdversion{1.1}
\docType{class}
\alias{concub-class}
\alias{getTable}
\alias{getTable,concub-method}
\alias{show,concub-method}
\alias{initialize,concub-method}


\title{Class \code{"concub"}}
\description{
An object of type concub
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("concub", ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{fact}:}{ A named list of named lists. Each item of the outer list represents the two or three categories. Each item of the inner lists represents a variable of the category. }
    \item{\code{population}:}{ A character vector containig the background population. As an alternative, an object with class \sQuote{eSet}, \sQuote{ExpressionSet}, or \sQuote{DGEList}; background population is then set to the outcome of \code{rownames(population)}. }
    \item{\code{keep.empty.vars}:}{ A boolean list with names being names of \code{fact}. }
    \item{\code{options}:}{ Additional options for individual categories. See Details. }
    \item{\code{approx}:}{ specifies the minimum expected value when an exact hypergeometric test (below) or the chi-squared approximation should be used. Defaults to 0. }
    \item{\code{null.model}:}{ A formula specifying the null-model of the test.  }
    \item{\code{test.result}:}{ A list to store test results. Filled up after \code{runConCub}. }
    \item{\code{test.result.filter}:}{ A list to store filtered test results. Filled up after \code{filterConCub}. }
    \item{\code{test.result.filter.heatmap}:}{ A list to store heatmaps for further manipulation. Filled up after \code{plotConCub}. }
  }
 The last three slots are not set by the user.
}

\details{
	Specifying the background population is crucial for the tests for association between factors. Usually the population is the set of all probe sets represented on a micro-array or the set of all genes in a genome. If an expression set is passed with the population-parameter, all probe sets beginning with the pattern \dQuote{AFFX} (Affymetrix quality control) are removed.
}

\section{Methods}{
	\describe{
		\item{getTable}{\code{signature(object = "concub")}:
			creates a table containing the results of the tests.

			Usage: \code{getTable(object, na.rm=TRUE, dontshow=list())}

			Arguments:
			\describe{
			\item{object}{an object of type concub}
			\item{na.rm}{logical. If \code{TRUE} (the default), rows with \code{NA} P-values (or odds ratios) are removed}
			}

			The resulting table is a data frame with 8 or 10 columns, depending on if a two- or three-way test was applied
			\enumerate{
				\item{'factor1':}{this column has the name of the first category}
				\item{'factor2':}{this column has the name of the second category}
				\item{'factor3':}{this column has the name of the third category}
				\item{n.'factor1':}{number of items in variable of first category}
				\item{n.'factor2':}{number of items in variable of second category}
				\item{n.'factor3':}{number of items in variable of third category}
				\item{p.value:}{(probably adjusted) p-value}
				\item{log2.odds.ratio:}{log2 of the sample odds ratio}
				\item{n.tags:}{number of items at position \eqn{x_{1,1,1}}{x_{1,1,1}} or \eqn{x_{1,1}}{x_{1,1}}}
				\item{tags:}{items at position \eqn{x_{1,1,1}}{x_{1,1,1}} or \eqn{x_{1,1}}{x_{1,1}} (e.g. gene identifiers)}
			}
			You have to run \code{filterConCub()} before you can get the table. If \code{filterConCub()} was not run, a warning is shown and \code{getTable} returns \code{NULL}.
		}

		\item{show}{\code{signature(object = "concub")} }
	}
}


\seealso{
\code{GOStats} to perform a simple two-way enrichment analysis
}
\examples{
showClass("concub")
}
\keyword{classes}
