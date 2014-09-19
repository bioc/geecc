\name{concubfilter-class}
\Rdversion{1.1}
\docType{class}
\alias{concubfilter-class}
\alias{show,concubfilter-method}
\alias{initialize,concubfilter-method}

\alias{p.value}
\alias{p.value<-}
\alias{p.value,concubfilter-method}
\alias{p.value<-,concubfilter-method}

\alias{minimum.l2or}
\alias{minimum.l2or<-}
\alias{minimum.l2or,concubfilter-method}
\alias{minimum.l2or<-,concubfilter-method}

\alias{test.direction}
\alias{test.direction<-}
\alias{test.direction,concubfilter-method}
\alias{test.direction<-,concubfilter-method}

\alias{skip.min.group}
\alias{skip.min.group<-}
\alias{skip.min.group,concubfilter-method}
\alias{skip.min.group<-,concubfilter-method}

\alias{skip.min.obs}
\alias{skip.min.obs<-}
\alias{skip.min.obs,concubfilter-method}
\alias{skip.min.obs<-,concubfilter-method}

\alias{skip.zeroobs}
\alias{skip.zeroobs<-}
\alias{skip.zeroobs,concubfilter-method}
\alias{skip.zeroobs<-,concubfilter-method}

\alias{drop.insignif.layer}
\alias{drop.insignif.layer<-}
\alias{drop.insignif.layer,concubfilter-method}
\alias{drop.insignif.layer<-,concubfilter-method}

\alias{drop.lowl2or.layer}
\alias{drop.lowl2or.layer<-}
\alias{drop.lowl2or.layer,concubfilter-method}
\alias{drop.lowl2or.layer<-,concubfilter-method}

\alias{drop.wrongdir.layer}
\alias{drop.wrongdir.layer<-}
\alias{drop.wrongdir.layer,concubfilter-method}
\alias{drop.wrongdir.layer<-,concubfilter-method}


\title{Class \code{"concubfilter"}}
\description{
An object of type concubfilter
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("concubfilter", ...)}.
}
\section{Slots}{
	For a more detailed description of some slots see below.
	\describe{
		\item{\code{nfact}:}{ a numeric giving the number of factors }
		\item{\code{names}:}{ a character vector giving the name of each factor }
		\item{\code{p.value}:}{ a numeric giving the P-value to be considered. Defaults to 0.1. }
		\item{\code{test.direction}:}{ a character giving the direction of association. Defaults to \dQuote{two.sided}. }
		\item{\code{minimum.l2or}:}{ a numeric giving the minimum absolute log2 odds ratio to be considered. Defaults to 0. }

		\item{\code{skip.min.group}:}{ a numeric giving the minimum number of tag a group is allowed to have. Defaults to 2. }
		\item{\code{skip.min.obs}:}{ a numeric giving the minimum number at the position of interest allowed. Defaults to 2. }
		\item{\code{skip.zeroobs}:}{ a logical. Defaults to TRUE. }

		\item{\code{drop.insignif.layer}:}{ A vector of logicals. By default, all positions are set to FALSE. }
		\item{\code{drop.wrongdir.layer}:}{ A vector of logicals. By default, all positions are set to FALSE. }
		\item{\code{drop.lowl2or.layer}:}{ A vector of logicals. By default, all positions are set to FALSE. }
	}
}

\details{
	The large number of different filter options (and corresponding getter and setter accessors) makes it necessary to maintain them in a special class. This differs from other packages like GOstats, where arguments for controlling the program and the results are stored in the same object.
}

\section{Methods and slot accessors}{
Several methods are implemented for class concubfilter. They can be roughly grouped into informative, basic, skip-test, and data-reduction methods.

Individual options can be accessed by the corresponding getter and setter methods, for instance
	\describe{
		\item{skip.zeroobs}{\code{signature(x = "concubfilter")}: get current setting to skip test in case of zero cell }
		\item{skip.zeroobs<-}{\code{signature(x = "concubfilter")}: set a new value to skip test in case of zero cell }
	}

\subsection{Informative methods}{
Currently only a single method is implemented.
	\describe{
		\item{show}{\code{signature(object = "concubfilter")}: a short summary about current filter settings }
	}
}

\subsection{Basic filters/threshoulds}{
	\describe{
		\item{p.value}{ The maximum P-value that should be taken into account }
		\item{minimum.l2or}{ The minumum absolute of log2 odds ratio that should be taken into account }
		\item{test.direction}{ The direction of the association. Can be \dQuote{two.sided}, \dQuote{greater} (test for over-representation), or \dQuote{less} (test for under-representation) }
	}
}


\subsection{Skip test}{
The following filters cause a skip of a test, i.e. the test is never run if at least one of the conditions is fullfilled.
	\describe{
		\item{skip.zeroobs}{ skip the test, if the position of interest (\eqn{x_{00}} or \eqn{x_{000}}) is zero, i.e. no tag from the population matches the conditions defined at the marginals. }
		\item{skip.mingroup}{ skip the test, if the groups considered at the position of interest are too small at all. }
		\item{skip.minobs}{ skip the test, if the position of interest contains less than \sQuote{this value} entries. }
	}
}

\subsection{Data reduction}{
The following filters reduce the amout of outcome of the tests. They are applied to both, the (2 or 3 dimensional) table containing the odds ratios and the table containing the corresponding P-values.
	\describe{
		\item{drop.insignif.layer}{ drop all layers in the tables where all \eqn{P}{P}-values are greater than the value defined in \sQuote{p.value}. }
		\item{drop.wrongdir.layer}{ drop all layers in the tables where all odds ratios are showing into the opposite direction as defined in  \sQuote{test.direction}. }
		\item{drop.lowl2or.layer}{ drop all layers in the tables where all absolutes of the log2 odds ratios are smaller than \sQuote{min.l2or}. }
	}
}

} % end of section methods

\examples{
showClass("concub")
}
\keyword{classes}
