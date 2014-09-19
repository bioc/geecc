\name{marioni}
\docType{data}
\alias{marioni}
\title{Affymetrix microarray gene expression data}
\description{The experiment aims to detect differentially expressed genes in Affymetrix micro-arrays and RNA-seq data in a comparative study. For this, samples from two tissues (liver and kidney) were compared.}
\usage{marioni}
\format{A \code{data.frame} containing gene expression values from an Affymetrix microarray, including P-values, log2-fold changes and alternative annotations}
\source{http://giladlab.uchicago.edu/data.html}
\references{
 Marioni, J. C. et al. (2008) RNA-seq: an assessment of technical reproducibility and comparison with gene expression arrays. Genome research
}
\examples{
data(marioni)
head(marioni[, 1:5])
}

\keyword{datasets}
