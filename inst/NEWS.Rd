\name{NEWS}
\title{News for Package 'geecc'}

\section{Changes in geecc version 1.1.11 (2016-08-xx)}{
  \itemize{
        \item major revision of \code{init} and \code{runConCub} resulting in significant reduction of memory allocation and runtime
    \item speed-up of set-operations
        \item corrected typos in some messages
  }
}


\section{Changes in geecc version 1.1.8 (2014-12-31)}{
  \itemize{
    \item explicite garbage collection avoids memory overflow of function 'runConCub'
    \item function \code{GO2level} supports definition of relations to be considered
  }
}


\section{Changes in geecc version 1.1.6 (2014-12-17)}{
  \itemize{
    \item extended 'show' for concub-objects
    \item bug fix for concub-objects that occured when 'population' was a subset of the categories
    \item originally private function 'pval2star' is now public

  }
}


\section{Changes in geecc version 1.1.5 (2014-11-12)}{
  \itemize{
    \item fixed bug in output of getTable (concub-object)
  }
}


\section{Changes in geecc version 1.1.4 (2014-10-30)}{
  \itemize{
    \item fixed typo in GO2level
    \item added function GO2offspring, which corresponds to _database_GO2ALLPROBES
  }
}


\section{Changes in geecc version 1.1.3 (2014-10-28)}{
  \itemize{
    \item bug fixed in GO2level
  }
}


\section{Changes in geecc version 1.1.2 (2014-10-27)}{
  \itemize{
    \item added function GO2level, which allows to group GO terms at a more general level
  }
}


\section{Changes in geecc version 1.1.1 (2014-10-24)}{
  \itemize{
    \item fixed bug in filterConCub that caused rearrangement of categories in plotConCub
    \item reduced required amount of memory
    \item functions GO2list and KEGG2list now accept dataframes
  }
}

