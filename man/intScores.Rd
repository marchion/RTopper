\name{intScores}
%\Rdversion{1.1}

\alias{intScores}

\docType{data}

\title{
  A list of genomic scores integrated across distinct data sets
  to be used to run the examples in the RTopper package
}

\description{
  A \code{list} containing a named \code{numeric} vector
  corresponding to the genomic score resulting from the integration
  across distinct data set.
  These integrated gene-to-phenotype scores are computed
  by the \code{\link{computeDrStat}} function.
  Can be used as input to \code{\link{runBatchGSE}}.
}

\usage{ data(intScores) }

\format{
  This object is a \code{list} of length one:
  
  \code{"integrated"}{: a \code{numeric} vector of length 500,
    corresponding to the integrated phenotype association scores
    computed for each of the 500 genes used in the examples;}
  
}

\source{
  Computed using the \code{\link{computeDrStat}} function
  from the TCGA data contained in \code{\link{dat}} and \code{\link{pheno}}.
}

\examples{
data(intScores)
class(intScores)
names(intScores)
str(intScores)
}

\keyword{datasets}
