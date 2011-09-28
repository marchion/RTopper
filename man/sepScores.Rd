\name{sepScores}
%\Rdversion{1.1}

\alias{sepScores}

\docType{data}

\title{
  A list of separate gene-to-phenotype association scores,
  obtained indipendently for each distinct data set
  to be used to run the examples in the RTopper package
}

\description{
  A \code{list} containing distinct named \code{numeric} vectors
  corresponding to the gene-to-phenotype association scores
  resulting from the sepatate analysis of each data set.
  These separate gene-to-phenotype scores are compurted
  by \code{\link{computeDrStat}} function.
  Can be used as input to \code{\link{runBatchGSE}}.
}

\usage{ data(sepScores) }

\format{
  This object is a \code{list} of length four, one element for data set:
  
  \code{"dat.affy"}{: a \code{numeric} vector of length 500,
    corresponding to the sepatate phenotype association scores
    computed for Affymetrix gene expression data;}
  
  \code{"dat.agilent"}{: a \code{numeric} vector of length 500,
    corresponding to the separate phenotype association scores
    computed for Agilent gene expression data;}
  
  \code{"dat.cnvHarvard"}{: a \code{numeric} vector of length 500,
    corresponding to the separate phenotype association scores
    computed for Harvard CNV data;}

  \code{"dat.cnvMskcc"}{: a \code{numeric} vector of length 500,
    corresponding to the separate phenotype association scores
    computed for the MSKCC CNV data;}
  
}


\source{
  Computed using the \code{\link{computeDrStat}} function
  from the TCGA data contained in \code{\link{data}}.
}

\examples{
data(sepScores)
class(sepScores)
names(sepScores)
str(sepScores)
}

\keyword{datasets}
