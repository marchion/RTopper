\name{gseResultsSep}
%\Rdversion{1.1}

\alias{gseResultsSep}

\docType{data}

\title{
  A list of separated gene set enrichment p-values
  to be used to run the examples in the RTopper package
}

\description{
  A \code{list} containing distinct named \code{numeric} vectors
  corresponding to the gene set enrichment p-value
  separately computed with \code{\link{runBatchGSE}}
  for each distinct data set.
  Note that for each data set there are two set of p-values,
  one for GO and one for KEGG.
  These separate p-values can be combined across data sets
  by the \code{\link{combineGSE}} function.
  Can be used as input to \code{\link{adjustPvalGSE}}.
}

\usage{ data(gseResultsSep) }

\format{
  This object is a \code{list} of length four:
    
  \code{"dat.affy"}{: a \code{list} of legnth two:
    \code{"go"} is a \code{numeric} vector of length 5,
    containing the p-values resulting from gene set enrichment
    analysis of 5 GO terms on Affymetrix gene expression data;
    \code{"keg"} is a \code{numeric} vector of length 5,
    containing the p-values resulting from gene set enrichment
    analysis of 5 KEGG pathways on Affymetrix gene expression data;}
  
  \code{"dat.agilent"}{: a \code{list} of legnth two:
    \code{"go"} is a \code{numeric} vector of length 5,
    containing the p-values resulting from gene set enrichment
    analysis of 5 GO terms on Agilent gene expression data;
    \code{"keg"} is a \code{numeric} vector of length 5,
    containing the p-values resulting from gene set enrichment
    analysis of 5 KEGG pathways on Agilent gene expression data;}
  
  \code{"dat.cnvHarvard"}{: a \code{list} of legnth two:
    \code{"go"} is a \code{numeric} vector of length 5,
    containing the p-values resulting from gene set enrichment
    analysis of 5 GO terms on Harvard CNV data;
    \code{"keg"} is a \code{numeric} vector of length 5,
    containing the p-values resulting from gene set enrichment
    analysis of 5 KEGG pathways on Harvard CNV data;}
  
  \code{"dat.cnvMskcc"}{: a \code{list} of legnth two:
    \code{"go"} is a \code{numeric} vector of length 5,
    containing the p-values resulting from gene set enrichment
    analysis of 5 GO terms on MSKCC CNV data;
    \code{"keg"} is a \code{numeric} vector of length 5,
    containing the p-values resulting from gene set enrichment
    analysis of 5 KEGG pathways on MSKCC CNV data;}

}

\source{
  Computed using the \code{\link{runBatchGSE}} function
  from the TCGA data contained in \code{\link{sepScores}}
  and \code{\link{fgsList}}.
}

\examples{
data(gseResultsSep)
class(gseResultsSep)
names(gseResultsSep)
str(gseResultsSep)
}

\keyword{datasets}
