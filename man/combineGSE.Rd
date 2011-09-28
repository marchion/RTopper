\name{combineGSE}

\alias{combineGSE}


\title{Combines GSE results}

\description{
  Combines GSE results obtained from a sepatate set of
  gene-to-phenotypes scores
}

\usage{
combineGSE(gseOut, method)
}

\arguments{
  \item{gseOut}{a list of lists containing the enrichment results
    to be combined. This is usually the the output of
    \code{\link{runBatchGSE}} obtained from a set of distinct
    genes-to-phenotype scores (usually one per
    genomic platform)}
  \item{method}{character, this argument specifies the method
    used to combine the enrichment results obtained from
    distinct genes-to-phenotype scores (usually one per
    genomic platform).
    Available options are the computation of the geometric
    or arithmetic  means, the use of the median, the selection of
    the minimun or the maximum enrichment score, and
    the random selection of a score (respectively
    \code{"geometricMean"}, \code{"mean"}, \code{"median"},
    \code{"min"}, \code{"max"}, and \code{"random"})}
}

\details{
  This function summarize enrichment results obtained from
  distinct gene-to-phenotypes scores  (usually one per
  genomic platform) by one of several alternative methods.
}

\value{
  The output is a list of lists containing integrated enrichment results
  for all  FGS collections
}

\references{
  Svitlana Tyekucheva, Luigi Marchionni, Rachel Karchin, and Giovanni Parmigiani.
  "Integrating diverse genomic data using gene sets."
  Manuscript submitted.
}

\author{ Luigi Marchionni \email{marchion@jhu.edu} }


\examples{

###load gse analysis results for separate gene-to-phenotype score
data(gseResultsSep)

###combine enrichment score results using geometric mean
gseABS.sep.geoMean <- combineGSE(gseResultsSep, method="geometricMean")

###combine enrichment score results using maximum value
gseABS.sep.max <- combineGSE(gseResultsSep, method="max")

}

%\note{ ~~further notes~~ 
% ~Make other sections like Warning with \section{Warning }{....} ~
%}

%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }

\keyword{ manip }

