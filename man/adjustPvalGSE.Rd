\name{adjustPvalGSE}

\alias{adjustPvalGSE}

\title{Multiple testing correction}

\description{
  This function is an interface to the \code{\link{mt.rawp2adjp}} function
  contained in the \code{multtest} package.
  \code{adjustPvalGSE} works on outputs from \code{\link{runBatchGSE}}
  and \code{\link{combineGSE}} returning adjusted p-values
}

\usage{
adjustPvalGSE(gseOut, proc = "BH", alpha = 0.05 , na.rm = FALSE)
}

\arguments{
  \item{gseOut}{list of lists, either the output from \code{\link{runBatchGSE}},
    or the output from \code{\link{combineGSE}}}
  \item{proc}{character, the method to be used for p-values adjusting.
    This parameter will be passed to the \code{mt.rawp2adjp}
    function from the \code{multtest} package. The available
    options include: \code{"Bonferroni"}, \code{"Holm"}, \code{"Hochberg"},
      \code{"SidakSS"}, code{"SidakSD"}, \code{"BH"} (the default),
      \code{"BY"}, \code{"ABH"}, and \code{"TSBH"}}
  \item{alpha}{numeric, the nominal type I error rate}
  \item{na.rm}{logical, the option for handling NA values in the
    list of raw p-values}
}

\details{
  The \code{adjustPvalGSE} function performs p-value adjusting
  for multiple testing correction on the list of lists resulting from
  enrichment analysis obtained using the \code{\link{runBatchGSE}}
  and \code{\link{combineGSE}} functions.
  This functions is based on the \code{\link{mt.rawp2adjp}} function
  contained in the \code{multtest} package.
}

\value{
  For each vector of p-value contained in the \code{gseOut} input object
  a data.frame is returned, containing original p-value and corrected
  p-values
}


\references{
  Svitlana Tyekucheva, Luigi Marchionni, Rachel Karchin, and Giovanni Parmigiani.
  "Integrating diverse genomic data using gene sets."
  Manuscript submitted.
}

\author{ Luigi Marchionni \email{marchion@jh.edu}}

%\note{ ~~further notes~~ 
% ~Make other sections like Warning with \section{Warning }{....} ~
%}

%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }

\examples{

###load gse analysis results for separate gene-to-phenotype score
data(gseResultsSep)

###adjust for multiple testing using the Benjamini and Hochberg method
gseABS.int.BH <- adjustPvalGSE(gseResultsSep)

###adjust for multiple testing using the Holm method
gseABS.int.holm <- adjustPvalGSE(gseResultsSep, proc = "Holm")

}

\keyword{manip}

