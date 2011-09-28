\name{dat}
%\Rdversion{1.1}

\alias{dat}

\docType{data}

\title{
  A test dataset for the RTopper package
}

\description{
  A small subset of pre-processed Glioblasoma Multiforme (GBM)
  genomic  data from The Cancer Genome Atlas (TCGA) project,
  encompassing Differential Gene Expression (DGE), and 
  Copy Number Variation (CNV).
  Can be used as input to \code{\link{convertToDr}}.  
}

\usage{ data(exampleData) }

\format{
  This object is a \code{list} of 4 \code{data.frames}
  containing genomic measurements
  obtained across distinct genomic scopes
  (copy number variation and gene expression),
  platforms (Affymetrix and Agiles),  and laboratories.
  In particular each \code{data.frame} consist of 500 gene measurements (by rows),
  for 95 distinct patients (by columns) from the following  4 distinct platforms:

  \code{"dat.affy"}{: DGE obtained using Affymetrix microarrays;}
  
  \code{"dat.agilent"}{: DGE obtained using Agilent microarrays;}
  
  \code{"dat.cnvHarvard"}{: CNV data obtained at Harvard;}
  
  \code{"dat.cnvMskcc"}{: CNV data obtained at Memorial Sloan Ketterng Cancer Center;}

}

\source{
   The Cancer Genome Atlas (TCGA) project
  \url{http://cancergenome.nih.gov/}
}

\references{
  The Cancer Genome Atlas (TCGA) Research Network.
  \dQuote{Comprehensive genomic characterization
    defines human glioblastoma genes and core pathways}.
  \emph{Nature}, 2008, October 23; 455(7216): 1061-1068
  
}

\examples{
data(exampleData)
ls()
class(dat)
names(dat)
sapply(dat,class)
sapply(dat,dim)
}

\keyword{datasets}
