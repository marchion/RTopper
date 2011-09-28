\name{pheno}
%\Rdversion{1.1}

\alias{pheno}

\docType{data}

\title{
  A test dataset for the RTopper package
}

\description{
  A \code{data.frame} with 2 columns containing the phenotypic class indicator
  for the 95 patients analyzed and used in the examples.
  Can be used as input to \code{\link{convertToDr}}.  
}

\usage{ data(exampleData) }

\format{
  This object is a \code{data.frame}
  with two columns:
  
  \code{"Sample"}{: the first column contains the patients identifiers;}

  \code{"Class"}{: the second columns contain a numeric indicator
    (0 or 1) corresponding to the phenotypic class of each patient;}
  
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
class(pheno)
colnames(pheno)
str(pheno)
}

\keyword{datasets}
