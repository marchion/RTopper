\name{fgsList}
%\Rdversion{1.1}

\alias{fgsList}

\docType{data}

\title{
  A list of Functional Gene Set (FGS) to be used to run the examples in the RTopper package
}

\description{
  A \code{list} containing distinct types of FGS (i.e. Gene Ontology,
  KEGG pathways). Each FGS is type is a \code{list}
  of named \code{chararacter} vectors, one for each FGS,
  containing the gene identifiers.  Vectors names describe the FGS.
  Can be used as input to \code{\link{runBatchGSE}}.  
}

\usage{ data(fgsList) }

\format{
  This object is a \code{list} of length two:
  
  \code{"go"}{: this is a \code{list} of 5 \code{character} vectors,
    corresponding to 5 distinct Gene Ontology (GO) terms. Genes
    annotated to each GO term are identified by their gene symbol;}
  
  \code{"kegg"}{: this is a \code{list} of 5 \code{character} vectors,
    corresponding to 5 distinct KEGG pathways. Genes
    annotated to each KEGG pathway are identified by their gene symbol;}
  
}

\source{
  The FGS were obtained from the \code{org.Hs.eg.db} pakage,
  (use the \code{\link[org.Hs.eg.db]{org.Hs.eg}} function
  to see the content);
  These FGS were annotated using data from
  \code{GO.db} and \code{KEGG.db} packages
  (use the \code{\link[GO.db]{GO}} and
  \code{\link[KEGG.db]{KEGG}} functions to see the content).
}

\examples{
data(fgsList)
class(fgsList)
names(fgsList)
str(fgsList)
}

\keyword{datasets}
