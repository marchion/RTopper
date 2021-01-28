\name{runBatchGSE}
\alias{runBatchGSE}


\title{ To perform GSE analysis over multipple experiments and functional themes}

\description{
  The  \code{runBatchGSE}  function enables performing Gene Set
  Enrichment analysis over multiple ranking statistics and multiple
  lists of gene sets.
  By default this function is an interface to the \code{\link{geneSetTest}}
  in the \code{limma} package, and most of the arguments passed to
  \code{runBatchGSE} are indeed passed to such lower
  level function.
  As an alternative the user can also define and pass to
  \code{runBatchGSE} a custom function, defining the ranking statistics
  and the gene set membership in the same way done
  for \code{\link{geneSetTest}} (see Details below).
}

\usage{
runBatchGSE(dataList, fgsList, ...)
}

\arguments{
   \item{dataList}{a list containing the gene-to-phenotype scores to be used
   as ranking statistics in the GSE analysis. This list is usually
   produced by running \code{computeDrStat}}
   \item{fgsList}{a list of FGS collection, in which each element is a list of character vectors,
   one for each gene set}
  \item{\dots}{additional arguments to be passed to lower level functions (see details below)}  
}


\details{
  This function performs enrichment analysis for all the gene-to-phenotype
  scores (argument \code{dataList}) passed to it over a list of F
  unctional Gene Set (FGS) (argument \code{fgsList}), returning
  a p-value for each FGS.
  Additional arguments can be bassed to this function to modify
  the way the enrichment test is performed, as follows:

  \code{absolute}{  logical, this specifies whether the absolute values of
    the ranking statistics  should be used in the test (the default
    being TRUE)}
  
  \code{gseFunc}{  a function to perform GSE analysis. If not specified
    the default is the \code{\link{geneSetTest}} function from the
    \code{limma} package. If a function is specified by the user, the
    membership of the analyzed genes to a FGS, and the ranking
    statistics must be defined in the same way this is done for
    \code{\link{geneSetTest}}, and the new function must
    return an integer (usually a p-value) (see the help for
    \code{\link{geneSetTest}})}

  The following main arguments are used by \code{\link{geneSetTest}}:
  
  \code{type}{  character, specifies the type of statistics used to rank
    the genes by \code{geneSetTest}: \code{'f'} for F-like statistics
    (default), \code{'t'} for t-like statistics, or \code{'auto'} for an
    educated guess}
  
  \code{alternative}{  character, defines the alternative with the
    following possible options: \code{'mixed'} (default),
    \code{'either'}, \code{'up'} or \code{'down'},
    \code{'two.sided'}, \code{'greater'}, or \code{'less'}}
  
  \code{ranks.only}{  logical, if \code{TRUE} (default) only ranks will be
    used by \code{geneSetTest}}
  
  \code{nsim}{  numeric, the number of randomly selected sets of genes to
    be used  in simulations to  compute the p-value}
}



\value{
  The output is a list of lists containing the set of enrichment results
  for all gene-to-phenotype scores and FGS collections used as input.
}

\references{
  Svitlana Tyekucheva, Luigi Marchionni, Rachel Karchin, and Giovanni Parmigiani.
  "Integrating diverse genomic data using gene sets."
  Manuscript submitted.
}

\author{ Luigi Marchionni \email{marchion@jhu.edu} }


\examples{

###require limma to run the example
require(limma)

###load integrated gene-to-phenotype scores
data(intScores)

###load separate gene-to-phenotype scores
data(sepScores)

###load list of functional gene sets
data(fgsList)

###run GSE analysis in batch with default parameters
gseABS.int <- runBatchGSE(dataList=intScores, fgsList=fgsList)

###run GSE analysis in batch with alternative parameters
gseABS.sep <- runBatchGSE(dataList=sepScores, fgsList=fgsList, 
        absolute=FALSE, type="t", alternative="up")

###run GSE analysis in batch passing an enrichment function
gseUP.int.2 <- runBatchGSE(dataList=intScores, fgsList=fgsList,
			   absolute=FALSE, gseFunc=wilcoxGST, alternative="up")

###define and use a new enrichment function
gseFunc <- function (selected, statistics, threshold) {
	diffExpGenes <- statistics > threshold
	tab <- table(diffExpGenes, selected)
	pVal <- fisher.test(tab)[["p.value"]]
	}
gseUP.sep.2 <- runBatchGSE(dataList=sepScores, fgsList=fgsList,
				 absolute=FALSE, gseFunc=gseFunc, threshold=7.5)

}

%\note{ ~~further notes~~ 
% ~Make other sections like Warning with \section{Warning }{....} ~
%}

\keyword{ manip }
