\name{computeDrStat}

\alias{computeDrStat}

\title{Computes gene-to-phenotype associations scores}

\description{
  \code{computeDrStat} Computes gene-to-phenotype associations scores,
  using as input the output from \code{convertToDr}.
  }
  
  \usage{
computeDrStat(data, columns = c(1:(ncol(data)-1)), method = "dev", integrate = TRUE)
}

\arguments{
  \item{data}{a list of data.frames containing genomic
    measurements. Each element of \code{dataIntersection}
    must account for the same set of patients(columns) and genes (rows)}
  \item{columns}{a data.frame indicating patients' phenotypic class}
  \item{method}{character, the number of genomic platforms}
  \item{integrate}{logical, wheter to integrate the gene-to-phenotype
    scores across platform or return separates scores for each platform}
}

\details{
  This function allows computing gene-to-phenotype association scores,
  using as input the gene-centered list produced by \code{computeDr}.
  The \code{computeDrStat} function works separately on each
  gene-centered data.frame  created by the \code{convertToDr} function,
  assuming that the phenotype information is stored in the last column
  named \code{"response"}. It is possible computing both separate
  association scores for each platform, as well as an integrated score,
  as specified by the \code{integrate} arguments.
  There are currently three methods available for obtaining the scores
  (see Tyekucheva et al, manuscript under review),
  as specified by the \code{methods} argument:
  
  \code{"dev"}{: this approach computes the score as the difference of deviances;}
  
  \code{"aic"}{: this approach computes the score as the Akaike
    information criterion  for model selection;}
  
  \code{"bic"}{: this approach computes the score as the penalized likelihood ratio;}
}

\value{
  A list of named vectors containing separate or integrared gene-to-phenotype
  association scores.
}

\references{
  Svitlana Tyekucheva, Luigi Marchionni, Rachel Karchin, and Giovanni Parmigiani.
  "Integrating diverse genomic data using gene sets."
  Manuscript submitted.
}

\author{ Luigi Marchionni \email{marchion@jhu.edu} }


%\note{ ~~further notes~~ 
% ~Make other sections like Warning with \section{Warning }{....} ~
%}

%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }

\examples{

###load data
data(exampleData)

###convert
dataDr <- convertToDr(dat, pheno, 4)

###compute the integrated score
bicStatIntegrated <- computeDrStat(dataDr, columns = c(1:4), method="bic", integrate = TRUE)

###compute separate scores for each genomic platform
bicStatSeparate <- computeDrStat(dataDr, columns = c(1:4), method="bic", integrate = FALSE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
