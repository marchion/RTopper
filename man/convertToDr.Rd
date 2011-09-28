\name{convertToDr}

\alias{convertToDr}

\title{Converts genomic data to a list suitable for computing
  gene-to-phenotype scores}

\description{
  \code{convertToDr} converts genomic data into a list further
  used for computing gene-to-phenotype association scores.}
  
\usage{
convertToDr(dataIntersection, response, nPlatforms = length(data))
}

\arguments{
  \item{dataIntersection}{a list of data.frames containing genomic
    measurements. Each element of \code{dataIntersection}
    must account for the same set of patients(columns) and genes (rows)}
  \item{response}{a data.frame indicating patients' phenotypic class}
  \item{nPlatforms}{numeric, the number of genomic platforms}
}

\details{
  This function coverts a list of data.frames containing distinct
  genomic measurtements performed on the same patients
  into a gene-centered used in further analyses for computing
  gene-to-phenotype scores.
  Data.frame in the input list (\code{dataIntersection}) must
  have the same dimentions, with columns being patients,
  and rows being genes. Column names identify the patients,
  while rownames identify the genes.
  The argument \code{response} is used to pass phenotypic information
  about samples to be analyzed. This is a simple two columns
  data.frame in which the first column correspond to patients
  identifiers, and the second column to the phenotypic response
  encoded as binary class (using the integers 0 and 1).
  The \code{nPlatforms} argument specifies the number of platforms
  that will be analyzed.
}

\value{
  A list of data.frames, one for each analyzed gene, summarizing
  all genomic measurements and phenotipic information across
  patients and platforms.
}

%\note{
%BLA
%}

\references{
  Svitlana Tyekucheva, Luigi Marchionni, Rachel Karchin,and Giovanni Parmigiani
  "Integrating diverse genomic data using gene sets."
  Manuscript submitted.
}

\author{ Luigi Marchionni \email{marchion@jhu.edu} }

%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }

\examples{

###load data
data(exampleData)

###convert
dataDr <- convertToDr(dat, pheno, 4)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
