\name{RTopper-package}

\alias{RTopper-package}

\alias{RTopper}

\docType{package}

\title{
  A package to perform run Gene Set Enrichment across genomic platforms
}

\description{
  Gene sets analysis considers whether genes sharing a biological
  property also behave in a related way in experimental data.
  This technique is commonly used in high throughput genomic
  analyses to assist results interpretation, and has been successfully
  applied in cancer genome projects for integrating information from
  multiple genome-wide assays. The \code{RTopper} package uses
  gene sets analysis to overcome the diversity of genomic data
  providing the statistical framework for integration across data types.
}

\section{RTopper package features}{
  \code{Rtopper} enables two gene set-based data integration
  approaches:
  \describe{
    \item{1}{ Integration+GSA: computing integrated gene-to-phenotype
      association scores, followed by conventional gene sets analysis;}
    \item{2}{ GSA+Integration: computing consensus significance score
      after all data types are analyzed individually;}
    \item{3}{ Use of alternative enrichment test: \code{RTopper} uses
      the Wilcoxon rank-sum test for enrichment testing, however
      alternative tests can be defined and  used;}
    \item{4}{ Multiple testing correction: \code{RTopper} enables
      adjustment of p-values obtained from enrichment analysis;}
  }
}

\details{
  \tabular{ll}{
    Package: \tab RTopper\cr
    Type: \tab Package\cr
    Version: \tab 0.1\cr
    Date: \tab 2010-11-12\cr
    License: \tab GPL version 3 or newer\cr
  }
}

\author{
Luigi Marchionni \email{marchion@jhu.edu}
}

\references{
  Svitlana Tyekucheva, Luigi Marchionni, Rachel Karchin, and Giovanni Parmigiani.
  "Integrating diverse genomic data using gene sets."
  Manuscript submitted.
}


\keyword{ package }

