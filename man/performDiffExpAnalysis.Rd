\name{performDiffExpAnalysis}
\alias{performDiffExpAnalysis}
\title{Performs a differential expression analysis via limma}
\usage{
performDiffExpAnalysis(MA, aveDups = FALSE, design = NULL, block = NULL,
  blockpairs = 2, ndups = 4, spacing = 1)
}
\arguments{
  \item{MA}{An object of type MAList}

  \item{aveDups}{The number of duplicated spots in the MA
  object}

  \item{design}{Design of the experiment (if null, a dye
  swap design is considered (1,-1,1,-1,...))}

  \item{block}{Blocking vector}

  \item{blockpairs}{For creation of a default blocking
  vector (default creation would be 1,1,2,2,...)}

  \item{ndups}{Number of duplicated spots}

  \item{spacing}{If spacing is not defined a sorted MA
  object is considered}
}
\value{
List containing a resulting Table from eBayes, the model
fit and the fit after eBayes
}
\description{
Performs a differential expression analysis via limma
}

