\name{modifyAnalysis}
\alias{modifyAnalysis}
\title{Modifies the analysis Template}
\usage{
modifyAnalysis(analysisTemplate, gprDir = gprDir,
  analysisPath = analysisPath, normMethod = normMethod, ndups = ndups,
  spacing = spacing, minDetectLim = minDetectLim,
  subsetByManualControl = subsetByManualControl,
  subsetByPMvsMM = subsetByPMvsMM, orderingVar = orderingVar,
  addInterestCand = addInterestCand)
}
\arguments{
  \item{analysisTemplate}{Analysis template (chip version)}

  \item{gprDir}{directory of gpr file}

  \item{analysisPath}{path were the html and quality
  control files should be written}

  \item{normMethod}{normalization method one wants to use
  (quantile or vsn)}

  \item{ndups}{number of duplicated spots on the chip}

  \item{spacing}{spacing between the duplicated spots}

  \item{minDetectLim}{minimum detection limit obtained from
  self-self hybridization experiment}

  \item{subsetByManualControl}{if data should be subsetted
  by spike in control (Chip 1 -> no spike in controls)}

  \item{subsetByPMvsMM}{if data should be subsetted by PM >
  MM condition}

  \item{orderingVar}{By which variable should spot ordering
  occur}

  \item{addInterestCand}{Candidate IDs that are added as
  labels in the plots and in the table}
}
\description{
Modifies the analysis Template
}

