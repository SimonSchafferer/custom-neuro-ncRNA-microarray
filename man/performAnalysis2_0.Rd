\name{performAnalysis2_0}
\alias{performAnalysis2_0}
\title{Perform an analysis on the Second 2.0 chip}
\usage{
performAnalysis2_0(gprDir = file.path(getwd(), "../"),
  analysisPath = file.path(getwd()), normMethod = "quantile", ndups = 4,
  spacing = 1, minDetectLim = 0.15, subsetByManualControl = TRUE,
  subsetByPMvsMM = TRUE, orderingVar = "Name",
  addInterestCand = c("900000", "20285", "20285_1", "20285_2"),
  executeKnitR = TRUE)
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
Perform an analysis on the Second 2.0 chip
}

