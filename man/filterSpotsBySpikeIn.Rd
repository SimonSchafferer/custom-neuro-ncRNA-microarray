\name{filterSpotsBySpikeIn}
\alias{filterSpotsBySpikeIn}
\title{Filtering of spot values by the condition PM > MM and spot intensity > spike in spot intensity}
\usage{
filterSpotsBySpikeIn(MA, normMethod = "", subsetByManualControl = FALSE,
  subsetByPMvsMM = FALSE, ndups = 4, spacing = 1, weights = NULL,
  design, writeLog = TRUE)
}
\arguments{
  \item{MA}{MAList object as defined by limma}

  \item{normMethod}{Optional normalization method name that
  attached to the log file name (e.g. _vsn_)}

  \item{subsetByManualControl}{Boolean, data may be
  subsetted by the spike in controls}

  \item{subsetByPMvsMMBoolean,}{data may be subsetted by
  restriction PM > MM}

  \item{ndups}{Number of duplicated spots}

  \item{spacing}{Placement of the duplicated spots (in
  default mode it is assumed that spots are sorted)}

  \item{weights}{Parameter for avedups function}

  \item{design}{Design of the experiment is considered when
  caluclating restrictions
  (subsetByManualControl,subsetByPMvsMM)}

  \item{writeLog}{Writes a log file of the filtering steps}
}
\value{
Filtered MA object
}
\description{
Filtering of spot values by the condition PM > MM and spot
intensity > spike in spot intensity
}

