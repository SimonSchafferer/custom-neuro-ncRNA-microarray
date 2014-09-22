\name{calcMeanValueOfSpots}
\alias{calcMeanValueOfSpots}
\title{Calculates the mean or median value of the spots on a chip by a given Identifier}
\usage{
calcMeanValueOfSpots(x, byType = "Name", method = "mean")
}
\arguments{
  \item{x}{An object of type EListRaw, RGList or EList}

  \item{byType}{An identifier that is used for grouping the
  values e.g. ID}

  \item{method}{The method name as character for spot
  averaging (currently 'mean' and 'median' are supported)
  this will then call colMeans or colMedians}
}
\value{
The averaged object
}
\description{
Calculates the mean or median value of the spots on a chip
by a given Identifier
}

