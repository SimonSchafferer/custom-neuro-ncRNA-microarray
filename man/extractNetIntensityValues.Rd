\name{extractNetIntensityValues}
\alias{extractNetIntensityValues}
\title{Extracts the Net Intensity Values by the identfier: Net intensity (mean) {Cy5} and Net intensity (mean) {Cy3}}
\usage{
extractNetIntensityValues(RG, deleteOther = TRUE)
}
\arguments{
  \item{RG}{Object of type RGList}

  \item{deleteOther}{If true, the values defined under
  RG$other are set to NULL}
}
\value{
RG object with Net Intensity Values as R and G
}
\description{
Extracts the Net Intensity Values by the identfier: Net
intensity (mean) {Cy5} and Net intensity (mean) {Cy3}
}

