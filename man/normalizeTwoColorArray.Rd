\name{normalizeTwoColorArray}
\alias{normalizeTwoColorArray}
\title{Normalizes a RG object by a given method (eather vsn or quantile normalization)
When choosing quantile normalization an offset of 0.1 for zero values is defined (methods used: printtiploess for within and Aquantile for between arrays))}
\usage{
normalizeTwoColorArray(RG, method = "quantile")
}
\arguments{
  \item{RG}{Object of type RGList}

  \item{method}{Choose between 'quantile' or 'vsn'}
}
\value{
Normalized RG object
}
\description{
Normalizes a RG object by a given method (eather vsn or
quantile normalization) When choosing quantile
normalization an offset of 0.1 for zero values is defined
(methods used: printtiploess for within and Aquantile for
between arrays))
}

