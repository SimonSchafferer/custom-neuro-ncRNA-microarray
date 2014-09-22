\name{plotMAofArrays}
\alias{plotMAofArrays}
\title{This function creates a ggplot object where MA plots of all Arrays are plotted, present in the MA object}
\usage{
plotMAofArrays(MA, subset = NULL, colorVar = "probe.type", ncol = 2,
  smooth = FALSE)
}
\arguments{
  \item{MA}{object of type MA}

  \item{subset}{Integer vector of the column indices
  defining the subset}

  \item{colorVar}{the factor for coloring the dots in the
  plot (could also be NULL) (this vector needs to be
  present in the data.frame MA$genes)}

  \item{ncol}{Number of columns that is used for facetting
  the plot}

  \item{smooth}{Adds a LOESS curve to the plot}
}
\description{
This function creates a ggplot object where MA plots of all
Arrays are plotted, present in the MA object
}

