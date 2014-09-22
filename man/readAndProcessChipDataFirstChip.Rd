\name{readAndProcessChipDataFirstChip}
\alias{readAndProcessChipDataFirstChip}
\title{Reads and processes the gpr files of two channel arrays. Attention changes commas to dots of the names in the first chip}
\usage{
readAndProcessChipDataFirstChip(currDir, columns = NULL,
  customGridNames = FALSE, specificationFile = "samplesInfo.txt",
  source = "genepix", na.strings = "<NO DATA>")
}
\arguments{
  \item{currDir}{directory were the file is read in}

  \item{columns}{custom column Names}

  \item{customGridNames}{ID are taken from gridName Table
  if set to TRUE}

  \item{specificationFile}{Text file in the directory of
  the gpr files were the gpr file is specified}

  \item{source}{Source of the Laser Scanner}

  \item{na.strings}{The Name of the spots without
  informatio or data (as defined by the scanning software)}
}
\value{
List containing an annotated data frame a target data.frame
and the RG object
}
\description{
Reads and processes the gpr files of two channel arrays.
Attention changes commas to dots of the names in the first
chip
}

