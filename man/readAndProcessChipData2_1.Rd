\name{readAndProcessChipData2_1}
\alias{readAndProcessChipData2_1}
\title{Reads and processes the gpr files of two channel arrays For The NRC5 gal file or Chip 2.1.}
\usage{
readAndProcessChipData2_1(currDir, specificationFile = "samplesInfo.txt",
  source = "genepix", na.strings = "<NO DATA>", columns, otherColumns,
  GALFile, ...)
}
\arguments{
  \item{currDir}{directory were the file is read in}

  \item{specificationFile}{Text file in the directory of
  the gpr files were the gpr file is specified}

  \item{source}{Source of the Laser Scanner}

  \item{na.strings}{The Name of the spots without
  informatio or data (as defined by the scanning software)}

  \item{columns}{custom column Names}

  \item{otherColumns}{Additional columns}

  \item{GALFile}{the GAL file of the chip}
}
\value{
List containing an annotated data frame a target data.frame
and the RG object
}
\description{
Reads and processes the gpr files of two channel arrays For
The NRC5 gal file or Chip 2.1.
}

