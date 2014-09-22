\name{addProbeTypeAndChipTypeToRGGenes}
\alias{addProbeTypeAndChipTypeToRGGenes}
\title{Adds the probe and chipType to the RG$genes table based on the database connection given}
\usage{
addProbeTypeAndChipTypeToRGGenes(RG, dbcon, oligoTable,
  anchorName = "unknown", noSignal = "noSig", refCol = "Name")
}
\arguments{
  \item{RG}{Object of type RGList}

  \item{dbcon}{Database connection}

  \item{oligoTable}{If no database connection is provided
  the oligoTable has to be provided containing the columns
  category,chipType,oligoTableID}

  \item{anchorName}{Name of the Anchor spots}

  \item{noSignal}{Name of the spots lacking a signal}
}
\value{
The 'genes' table of the RG object provided
}
\description{
Adds the probe and chipType to the RG$genes table based on
the database connection given
}

