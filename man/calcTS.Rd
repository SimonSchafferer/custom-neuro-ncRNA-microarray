\name{calcTS}
\alias{calcTS}
\title{Calculates the optimal table size for displaying a data.frame as googleVis table in the html file, when using markdown}
\usage{
calcTS(tDF, maxRowsToDisplay = 20, charMultiplicationFactorInPixels = 8,
  tableHeightInPixels = 32)
}
\arguments{
  \item{data.frame}{that should represent the table to
  display}

  \item{maxRowsToDisplay}{Maximum number of rows that
  should be displayed before paging occurs (default 20)}

  \item{charMultiplicationFactorInPixels}{Size of a char in
  the browser in pixels (default 8)}

  \item{tableHeightInPixels}{Height of the table rows in
  pixels (default 32)}
}
\value{
GvisTable object that can be called with plot() in a code
chunk
}
\description{
Calculates the optimal table size for displaying a
data.frame as googleVis table in the html file, when using
markdown
}

