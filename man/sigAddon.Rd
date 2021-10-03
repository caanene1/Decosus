\name{sigAddon}
\alias{sigAddon}

\title{
Process new signatures for inclusion into the decosus tool
}
\description{
Function to process signatures from the CellMarker database.
}
\usage{
sigAddon(path=url.inl, is.web = T, min.set=3)
}

\arguments{
  \item{path}{
The path to the CellMarker database. The path can be a web link or local.
The function expects this file to contain minimally the "cellName" and "geneSymbol" column.
}
  \item{is.web}{
Boolen indicating if the path is a web url. Defaults to TRUE.
}
  \item{min.set}{
Integer indicating min number of genes in the retained signatures.
}

}

\details{
This function provides a method to integrate the CellMarker signature to the analysis.
}

\value{
A list of data dataframes.
 \item{1}{Table of the processed signature}
 \item{2}{Raw table used, either from web or local}
}

\references{
Depends on:
 base R functions.
}

\author{
Chinedu A. Anene, PhD
}
\note{
This function expects your data to come from http://bio-bigdata.hrbmu.edu.cn/CellMarker/download.jsp
Example url is "http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt".

}

\seealso{
The output of this function works well for the Decosus tool.
}

\keyword{ ~Preprocessing }
