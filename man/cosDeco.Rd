\name{cosDeco}
\alias{cosDeco}

\title{
Create consensus deconvolution.
}
\description{
Function to apply seven deconvolution/signature methods.
}
\usage{
cosDeco <- function(x=df, rnaseq=T, plot=TRUE, ext=FALSE,
                    sig=NULL, anno.1=NULL, anno.2=NULL,
                    cp=NULL, free=FALSE)
}

\arguments{
  \item{x}{
The data frame with samples in column/names and rows/rownames as genes.
}
  \item{rnaseq}{
Boolen indicating if the x is   the platform. Defaults to TRUE.
}
  \item{plot}{
Boolen to plot or not plot the evaluation correlation plots. Only need for development of the tool.
}
  \item{ext}{
Boolen indicating an extension to the in-built signatures, triggered by sig, anno.1, anno.2.
}
  \item{sig}{
Data frame of the extension signature. Needs at least anno.1 for conesus.
}
  \item{anno.1}{
Data frame of the annotation for signature.This is used for the relative conseus.
}
  \item{anno.2}{
Data frame of the annotation for signature.This is used for the absolute conseus.
}
  \item{cp}{
Data frame of the colour to use in the correlation plot. Relevant if plot is true, but not compulsory.
}
  \item{free}{
Boolen indicating if you want to remove methods with restrive licences. Currently affects EPIC.
}
  \item{scale.i}{
  Boolen indicating if you want to individually scale outputs from each methods before combining.
}
  \item{agg.method}{
  String of the aggregation method to use, either "mean" or "geomean".
}
  \item{mini.output}{
Boolen indicating if you want to extract just the final reuslts or everythin.
}
}

\details{
A Consesus approache to cell proportions/deconvulations builds on the strength of aggregation and deconvolution methods. This generates robust estimates of cell proprotions/contents in a gene expression dataset. It provides a flexible paltform to add more signatures.
}

\value{
A DecoCell class object, with input, intermediate results and a list "res.final" of two dataframes and two correlation pdfs in the working directory.
 \item{main_samples}{Table of consesus values for comparing samples }
 \item{main_cells}{Table of consesus values for comparing cells }
  \item{raw_results}{Table of raw outputs of each method}
}

\references{
Depends on:
xcell
MCP
Danaher
Davoli
Rooney
quanTISeq
EPIC
 Uses base R functions.
}

\author{
Chinedu A. Anene, PhD,
Emma Taggart
}

\note{
The packages associated with the methods are configured to install automatically. However, you can installl them directry. "devtools::install_github('dviraran/xCell', force = TRUE)" "devtools::install_github("ebecht/MCPcounter",ref="master", subdir="Source")"  "devtools::install_github("GfellerLab/EPIC", build_vignettes = TRUE)". You need to install devtools.

Installation
From [github](https://github.com/caanene1) with:
devtools::install_github('caanene1/Decosus')
From source in your working directory with:
install.packages("Decosus_0.1.2.tar.gz", repos=NULL, type="source")

Please, read the github page for how to extened the signature and the required columns.
}

\seealso{
The output of this function works well for the deconvolution of any tissue type. For cancer specific consesus, see the methods described by PMID: 31510660  or PMID: 31641033.
}

\keyword{ ~Deconvolution }
\keyword{ ~Constrained least square reqiression }
\keyword{ ~Immune-signature }
\keyword{ ~Descosus }


