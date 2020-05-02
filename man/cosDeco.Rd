\name{cosDeco}
\alias{cosDeco}

\title{
Creates consesus from published cell deconvolution/aggregation methods.
}
\description{
Function to consolidate seven deconvolution method and signatures.
}
\usage{
cosDeco(x = df, platform = "Array")
}

\arguments{
  \item{x}{
The data frame with samples in column/names and rows/rownames as genes.
}
  \item{platform}{
String value indicating the platform. Defaults to Array.
}

}
\details{
A Consesus approache to cell proportions/deconvulations builds on the strength of aggregation and deconvolution methods. This generates robust estimates of cell proprotions/contents in a gene expression dataset.

}
\value{
A list of three dataframes and two correlation pdfs in the working directory.
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
 Uses base R functions and ??.
}

\author{
Chinedu A. Anene, PhD
}
\note{
The packages associated with the methods are configured to install automatically. However, you can installl them directry. "devtools::install_github('dviraran/xCell', force = TRUE)" "devtools::install_github("ebecht/MCPcounter",ref="master", subdir="Source")"  "devtools::install_github("GfellerLab/EPIC", build_vignettes = TRUE)". You need to install devtools.

Installation
From [github](https://github.com/caanene1) with:
devtools::install_github('caanene1/Decosus')
From source in your working directory with:
install.packages("Decosus_0.1.2.tar.gz", repos=NULL, type="source")
}

\seealso{
The output of this function works well for the deconvolution of any tissue type. For cancer specific consesus, see the methods described by PMID: 31510660  or PMID: 31641033.
}

\keyword{ ~Deconvolution }
\keyword{ ~Constrained least square reqiression }
\keyword{ ~Immune-signature }
\keyword{ ~Descosus }


