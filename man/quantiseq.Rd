\name{quantiseq}
\alias{quantiseq}

\title{
Use quanTIseq to deconvolute a gene expression matrix.
}
\description{
Function to apply the qunatiseq deconvolution method as described by bioRxiv 223180. https://doi.org/10.1101/223180.
}
\usage{
quantiseq(mix.mat, arrays = FALSE, mRNAscale = TRUE, method = "lsei")
}

\arguments{
  \item{mix.mat}{
The expression matrix to deconvolute, with samples in column and genes in rows.
}

 \item{arrays}{
The platform used to generate the expression dataset.
}

 \item{mRNAscale}{
Option to scale expression matrix, defults to TRUE.
}

 \item{method}{
Option to select the right method, defults to "lsei".
}

}
\details{
Immune devconvolution using the qunatiseq mnethod described by bioRxiv 223180. https://doi.org/10.1101/223180.

}
\value{
Cell proportions.
 \item{t_results}{Table of deconvolution across samples.}
}
\references{
 https://doi.org/10.1101/223180.
 https://github.com/FFinotello/quanTIseq

}
\author{
Chinedu A. Anene, PhD
}
\note{
The user need to ensure that colnames and rownames are set with samples and genes, respectively.
}

\seealso{
The output of this function is incorporated directly to the Decosus package for consesus.
}

\keyword{ ~Deconvolution }
\keyword{ ~Constrained least square reqiression }
\keyword{ ~Immune-signature }
\keyword{ ~Descosus }
