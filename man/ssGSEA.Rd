\name{sigAddon}
\alias{sigAddon}

\title{
Process new signatures for inclusion into the decosus tool
}
\description{
Function to process signatures from the CellMarker database.

Understanding the Alpha Parameter
	•	Alpha (α) controls how much weight is given to the ranking of genes within the gene set:
	•	α = 0: All genes in the gene set are weighted equally (no bias based on expression rank).
	•	α > 0: Higher-ranked genes receive more weight. Typically, α = 0.25 is used as a default.
	•	Larger α values further emphasize highly ranked genes and may reduce the influence of lower-ranked genes in the gene set.

Factors for Selecting Alpha
	1.	Data Distribution:
	•	If your expression data has a wide dynamic range (e.g., RNA-Seq data), a higher alpha (0.5 or 1.0) may better capture the most highly expressed genes.
	•	For microarray data with a narrower range, a lower alpha (0.25) may work better.
	2.	Gene Set Size:
	•	Smaller gene sets may benefit from higher alpha values to emphasize a few critical genes.
	•	Larger gene sets might perform better with a lower alpha to balance the contributions across more genes.
	3.	Biological Context:
	•	If the biological process you’re studying is influenced by highly expressed or key driver genes, consider a higher alpha.
	•	If the process is governed by the collective effect of multiple genes, a lower alpha may be more appropriate.
	4.	Sensitivity Analysis:
	•	Test multiple alpha values (e.g., 0, 0.25, 0.5, 1.0) to observe how enrichment scores change.
	•	Assess the stability of your results across different values and select the alpha that provides the most consistent and biologically meaningful results.


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
