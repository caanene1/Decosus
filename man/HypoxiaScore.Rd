\name{HypoxiaScore}
\alias{HypoxiaScore}

\title{
Creates consesus hypoxia score from published gene expression signatures.
}
\description{
Function to consolidate seven hypoxia signatures.
}
\usage{
HypoxiaScore(x = df)

}

\arguments{
  \item{x}{
The data frame with samples in column/names and rows as genes.
}

}

\details{
A Consesus approache to hypoxia score using seven published gene expression signatures.
This generates robust estimate of hypoxia.

}

\value{
Columns of hypoxia score for each signature and column HS for consesus
 \item{result}{Table of hypoxia score and consesus values }
}

\references{
Depends on:
Buffa et al 2010
Winter et al 2007
Ragnum et al., 2015
Eustace et al., 2013
Sorensen et al., 2010
Elvidge et al., 2006
Hu et al., 2009
Uses base R functions.
}

\author{
Chinedu A. Anene, PhD
}





