library(BioML)
?deLimma

comps <- list("")

Res <- deLimma(count = GSE76801, samp = samp, minCpm = 1, minSample = 3, 
               type = "RNA-Seq")
installed.packages()
BiocManager::install("edgeR")

update.packages(checkBuilt=TRUE, ask=FALSE)
Sys.getenv("PATH")
Sys.which("stats.dll")
