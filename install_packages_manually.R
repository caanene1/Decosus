
# List of packages grouped by source
cran_packages <- c("MASS", "readxl", "data.table", "corrplot", "devtools")
bioc_packages <- c("GSEABase", "GSVA",  "preprocessCore")

  # Set CRAN mirror if not already set
if (is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
}

# Install missing CRAN packages
missing_cran <- setdiff(cran_packages, installed.packages()[, "Package"])
if (length(missing_cran) > 0) {
    install.packages(missing_cran)
}

# Install missing Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

missing_bioc <- setdiff(bioc_packages, installed.packages()[, "Package"])
if (length(missing_bioc) > 0) {
    BiocManager::install(missing_bioc)
}


## Install github source
devtools::install_github('dviraran/xCell', force = TRUE)
devtools::install_github("ebecht/MCPcounter",ref="master", subdir="Source")
devtools::install_github("GfellerLab/EPIC", build_vignettes = TRUE)
devtools::install_github("cran/limSolve")
