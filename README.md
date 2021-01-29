# Decosus
Decosus consolidates the results of 7 independently published deconvolution software to generate a robust estimation of cell composition in heterogeneous tissue from bulk expression data.

# Description
The inability to identify the phenotype and abundance of cells is a major limitation of bulk expression data sets. Deconvolution methods use gene signatures from purified cell populations to estimate the presence of different cells. Decosus combines the results from 7 independent deconvolution methods in order to more accurately estimate the abundance of key immune and stromal cells than by using one method alone. A unique advantage of Decosus is the ability to adapt the tool to facilitate the type of downstream analyses required, thus two versions of Decosus are possible: 
'''
 - Sample: includes data from all 7 methods to give the most comprehensive overview of cell composition, and allows across-sample comparison for each cell type.
 - Cell: only uses methods which enable estimation of cell proportion, in order to  compare the abundance of one cell type to another, as well as across-sample comparison. 
Further, the tool can use different cell mapping depending on user needs.
'''


# Dependencies
Decosus depends on certain published signatures and methods. The packages associated with the methods are configured to install automatically. However, you can install them directly if you have problems getting them on your machine.
Use: 
  devtools::install_github('dviraran/xCell', force = TRUE)
  devtools::install_github("ebecht/MCPcounter",ref="master", subdir="Source")
  devtools::install_github("GfellerLab/EPIC", build_vignettes = TRUE)
You may need to install devtools.

# Installation
You can install the released version of Decosus from [github](https://github.com/caanene1) with:

``` r
devtools::install_github('caanene1/Decosus')
```
or 
``` r
devtools::install_github('BioInforCore-BCI/Decosus')
```

# Example
To use Decosus on an expression matrix, simply run: 
```{r example}
results <- Decosus::cosDeco(x = df, platform = "Array", map="Normal", plot.corr=FALSE)
```

# Additional usage note
Run array or bulk RNAseq data as a data frame through the tool to generate 3 outputs:
 - Main_samples for across-sample comparison version of the tool
 - Main_cells for the within-sample comparison version of the tool
 - Raw_results for the results of each individual method in one table

Each of the contributing methods utilises a different combination of cells; when a particular cell type appears in 2 or more of the contributing methods, a consensus for that cell type is generated and this cell is labelled cell_name_consensus in the output tables.

When a cell type is unique to one of the contributing methods, it is labelled cell_name_method_name.

All consensus cells and all unique cells for each version are included in the results table, but if needed, irrelevant or non-consensus cells can be removed from the results once generated.

As the Cells version of the tool uses fewer contributing methods, the main_cells output table contains fewer cell types than main_samples. 
