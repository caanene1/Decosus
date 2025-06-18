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

Or using the install_packages_manually.R file in the root of the repo to achieve this. 
NOTE: The system needs to have GNU Fortran compiler to setup limSolva. 
In most cases you will have this already, Macs maybe any issue so head over to https://mac.r-project.org/tools/.
limSolva is currently archived by CRAN so it is installing directly from that archive. I will attempt to remove the limSolve dependecy from the quantiseq implementationif possible.

# Installation
Decosus is tested for the current version of R 4.5.1
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
results <- Decosus::cosDeco(x=df, rnaseq=T, plot=TRUE, ext=FALSE,
                    sig=NULL, anno.1=NULL, anno.2=NULL,
                    cp=NULL, free=FALSE)

# Use below to get more detials on the arguments
?cosDeco
```
In the example code above, df is the gene expression profiles in your samples and it must have just one column of gene ID or symbols as shown below.

| gene  | sample1 | sample4 | sample3 |
| --- | --- | --- | --- |
| TGF | 30 | 20 | 1 | 
| PDL1  |  100  | 11  |  300  |
| CD21  |  10  | 6  |  5 |


# Additional usage note
Run array or bulk RNAseq data as a data frame through the tool to generate 3 outputs:
 - Main_samples for across-sample comparison version of the tool
 - Main_cells for the within-sample comparison version of the tool
 - Raw_results for the results of each individual method in one table

Each of the contributing methods utilises a different combination of cells; when a particular cell type appears in 2 or more of the contributing methods, a consensus for that cell type is generated and this cell is labelled cell_name_consensus in the output tables.

When a cell type is unique to one of the contributing methods, it is labelled cell_name_method_name.

All consensus cells and all unique cells for each version are included in the results table, but if needed, irrelevant or non-consensus cells can be removed from the results once generated.

As the Cells version of the tool uses fewer contributing methods, the main_cells output table contains fewer cell types than main_samples. 

# Extension of signatures
To extend the signatures as described in the manuscript, set the `ext` argument to TRUE, and provide the signature and annotation using the `sig`, `anno.1`, and `anno.2` arguments.

A template file for creating these extended signatures is available in the /template_for_adding_extension folder. The most important requirement is to ensure that the column headings remain consistent, so the tool can correctly identify and process the signatures.

The /template_for_adding_extension folder also contains a subfolder /01_inbuilt-signatures, which provides details of the original cell mappings.

If you want the extended signatures to be used in generating the consensus values, ensure that the cell_type names in your extension files (anno.1 and anno.2) match those in the Map_1 and Map_2 sheets, respectively.

##### Running with extension
``` 
df <- read.csv("rna_expression.csv")
sig_ext <- read.csv("sig.csv")
#### Optionally
anno.1_ext <- read.csv("anno.1.csv")
anno.2_ext <- read.csv("anno.1.csv")
#### Run decusus with extension
res <- cosDeco(x=df, rnaseq=T, plot=TRUE, ext=TRUE,
                    sig=sig_ext, anno.1=anno.1_ext, anno.2=anno.2_ext,
                    cp=NULL, free=FALSE) 
```

