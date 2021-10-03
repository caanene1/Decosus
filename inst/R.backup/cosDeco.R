#' @title Create consensus deconvolution from "xcell, MCP, Danaher, Davoli, Rooney, quanTISeq, EPIC".
#'
#' @description Function to apply seven deconvolution/signature methods as described by the respective methods.
#'
#' @param x, platform, map, plot.corr: defaults: df, "Array", "Normal", FALSE
#'
#' @return t_results
#'
#' @keywords http://bio-bigdata.hrbmu.edu.cn/CellMarker/download.jsp
#'
#' @examples See the original methods.
#'
#' @export
#'
cosDeco <- function(x = df, platform = "Array") {
  plot.corr=FALSE

  #### Set variables ####
  # Negation
  `%notin%` <- Negate(`%in%`)

  # Curated Signatures
  y <- readxl::read_xlsx(system.file("extdata", "Signatures.xlsx",
                                     package = "Decosus"),
                         sheet = "Curated")
  print(dim(y))


  # if(in.CM == T){
  #  check <- dim(CM)[2]
  #  if (check == 3 & names(CM) %in% c("Type", "Gene", "Group")){
  #   y <- rbind(y, CM)
  # } else {
  #    print("Fool someone else, get the right data")
  #   print("Cm need to be 3 columns named Type, Gene, Group")
  #  } }


  x_class <- sapply(x, class)

  # Get the single factor column as gene names
  col_factor <- x_class[x_class %in%  c("character", "factor")]
  if (!length(col_factor) == 1) {
    stop("Expression matrix must have a column of hgnc symbols")
  }

  # Get and rename symbol column
  HGNC <- which(colnames(x) == names(col_factor))
  names(x)[HGNC] <- "Gene"

  # Get the expression of gene sets
  x_sub <- x[x$Gene %in% y$Gene, ]
  #### End of variable handling ####

  #### Handling of gene sets ####
  # Get the unique sources
  soruce <- sapply(unique(y["Group"]),
                   function(x) paste(as.character(x)))
  print(soruce)

  # Separate by sources
  gene_sets <- lapply(soruce,
                      function(x){subset(y, Group == x)})

  # Merge set with expression corresponding matrix
  x_sub_sets <- lapply(gene_sets,
                       function(x){merge(x, x_sub, by = "Gene",  all = F)})
  #### Ends of gene sets handling ####

  #### Aggregation ####
  #### Anthony Look art this again +++++
  # Force numeric for NA values
  mean_convert <- function(x) {
    mean(as.numeric(as.character(x))) }

  # Subset only numeric variables
  numeric_columns <- function(x) {
    num_col <- sapply(x, is.numeric)
    return(x[ , num_col]) }

  # Aggregate expression values by Cell (Abundance)
  aggregate_expr <- function(x) {
    x_num <- numeric_columns(x)
    print(dim(x_num))
    aggregated <- aggregate(x_num, by = list(x$Type),
                            FUN = mean_convert)
    names(aggregated)[1] <- "Cell"
    aggregated$Source <- unique(x$Group)
    return(aggregated) }
  #### End of aggregation ####

  # Cell proportions curated signatures
  deconvoluted <- lapply(x_sub_sets, aggregate_expr)
  #
  gc()

  ############## Cell deconvolution methods ###################
  # Assign new data-frame
  x_immdeconv <- x
  # set row_names
  rownames(x_immdeconv) <- x_immdeconv$Gene
  x_immdeconv <- x_immdeconv[,-which(names(x_immdeconv) %in% "Gene")]




  ## xCell ##
  # TO DO: Consider including xCell different way.
  library(xCell)
  # requireNamespace("xCell")
  if (platform != "array") {
    xcell <- as.data.frame(xCell::xCellAnalysis(x_immdeconv, rnaseq = TRUE))
  } else {
    xcell <- as.data.frame(xCell::xCellAnalysis(x_immdeconv, rnaseq = FALSE))
  }
  #
  xcell$Source <- "xcell"
  ## End of xCell ##
  ##
  if (length(colnames(xcell)) == 0) {
    stop("xcell analysis failed consider installing the right version at rdrr.io/github/alex-pinto/xc/")
  }



  ##
  ## MCP ##
  MCP <- as.data.frame(MCPcounter::MCPcounter.estimate(x_immdeconv,
                       featuresType = "HUGO_symbols",
                       genes = read.table(system.file("MCPcounter", "genes.txt", package = "Decosus"),
                                          sep = "\t", stringsAsFactors = FALSE,
                                          header = TRUE, colClasses = "character",
                                          check.names = FALSE)))
  MCP$Source <- "MCP"
  ## End of MCP ##

  # Suppress warning for EPIC
  options(warn=-1)
  ## EPIC ##
  Epic <- as.data.frame(t(EPIC::EPIC(bulk = x_immdeconv)[[2]]))
  options(warn=0)
  Epic$Source <- "EPIC"
  ## End of EPIC ##

  ## quanTISeq ##
  if (platform == "arrays") {
    arrays = TRUE
  } else {
    arrays = FALSE
  }
  #
  quanTISeq <- Decosus::quantiseq(x_immdeconv, arrays = arrays,
                         mRNAscale = TRUE, method = "lsei")
  quanTISeq$Source <- "quanTISeq"
  ## End of quanTISeq ##

  ####### Sort the names and bind it for processing #######
  # Cell names
  fun_insert_cell_names <- function(d) {
    d1 <- cbind(rownames(d), d)
    names(d1)[1] <- "Cell"
    rownames(d1) <- NULL
    return(d1) }
  #
  deconvoluted2 <- lapply(list(xcell, MCP, Epic, quanTISeq),
                          fun_insert_cell_names)
  ############## End of cell deconvolution methods ###################

  ############## Building the conseus ###################
  all_results <- append(deconvoluted, deconvoluted2)
  ########
  ######## Remove garbage, for low power computers ########
  rm(y, x_class, col_factor, HGNC, x_sub, soruce, gene_sets, x_sub_sets,
     deconvoluted, x_immdeconv, xcell, MCP, Epic, quanTISeq, deconvoluted2)


  ################# ############### ############## #######################
  # Bind to frame
  all_results = as.data.frame(data.table::rbindlist(all_results))
  all_results$ID <- paste(all_results$Cell, all_results$Source, sep = "_")

  ##### Load and process annotation #######
  fun_l_anno <- function(k, xlSheet) {
  anno <- readxl::read_xlsx(system.file("extdata", "Signatures.xlsx",
                                        package = "Decosus"), sheet = xlSheet)
  ## Create match IDs
  anno$ID <- paste(anno$Source_cell_type, anno$Source, sep = "_")
  ## Bind
  a_res <- merge(k, anno, by = "ID")
  #
  sub_a_res <- subset(a_res, !is.na(a_res$cell_type))
  return(list(a_res, sub_a_res)) }


  ### Updated for the modified cell map 29/01/2021
  for_samples <- lapply(list("Map_1", "Map_1b", "Map_1c"),
                        function(x) {
                          fun_l_anno(k = all_results, xlSheet = x)
                          })

  #
  for_cells <- lapply(list("Map_2", "Map_2b", "Map_2c"),
                        function(x) {
                          fun_l_anno(k = all_results, xlSheet = x)
                        })

  ################# ############### ############## ######################
  ########
  ### Aggregate the methods ###
  aggregate_cell <- function(x) {
    # See fun, numeric_columns, mean_convert, above
    x_num <- numeric_columns(x)
    aggregated <- aggregate(x_num, by = list(x$cell_type),
                            FUN = mean_convert)
    names(aggregated)[1] <- "Cell"
    aggregated$Source <- "consensus"
    return(aggregated) }

  # Function consensus available
  avilab <- function(f) {
  sub <- f
  con_score <- aggregate_cell(sub)
  rownames(con_score) <- paste(con_score$Cell, con_score$Source, sep = "_")
  con_score <- con_score[ , which(names(con_score) %notin% c("Cell","Source"))]
  # Match the sub_let
  rownames(sub) <- sub$ID
  sub <- sub[ , which(names(sub) %notin%
                        c("Source.x","Source.y",
                          "Source_cell_type",
                          "full_annotation",
                          "cell_type", "Cell", "ID"))]
  #
  con_avil <- rbind(sub, con_score)
  outavil <- list(con_score, con_avil)
  return(outavil)
  }

  # Function bind consensus and non-consensus
  cons_nonCons <- function(u, z) {
    sub_non <- subset(u, is.na(u$cell_type))
    rownames(sub_non) <- sub_non$ID
    sub_non <- sub_non[ , which(names(sub_non) %notin%
                                  c("Source.x", "Source.y",
                                    "Source_cell_type",
                                    "full_annotation",
                                    "cell_type", "Cell", "ID"))]
    res_f <- rbind(z, sub_non)
    return(res_f) }



  ########### Move this function to a location to make it first class
  ########### Then call it here, then use it in the pro.res function
  # Function plot correlation
  cor_pp <- function(p, pdf.name) {
  # Set the custom colour scheme
  col <- colorRampPalette(c("steelblue4", "skyblue1", "white", "tomato","darkred"))
  # Preprocess
  cor_mcons <- cor(as.data.frame(t(p)))
  # Cell type colour
  coCo <- data.frame(Name = colnames(cor_mcons))
  b_color <- readxl::read_xlsx(system.file("extdata", "Signatures.xlsx", package = "Decosus"),
                                   sheet = "Colour")
  #
  coloo <- merge(coCo, b_color, by = "Name", sort = F)
  mycolors <- coloo$Colour
  names(mycolors) <- coloo$Name
  #
  ord_mycolors <- mycolors[corrplot::corrMatOrder(cor_mcons,
                                                  order = "hclust",
                                                  hclust.method = "complete")]

  # Prints PDF to working directory
  pdf(file = pdf.name)
  corrplot::corrplot(cor_mcons,
                     col = col(100),
                     method = "circle",
                     outline = FALSE,
                     pch.cex = 0.03,
                     tl.cex = 0.25,
                     mar = c(3,3,3,3),
                     cl.ratio = 0.09,
                     cl.cex = 0.6,
                     tl.col = ord_mycolors,
                     type = "lower",
                     order = "hclust",
                     hclust.method = "complete",
                     tl.srt = 60)
  dev.off()
  }

  # Plot Here for the call above
  if(plot.corr){
  cor_pp(Sample_c[[2]], "Sample_consensus.pdf")
  cor_pp(Cell_c[[2]], "Cell_consensus.pdf")
  }
  ######### End of plot correlations ##########

  ### Build output
  {
    main_samples <- lapply(for_samples, function(x) {
      cons_nonCons(u = x[[1]], z = avilab(x[[2]])[[1]])
    })

    names(main_samples) <- c("Map_1", "Map_1b", "Map_1c")

    main_cells <- lapply(for_cells, function(x) {
      cons_nonCons(u = x[[1]], z = avilab(x[[2]])[[1]])
    })
    names(main_cells) <- c("Map_2", "Map_2b", "Map_2c")
  }

  output <- list(main_samples = main_samples,
                 main_cells = main_cells,
                 raw_results = all_results)
  return(output)
}
