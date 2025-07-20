# Class
base.class = "data.frame"
`%notin%` <- Negate(`%in%`)

mdl.p <- list(xCell = base.class,
              MCP = base.class,
              quanTISeq = base.class,
              EPIC = base.class,
              other = base.class)

setClass("DecoCell", slots = list(e.data = base.class,
                                  p.data = "list",
                                  s.data = "list",
                                  sig1 = base.class,
                                  sig2 = base.class,
                                  sig3 = base.class,
                                  anno1 = "list",
                                  anno2 = "list",
                                  anno3 = "list",
                                  c.data = base.class,
                                  consesus = "list",
                                  res.final = "list"))

setClass("pDeco", slots = mdl.p, contains = "DecoCell")
setClass("sDeco", slots = list(order = base.class), contains = "DecoCell")


# Generics
setGeneric("addExpression", function(x, y) standardGeneric("addExpression"))

setGeneric("addSignature", function(x, ext, sig, anno.1, anno.2) standardGeneric("addSignature"))

setGeneric("deconvolveCell", function(x, y, rnaseq, free, scale.i) standardGeneric("deconvolveCell"))

setGeneric("deconvoleConsensus", function(x, meth) standardGeneric("deconvoleConsensus"))


# Methods
setMethod("addExpression", "DecoCell", function(x, y) {

  if(!is(object = x, class2 = 'DecoCell')) {

    stop('object must be of class DecoCell, pDeco or sDeco\n')
  }

  y_class <- sapply(y, class)

  col_factor <- y_class[y_class %in%  c("character", "factor")]
  if (!length(col_factor) == 1) {
    stop("ERROR: Expression matrix must have JUST one column of hgnc symbols")
  }

  HGNC <- which(colnames(y) == names(col_factor))
  names(y)[HGNC] <- "Gene"

  y <- y[!duplicated(y$Gene),]
  rownames(y) <- y$Gene

  x@e.data <- y
  x
})

setMethod("addSignature", "DecoCell", function(x, ext=TRUE, sig=NULL,
                                               anno.1=NULL, anno.2=NULL) {

  if(!is(object = x, class2 = 'DecoCell')) {

    stop('object must be of class DecoCell, pDeco or sDeco\n')
  }

  if(nrow(x@sig1) == 0){
    x@sig1 <- as.data.frame(
      readxl::read_xlsx(system.file("extdata", "Signatures.xlsx",
                                    package = "Decosus"),
                        sheet = "Curated"))

    x@anno1 <- list(as.data.frame(readxl::read_xlsx(system.file("extdata", "Signatures.xlsx",
                                                                package = "Decosus"),
                                                    sheet = "Map_1")),

                    as.data.frame(readxl::read_xlsx(system.file("extdata", "Signatures.xlsx",
                                                                package = "Decosus"),
                                                    sheet = "Map_2"))) }


  if(ext){
    if(sum(ifelse(names(sig) %in% names(x@sig1), 0, 1)) != 0 |
       sum(ifelse(names(anno.1) %in% names(anno.2), 0, 1)) != 0 |
       sum(ifelse(names(anno.1) %in% names(x@anno1[[1]]), 0, 1)) != 0){

      print("WARNING: Wrong signature or annotation file. \n
              The ext argument is not used.")
      print("Extension must contain:
              --Signature -Type -Gene -Group \n
              --Annotation Source Source_cell_type cell_type full_annotation")
      print("See https://github.com/caanene1/Decosus")

    } else {

      if(nrow(x@sig2) == 0){
        x@sig2 <- sig[c("Type", "Gene", "Group")]
        x@anno2 <- list(anno.1[c("Source", "Source_cell_type",
                                 "cell_type", "full_annotation")],
                        anno.2[c("Source", "Source_cell_type",
                                 "cell_type", "full_annotation")])

      } else if(nrow(x@sig3) == 0) {
        x@sig3 <- sig[c("Type", "Gene", "Group")]
        x@anno3 <- list(anno.1[c("Source", "Source_cell_type",
                                 "cell_type", "full_annotation")],
                        anno.2[c("Source", "Source_cell_type",
                                 "cell_type", "full_annotation")])

      } else {
        stop("Only three slots allowed for signatures. \n
         Remove old slots and try again.")
      }
    }
  }
  x
})

setMethod("deconvolveCell", "DecoCell", function(x, y, rnaseq, free=FALSE,
                                                 scale.i=T) {

  if(!is(object = x, class2 = 'DecoCell')) {
    stop('object must be of class DecoCell, pDeco or sDeco\n')
  }


  if (nrow(x@e.data) == 0){
    stop("Error: No expression data found. Call addExpression() first.")
  }

  xx <- x@e.data[,-which(names(x@e.data) %in% "Gene")]

  ##
  if(scale.i){
    df.scale <- function(df="data") {
      ## Scale the outputs individually
      ## Don't centre to avoid negative values

      out <- t(apply(df, 1, function(x) scale(x, center=FALSE)))
      out <- as.data.frame(out)

      names(out) <- names(df)

      return(out)
    }
    ##
  }


  if(y == "p"){
    library(xCell)
    x@xCell <- as.data.frame(xCell::xCellAnalysis(xx, rnaseq=rnaseq))
    if (length(colnames(x@xCell)) == 0) {
      print("xcell failed, see rdrr.io/github/alex-pinto/xc/") }


    x@MCP <- as.data.frame(MCPcounter::MCPcounter.estimate(xx, featuresType = "HUGO_symbols",
                                                           genes = read.table(system.file("MCPcounter", "genes.txt",
                                                                                          package = "Decosus"),
                                                                              sep = "\t", stringsAsFactors = FALSE,
                                                                              header = TRUE, colClasses = "character",
                                                                              check.names = FALSE)))


    if(!free){
      options(warn=-1)
      x@EPIC <- as.data.frame(t(EPIC::EPIC(bulk = xx)[[2]]))
      options(warn=0)

      if(scale.i){
        x@EPIC <- df.scale(x@EPIC)
      }

      x@EPIC$Source <- "EPIC"

    }

    x@quanTISeq <- Decosus::quantiseq(xx, arrays = !rnaseq,
                                      mRNAscale = TRUE,
                                      method = "lsei")



    if(scale.i){
      x@xCell <- df.scale(x@xCell)
      x@MCP <- df.scale(x@MCP)
      x@quanTISeq <- df.scale(x@quanTISeq)
      }


    x@xCell$Source <- "xcell"
    x@MCP$Source <- "MCP"
    x@quanTISeq$Source <- "quanTISeq"


    x@p.data <- lapply(list(x@xCell, x@MCP, x@EPIC, x@quanTISeq), function(z) {
      if(nrow(z) >= 0){
        z <- cbind(rownames(z), z)
        names(z)[1] <- "Cell"
        rownames(z) <- NULL
      }
      return(z) })

  } else if(y == "s") {

    signature <- do.call(rbind, list(x@sig1, x@sig2, x@sig3))

    x_sub <- x@e.data[x@e.data$Gene %in% signature$Gene, ]
    soruce <- sapply(unique(signature["Group"]),
                     function(k) paste(as.character(k)))

    gene_sets <- lapply(soruce,function(k){
      subset(signature, Group == k)})

    x_sub_sets <- lapply(gene_sets, function(k){
      merge(k, x_sub, by = "Gene",  all = F)})

    mean_convert <- function(x) {
      mean(as.numeric(as.character(x))) }

    numeric_columns <- function(x) {
      num_col <- sapply(x, is.numeric)
      return(x[ , num_col]) }

    aggregate_expr <- function(x, sca=scale.i) {
      x_num <- numeric_columns(x)
      aggregated <- aggregate(x_num, by = list(x$Type),
                              FUN = mean_convert)
      names(aggregated)[1] <- "Cell"

      if(sca){
        aggregated <- cbind(aggregated[1], df.scale(aggregated[-1]))
      }

      aggregated$Source <- unique(x$Group)
      return(aggregated) }

    x@s.data <- lapply(x_sub_sets, aggregate_expr)

  } else {
    print("Y argument is not recorgnised")
    # Place holder for future updates
  }

  x
})

setMethod("deconvoleConsensus", "DecoCell", function(x, meth="geomean") {

  if(!is(object = x, class2 = 'DecoCell')) {
    stop('object must be of class DecoCell, pDeco or sDeco\n')
  }

  if(length(x@p.data) + length(x@s.data) == 0){
    stop("Please, run deconvovleCell first")
  }

  x@c.data <- do.call(rbind, append(x@p.data, x@s.data))
  x@c.data$ID <- paste(x@c.data$Cell, x@c.data$Source, sep="_")

  process.anno <- function(k, anno) {
    anno$ID <- paste(anno$Source_cell_type, anno$Source, sep = "_")
    a_res <- merge(k, anno, by = "ID")
    #
    sub_a_res <- subset(a_res, !is.na(a_res$cell_type))
    return(list(a_res, sub_a_res)) }

  anno.count <- length(x@anno2) + length(x@anno3)

  if(anno.count == 4){
    annotate1 <- do.call(rbind, list(x@anno1[[1]], x@anno2[[1]], x@anno3[[1]]))
    annotate2 <- do.call(rbind, list(x@anno1[[2]], x@anno2[[2]], x@anno3[[2]]))
  } else if(anno.count == 2){
    annotate1 <- do.call(rbind, list(x@anno1[[1]], x@anno2[[1]]))
    annotate2 <- do.call(rbind, list(x@anno1[[2]], x@anno2[[2]]))
  } else {
    annotate1 <- x@anno1[[1]]
    annotate2 <- x@anno1[[2]]
  }

  for_samples <- process.anno(k=x@c.data, anno = annotate1)
  for_cells <- process.anno(k=x@c.data, anno = annotate2)

  ######## Methods
  mean_convert2 <- function(x) {
    mean(as.numeric(as.character(x))) }

  geo.mean <- function(x, na.rm=TRUE, zero.prop = FALSE){

    if (any(x < 0, na.rm = TRUE)) {
      return(NaN)
    }
    if (zero.prop) {
      if (any(x == 0, na.rm = TRUE)) {
        return(0)
      }
      exp(mean(log(x), na.rm = na.rm))
    } else {
      exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
    }
  }
  ########

  numeric_columns2 <- function(x) {
    num_col <- sapply(x, is.numeric)
    return(x[ , num_col]) }

  aggregate_cell <- function(x, me=meth) {
    x_num <- numeric_columns2(x)

    if(me == "mean") {
      aggregated <- aggregate(x_num, by = list(x$cell_type),
                              FUN = mean_convert2)

    } else if (me == "geomean") {
      aggregated <- aggregate(x_num, by = list(x$cell_type),
                              FUN = geo.mean)
    } else {
      stop("Provide correct method argument, either mean or geomean")
    }

    names(aggregated)[1] <- "Cell"
    aggregated$Source <- "consensus"
    return(aggregated) }

  avilab <- function(f) {
    sub <- f
    con_score <- aggregate_cell(sub)
    rownames(con_score) <- paste(con_score$Cell, con_score$Source, sep = "_")
    con_score <- con_score[ , which(names(con_score) %notin% c("Cell","Source"))]

    rownames(sub) <- sub$ID
    sub <- sub[ , which(names(sub) %notin%
                          c("Source.x","Source.y",
                            "Source_cell_type",
                            "full_annotation",
                            "cell_type", "Cell", "ID"))]
    #
    con_avil <- rbind(sub, con_score)
    outavil <- list(con_score, con_avil)
    return(outavil) }

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


  {
    main_samples <- cons_nonCons(u=for_samples[[1]],
                                 z=avilab(for_samples[[2]])[[1]])

    main_cells <- cons_nonCons(u=for_cells[[1]],
                               z=avilab(for_cells[[2]])[[1]])
  }

  x@res.final <- list(main_samples = main_samples,
                      main_cells = main_cells)
  x
})



#' @title Correlation plot for data frame
#'
#' @description Function to plot clean and coloured correlation plot.
#'
#' @param p, cp, pdf.name:
#'
#' @return plot saved to working directory.
#'
#' @keywords plotting
#'
#' @examples cor.pp(p, cp, pdf.name)
#'
#' @export
#'
cor.pp <- function(p, cp=NULL, pdf.name) {

  col <- colorRampPalette(c("steelblue4", "skyblue1",
                            "white", "tomato","darkred"))

  cor_mcons <- cor(as.data.frame(t(p)))
  cor_mcons[is.na(cor_mcons)] <- 0
  coCo <- data.frame(Name = colnames(cor_mcons))

  if(is.null(cp)){
    ord_mycolors <- NULL

  } else {

    if(nrow(cp) >= 1) {
    coloo <- merge(coCo, cp, by = "Name", sort = F)
    mycolors <- coloo$Colour
    names(mycolors) <- coloo$Name
    ord_mycolors <- mycolors[corrplot::corrMatOrder(cor_mcons,
                                                    order = "hclust",
                                                    hclust.method = "complete")]
    }
  }


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


#' @title Create consensus deconvolution
#'
#' @description
#' Applies seven deconvolution or gene signature-based estimation methods
#' to derive consensus estimates of cell-type or signature proportions from expression data.
#' Combines results for robustness and flexibility in downstream analyses.
#'
#' @param df Expression data frame (genes × samples), only one column of none numeric allowed.
#'          This is the main input.
#' @param rnaseq Logical flag indicating whether the input expression matrix is RNA-seq data.
#'               Affects normalization and method selection. Default is `TRUE`.
#' @param sig Data frame of extension signature matrix to be combined with the in-built signatures.
#'             It must have three columns named "Type	Gene	Group"
#' @param anno.1 Map_1 data frame for the signatures in the extension given as "sig" (relative consensus).
#' @param anno.2 Map_2 data frame for the signatures in the extension given as "sig" (absolute consensus).
#' @param cp Logical; whether to generate correlation plots during processing.
#' @param plot Logical; whether to display plot. Default is `FALSE`.
#' @param free Logical; if `TRUE`, allows unconstrained (or relaxed) deconvolution when applicable.
#'             This may apply to certain optimization steps.
#' @param mini.output Logical; if `TRUE`, returns a simplified output with only essential results.
#'                    Default is `FALSE`.
#'
#' @return A list containing:
#' \describe{
#'   \item{consensus}{Consensus deconvolution result (matrix of cell-type proportions)}
#'   \item{methods}{Results from each individual method}
#'   \item{params}{Metadata and parameters used}
#' }
#'
#' @keywords deconvolution, gene expression, cell type proportion, signature matrix, consensus
#'
#' @examples
#' # Basic usage with expression matrix:
#' cosDeco(x = exp)
#'
#' # Using a custom signature matrix and plotting enabled:
#' cosDeco(x = exp, sig = my_sig, plot = TRUE)
#'
#' @export
cosDeco <- function(x=df, rnaseq=T, ext=FALSE, sig=NULL, anno.1=NULL,
                    anno.2=NULL, cp=NULL, plot=TRUE, free=FALSE, scale.i=T,
                    agg.method="mean", mini.output=TRUE) {

  output <- new("pDeco")
  output <- addExpression(output, x)
  output <- addSignature(output, ext=ext, sig=sig,
                         anno.1=anno.1, anno.2=anno.2)

  output <- deconvolveCell(output, "p",  rnaseq=rnaseq, free=free, scale.i=scale.i)
  output <- deconvolveCell(output, "s",  rnaseq=rnaseq, scale.i=scale.i)
  output <- deconvoleConsensus(output, meth=agg.method)


  if(plot){
    if(is.null(cp)){
      cp.in <- readxl::read_xlsx(system.file("extdata", "Signatures.xlsx",
                                          package="Decosus"),
                              sheet = "Colour")
    } else {
      cp.in <- cp
    }

    cor.pp(output@res.final[["main_samples"]], cp=cp.in,
           pdf.name = "Sample_consensus.pdf")

    cor.pp(output@res.final[["main_cells"]], cp=cp.in,
           pdf.name = "Cell_consensus.pdf")
  }


  if(mini.output){
    return(output@res.final)
  } else {
    return(output)
  }

}
