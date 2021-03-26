#' @title Create consensus hypoxia score based on Buffa, Winter, Ragnum, Eustace, Sorensen, Elvidge, Hu
#'
#' @description Function to apply seven hypoxia signature and generate a hypoxia score HS
#'
#' @param x
#'
#' @return HS and individual methods
#'
#' @keywords Hypoxia
#'
#' @examples See the original methods.
#'
#' @export
#'
HypoxiaScore <- function(x = df) {

  y <- readxl::read_xlsx(system.file("extdata", "Hypoxia_gene.xlsx",
                                     package = "Decosus"),
                         sheet = "marker")

  sig <- unique(y$Soruce)

  x_class <- sapply(x, class)
  col_factor <- x_class[x_class %in%  c("character", "factor")]
  if (!length(col_factor) == 1) {
    stop("Expression matrix must have a column of hgnc symbols")
  }

  symb <- which(colnames(x) == names(col_factor))
  names(x)[symb] <- "Gene"

  {
    res <- data.frame("ID" = names(x[- symb]))
    for(i in sig) {
      gg <- y[y$Soruce %in% i, ][["Gene"]]
      dgg <- x[x$Gene %in% gg, ]
      res[[i]] <- apply(dgg[- symb], 2, median) }
  }

  res$HS <- apply(res[-1], 1, mean)

  return(res)

  ## HS = Sum(Median(S)) / n
}

