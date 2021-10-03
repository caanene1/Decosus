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

  GeoMean <- function(x, na.rm=TRUE, zero.prop = FALSE){

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

      if (i == "house_keep") {
        # Calculate geometric mean
        # Ref:
        # PMID: 23810203, PMID: 28299307
        # PMID: 28275002
        res[[i]] <- apply(dgg[- symb], 2, GeoMean)

         } else {

        res[[i]] <- apply(dgg[- symb], 2, median)
      }
    }

    }

  res$HS <- apply(res[c(sig[sig != "house_keep"])], 1, mean)
  res$HS_N <- (res$HS / res$house_keep)

  return(res)

  ## HS = Sum(Median(S)) / n
}

