#' @title Process new signatures for inclusion into the analysis
#'
#' @description Function to process signatures from the CellMarker database
#'
#' @param path, is.web, min.set
#'
#' @return Processed table similar to the table at signature.xlsx
#'
#' @keywords CellMarker
#'
#' @examples See the original methods.
#'
#' @export
#'
sigAddon <- function(path = url.in, is.web = T, min.set=3){

  if(is.web){
    print("Web import!")
  } else {
    print("Local import!")
  }

  in.file <- read.delim(url.in)
  in.file$Group <- paste(in.file$cancerType, in.file$tissueType,
                         sep = "_")

  in.file <-in.file[c("cellName", "geneSymbol")]
  uni.cell <- unique(in.file$cellName)

  i <- uni.cell[1]
  out.file <- data.frame(Type=character(),
                         Gene=character())

  for(i in uni.cell){
    c.ii <- in.file[in.file$cellName == i, ]
    c.ii <- c.ii$geneSymbol
    c.ii <- lapply(c.ii, function(z){
      z <- strsplit(z, " +")
    })
    c.ii <- as.character(unlist(c.ii))

    c.ii <- data.frame(Type = i,
                       Gene = c.ii)

    if(nrow(c.ii >= min.set)){
      out.file <- rbind(out.file, c.ii)
    }


  }

  out.file$Gene <- sub("\\[|\\]", "", out.file$Gene)
  out.file$Gene <- sub(",", "", out.file$Gene)
  out.file$Group <- "CM"

  filterF <- paste(out.file$Type, out.file$Gene, sep = "_")
  out.file <- out.file[!duplicated(filterF), ]

  return(list(out.file, in.file))
}


