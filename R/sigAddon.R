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

  in.file <- read.delim(path)
  in.file <- in.file[!in.file$markerResource %in% "Single-cell sequencing", ]
  in.file <-in.file[c("cellName", "geneSymbol")]

  # i <- uni.cell[1]
  out.file <- data.frame(Type=character(),
                         Gene=character())

  for(i in unique(in.file$cellName)){
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
  out.file$Group <- "CMdb"
  out.file$filter <- paste(out.file$Type, out.file$Gene)
  out.file <- out.file[!duplicated(out.file$filter), ]
  out.file$is.gene <- ifelse(grepl("[a-z]", out.file$Gene), "no", "yes")
  out.file <- out.file[out.file$is.gene %in% "yes", ]
  out.file <- out.file[!out.file$Gene %in% c("NA", "II", "I", "MHC"), ]
  ##
  tcount <- as.data.frame(table(out.file$Type))
  tcount <- tcount[tcount$Freq >= min.set,]
  ##
  out.file <- out.file[out.file$Type %in% tcount$Var1, ]


  return(list(out.file[c("Type", "Gene", "Group")], tcount, in.file))
}
