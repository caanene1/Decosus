## The code uses fmatch so install the package
## Don't load it as this is handled automatically
install.packages("fastmatch")
# You should now make the function avaialbale and test the dummy first

#' @title ssGSEA
#'
#' @description Function to apply ssGSEA
#'
#' @param exp Expression vector to score i.e one sample.
#'
#' @param gs The gene set (gs) to use to score the sample given in exp.
#'
#' @param dir The expected direction of the genes in gs i.e 1 is positive and -1 is negative.
#'
#' @param a Alpha to control weight given to the ranking of genes within gs.
#'          See readme for explantion our how to select the best alpha for your case. 0.25 is reasonable.
#'
#' @return Enrichment score of gs in exp using a.
#'
#' @keywords ssGSEA Enrchiment
#'
#' @examples See the readme file or get help with ?ssGSEA.
#'
#' @export
#'
ssGSEA <- function(exp, gs, dir=NULL, a=0.25) {

  # Sort the expression decreasing order
  sorted <- sort(exp, decreasing = T)
  genes <- names(sorted)

  # Check membership
  lookup <- fastmatch::fmatch(genes, gs)
  in_gs <- !is.na(lookup)

  # Rank and calculate weights
  r <- seq_along(sorted)
  w <- abs(r[in_gs]) ^a
  w_sum <- sum(w)

  # Normalise w_sum to avoid zero division
  w_sum <- ifelse(w_sum == 0, 1, w_sum)
  w <- w / w_sum

  # Set running sum vector
  r_sum <- numeric(length(genes))

  # Check if direction is provided and compute with direction
  if(!is.null(dir)) {
    # Match direction
    dir_in <- dir[which(gs %in% genes[in_gs])]
    # Assign unknown direction to postive
    dir_in[is.na(dir_in)] <- 1

    # Apply running sum with direction
    r_sum[in_gs] <- w * dir_in

    } else {
      # No direction defualt
      r_sum[in_gs] <- w
    }

  # Normalise
  r_sum[!in_gs] <- -1/ (length(genes) - sum(in_gs))
  # Compute enrichment
  es <- max(cumsum(r_sum))

    return(es)
  }



#' @title mssGSEA
#'
#' @description Function to apply ssGSEA on a matrix
#'
#' @param mexp Matrix of expression one sample per column and rownames as genes.
#'
#' @param lgs List og gene sets (gs) to use to score the sample given in mexp.
#'
#' @param ldir List of the expected direction of the genes in gs i.e 1 is positive and -1 is negative.
#'
#' @param a Alpha to control weight given to the ranking of genes within gs.
#'          See readme for explantion our how to select the best alpha for your case. 0.25 is reasonable.
#'
#' @return Enrichment scores per sample per gene list  using a.
#'
#' @keywords ssGSEA Enrchiment
#'
#' @examples See the readme file or get help with ?ssGSEA.
#'
#' @export
#'
mssGSEA <- function(mexp, lgs, ldir=NULL, a=0.5) {
 # Check data are correct
  # Check if input_data is a matrix
  if (!is.matrix(mexp) | !is.list(lgs)) {
    stop("Error: mexp must be a matrix. \n
         And lgs must be a list.")
  }

  if(!is.null(ldir)){
    if (!is.list(ldir)) {
      stop("Error: If ldir is provided, then it must be a list. \n
         The names in the list must match the names of lgs list.")
    }
  }

  # Run the implementation
  res <- lapply(names(lgs), function(gs) {
    # Set the direction
           if(!is.null(ldir)){
             dir_in <- ldir[[gs]]
           } else {
             dir_in <- ldir
           }

            lapply(colnames(mexp), function(x) {
              exp <- mexp[, x]
              ssGSEA(exp, lgs[[gs]], dir = dir_in, a=a)
    })
  })

  # Create the final results table
  names(res) <- names(lgs)
  res <- as.data.frame(do.call(cbind, res))
  res[] <- lapply(res, function(x) as.numeric(as.character(x)))
  rownames(res) <- colnames(mexp)

  return(res)
}



