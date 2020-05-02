#' @title Use quanTIseq to deconvolute a gene expression matrix.
#'
#' @description Function to apply the qunatiseq deconvolution method as described by bioRxiv 223180. https://doi.org/10.1101/223180.
#'
#'
#' @param mix.mat, arrays, mRNAscale, method,  defults mat, FALSE, TRUE, "LSEI
#'
#' @return t_results
#'
#' @keywords
#'
#' @examples See the original qunatiseq deconvolution method, Source code from https://github.com/FFinotello/quanTIseq
#'
#' @export
#'
quantiseq <- function(mix.mat, arrays = FALSE, mRNAscale = TRUE, method = "lsei") {
  # For the parameters,
  # use method = c("lsei"-limSolve or these "hampel", "huber", "bisquare" -robust regressions)
  # Always use lsei for quantiseq
  # set arrays = TRUE if it is array
  ## End of the parameter parts
  message("\nSetting helper functions for quanTIseq \n")
  #######################################################################################
  ### QuantiSEq helpers, see Source code from https://github.com/icbi-lab/immunedeconv
  fixMixture <- function(mix.mat, arrays = FALSE) {
    # Map gene names
    mix.mat <- mapGenes(mix.mat)
    # Un-log data in log2 base
    if (max(mix.mat) < 50) {
      mix.mat <- 2^mix.mat
    }
    # Quantile normalization
    if (arrays) mix.mat <- makeQN(mix.mat)
    # TPM normalization
    mix.mat <- t(t(mix.mat)*1e6/apply(mix.mat, 2, sum))
    return(mix.mat)
  }

  makeQN <- function(mix.mat) {
    cnames <- colnames(mix.mat)
    rnames <- rownames(mix.mat)
    mix.mat <- preprocessCore::normalize.quantiles(as.matrix(mix.mat))
    colnames(mix.mat) <- cnames
    rownames(mix.mat) <- rnames
    return(mix.mat)
  }

  mapGenes <- function(mydata) {
    HGNC <- read.csv(system.file("quantiseq", "HGNC_genenames_20170418.txt", package = "Decosus"),
                     header = TRUE, sep = "\t")
    curgenes <- rownames(mydata)
    newgenes <- rep(NA, length(curgenes))
    newgenes2 <- rep(NA, length(curgenes))
    ind <- match(curgenes, HGNC$ApprovedSymbol)

    # Current symbols and withdrawn ones
    genes.ind.notNA <- which(!is.na(ind))
    for (i in genes.ind.notNA) {
      genei <- curgenes[i]
      if (HGNC$Status[ind[i]] == "Approved") {
        newgenes[i] <- curgenes[i]
      } else if (HGNC$Status[ind[i]] == "EntryWithdrawn") {
        next
      } else {
        Wstring <- "symbolwithdrawn,see"
        newsymbol <- gsub(Wstring, "", HGNC$ApprovedName[ind[i]])
        newgenes2[i] <- newsymbol
      }
    }

    # Not found as symbols
    genes.ind.NA <- which(is.na(ind))
    for (i in genes.ind.NA) {
      genei <- curgenes[i]

      # Previos symbol?
      ind1 <- grep(genei, HGNC$PreviousSymbols)
      for (i1 in ind1) {
        array1 <- unlist(strsplit(as.character(HGNC$PreviousSymbols[i1]), ","))
        flag1 <- length(which(array1 == genei)) > 0
        if (flag1) {
          newsymbol <- as.character(HGNC$ApprovedSymbol[i1])
          newgenes2[i] <- newsymbol
        }
      }

      # Synonym?
      ind2 <- grep(genei, HGNC$Synonyms)
      for (i2 in ind2) {
        array2 <- unlist(strsplit(as.character(HGNC$Synonyms[i2]), ","))
        flag2 <- length(which(array2 == genei)) > 0
        if (flag2) {
          newsymbol <- as.character(HGNC$ApprovedSymbol[i2])
          newgenes2[i] <- newsymbol
        }
      }

    }
    newgenes2[which(newgenes2 %in% setdiff(newgenes,NA))] <- NA
    ind <- intersect(which(is.na(newgenes)),
                     which(!is.na(newgenes2)))
    newgenes[ind] <- newgenes2[ind]
    mydata <- mydata[which(!is.na(newgenes)),]
    newgenes <- newgenes[which(!is.na(newgenes))]
    # Take the median if duplicates are present
    outdata <- aggregate(mydata, by = list(newgenes), FUN = median)
    rownames(outdata) <- outdata[,1]
    outdata <- outdata[,-1, drop = FALSE]
    outdata <- as.data.frame(outdata)
    return(outdata)
  }

  quanTIseq <- function(currsig, currmix, scaling, method) {
    method <- match.arg(method, c("lsei", "hampel", "huber", "bisquare"))

    cgenes <- intersect(rownames(currsig), rownames(currmix))
    currsig <- as.matrix(currsig[cgenes,])
    currmix <- as.matrix(currmix[cgenes,])
    if (method == "lsei") {
      # Run deconvolution with constrained least squares
      G <- matrix(0, ncol = ncol(currsig), nrow = ncol(currsig))
      diag(G) <- 1
      G <- rbind(G, rep(-1, ncol(G)))
      H <- c(rep(0,ncol(currsig)),-1)
      results <- apply(currmix, 2, DClsei,
                       A = currsig, G = G, H = H,
                       scaling = scaling)
    } else {
      # Run deconvolution with robust regression
      results <- apply(currmix, 2, DCrr, A = currsig,
                       method = method, scaling = scaling)
    }

    #if (nrow(results)!=ncol(currmix))
    results <- t(results)
    return(results)

  }

  DClsei <- function(b, A, G, H, scaling) {
    sc <- norm(A, "2")
    A <- A/sc
    b <- b/sc
    res <- limSolve::lsei(A = A, B = b, G = G, H = H, verbose = FALSE)
    est <- res$X
    est.sum <- sum(est)
    est <- est/scaling
    est <- est/sum(est)*est.sum
    est <- c(est, pmax(0, 1 - sum(est)))
    names(est)[length(est)] <- "Other"
    return(est)
  }

  DCrr <- function(b, A, method, scaling){
    # Robust regression
    m <- paste0("psi.", method)
    if (m == "psi.hampel") {
      bres <- MASS::rlm(b ~ A, psi = m, a = 1.5, b = 3.5, c = 8, maxit = 1e3)
    } else {
      bres <- MASS::rlm(b ~ A, psi = m, maxit = 1e3)
    }
    est <- bres$coefficients
    # Remove intercept
    est <- est[-1]
    # Set negative values to 0
    est[est < 0] <- 0
    # Normalize total to 1=100%
    est <- est/sum(est)
    # Scale by mRNA content
    est.sum <- sum(est)
    est <- est/scaling
    est <- est/sum(est)*est.sum
    names(est) <- gsub("^A", "", names(est))
    return(est)
  }
  #######################################################################################

  message("\nRunning quanTIseq deconvolution module\n")

  if (is.numeric(mix.mat[[1,1]]) != TRUE) {
    stop("Wrong input format for the mixture matrix! Please follow the instructions of the documentation.")
  }

  # Load signature
  sig.mat <- read.table(system.file("quantiseq", "TIL10_signature.txt", package = "Decosus"),
                        header = TRUE, sep = "\t", row.names = 1)

  # Load normalization factors (set all to 1 if mRNAscale==FALSE)
  if (mRNAscale) {
    mRNA <- read.table(system.file("quantiseq", "TIL10_mRNA_scaling.txt", package = "Decosus"),
                       sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    colnames(mRNA) <- c("celltype", "scaling")

    mRNA <- as.vector(as.matrix(mRNA$scaling[match(colnames(sig.mat), mRNA$celltype)]))
  } else {
    mRNA <- rep(1, ncol(sig.mat))
  }

  # Preprocess mixture matrix
  message(paste0("Gene expression normalization and re-annotation (arrays: ",
                 arrays, ")\n"))
  mix.mat <- fixMixture(mix.mat, arrays = arrays)

  # Remove noisy genes
  lrmgenes <- as.vector(read.table(system.file("quantiseq", "TIL10_rmgenes.txt", package = "Decosus"),
                                   header = FALSE, sep = "\t")[,1])
  #
  n1 <- nrow(sig.mat)
  sig.mat <- sig.mat[!rownames(sig.mat) %in% lrmgenes,, drop = FALSE]
  n2 <- nrow(sig.mat)
  message(paste0("Removing ", n1 - n2, " noisy genes\n"))
  # End of removing noisy genes

  # Signature genes present in the mixture
  ns <- nrow(sig.mat)
  us <- length(intersect(rownames(sig.mat), rownames(mix.mat)))
  perc <- round(us*100/ns, digits = 2)
  message(paste0("Signature genes found in data set: ",
                 us, "/", ns, " (", perc, "%)\n"))

  # Run deconvolution
  message(paste0("Mixture deconvolution (method: ", method, ")\n"))
  results1 <- quanTIseq(sig.mat, mix.mat,
                        scaling = mRNA,
                        method = method)
  ## Sort Tregs and T.cells.CD4
  if ("Tregs" %in% colnames(sig.mat) && "T.cells.CD4" %in% colnames(sig.mat)
      && method %in% c("lsei")) {
    ####
    minTregs <- 0.02
    i <- which(colnames(sig.mat) == "T.cells.CD4")
    results2 <- quanTIseq(sig.mat[,-i],
                          mix.mat,
                          scaling = mRNA[-i],
                          method = method)

    ind <- which(results1[,"Tregs"] < minTregs)

    if (length(ind) > 0) {

      results1[ind,"Tregs"] <- (results2[ind,"Tregs"] + results1[ind,"Tregs"])/2
      results1[ind,"T.cells.CD4"] <- pmax(0, results1[ind, "T.cells.CD4"] -
                                            (results2[ind,"Tregs"] + results1[ind,"Tregs"])/2)
    }

  }
  results <- results1
  results <- results/apply(results, 1, sum)
  t_results <- as.data.frame(t(results))
  return(t_results)
}




