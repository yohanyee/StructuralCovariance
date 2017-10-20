#' Successive normalization of rectangular array
#' 
#' A method of successively normalizing both rows and columns of a matrix, a la Brad Efron and further described by Olshen et al.
#' See for more details: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2868388/ 
#' Note that the proof on convergence contained within aforementioned is flawed according to the 2013 paper published by the same authors
#' 
#' @param X A rectangular array
#' @param type Normalization method. This argument is passed as \code{type} to the \code{norm()}
#' @param tol Tolerance for convergence
#' @param verbose Be verbose 
#' @param na.set \code{norm()} requires that no elements in the rectangular array be \code{NA}. Replace \code{NA}s with this value
#' @return Normalized array
#' 
#' @export
norm_successive <- function(X, type="F", tol=1e-8, verbose=TRUE, na.set=0) {
  X[is.na(X)] <- na.set
  norm_types <- list(f="Frobenius norm", o="One norm", `1`="One norm", i="Infinity norm", m="Maximum modulus", `2`="Spectral (2)-norm")
  norm_difference <- 1
  norm_before <- norm(X, type = type)
  i <- 1
  while(norm_difference >= tol) {
    # Iterate
    rscale <- t(scale(t(X), center = TRUE, scale=TRUE)) # Row scale
    cscale <- scale(rscale, center = TRUE, scale=TRUE) # Column scale
    
    # Compute norm differences
    norm_after <- norm(cscale, type=type)
    norm_difference <- abs(norm_after - norm_before)
    norm_before <- norm_after
    X <- cscale
    
    # Print progress
    if (verbose) {
      print(paste("Iteration:", i, "|", norm_types[[tolower(type)]], "difference:", norm_difference))
      i <- i + 1
    }
  }
  return(X)
}

# Normalize 
#' @export
norm_rowsum <- function(X) {
  return(X/rowSums(X))
}

# Normalize 
#' @export
norm_colsum <- function(X) {
  return(X/colSums(X))
}

#' @export
norm_fullsum <- function(X) {
  return(X/sum(X))
}

#' @export
norm_rowregression <- function(X) {
  X_reg <- construct_like_matrix(X)
  X_rsum <- rowSums(X)
  for (i in 1:dim(X)[2]) {
    X_reg[,i] <- residuals(lm(X[,i] ~ X_rsum))
  }
  return(X_reg)
}

#' @export
norm_colregression <- function(X) {
  X_reg <- construct_like_matrix(X)
  X_csum <- colSums(X)
  for (i in 1:dim(X)[1]) {
    X_reg[i,] <- residuals(lm(X[i,] ~ X_csum))
  }
  return(X_reg)
}