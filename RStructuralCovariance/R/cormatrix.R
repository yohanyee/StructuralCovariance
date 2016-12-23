#' Construct a correlation matrix by correlating volumes
#' 
#' Creates a correlation matrix by computing pairwise correlations between volumes (column vectors of a \code{data.frame})
#' 
#' @param v volume \code{data.frame} consisting of individuals as rows and volumes as columns
#' @param method character string giving the type of correlation coefficient being computed
#' @param p.adjust character string giving method to use when correcting p values for multiple comparisons
#' @param normalize logical value denoting whether volumes should be normalized (rows summed to 1)
#' @param resample logical value denoting whether to resample individuals from the volume \code{data.frame}
#' @param resample_N integer value between 3 and the number of rows of the volume \code{data.frame} 
#' @param resample_replace logical value denoting whether to resample with replacement
#' @param quiet logical value denoting whether progress should be printed
#' 
#' @return a cormatrix object
#' @export
correlate_volumes <- function(v, cor.method="pearson", p.adjust.method="fdr", normalize=TRUE, resample=FALSE, resample_N=NULL, resample_replace=TRUE, quiet=FALSE) {
  mice_num <- dim(v)[1]
  struc_num <- dim(v)[2]
  struc_names <- colnames(v)
  if (resample) {
    if (is.null(resample_N)) {
      resample_N <- mice_num
    }
    resample_N <- max(3,min(mice_num, resample_N))
    mice_num <- resample_N
    indices <- sample(1:mice_num, size = resample_N, replace = resample_replace)
    v <- v[indices,]
  }
  if (normalize) {
    v <- v/rowSums(v)
  }
  
  mtx <- matrix(NA, nrow=struc_num, ncol=struc_num, dimnames = list(struc_names, struc_names))
  mtx_p <- matrix(NA, nrow=struc_num, ncol=struc_num, dimnames = list(struc_names, struc_names))
  
  # Compute diagonals
  for (i in 1:struc_num) {
    ct <- cor.test(v[,i], v[,i], cor.method=method)
    mtx[i,i] <- ct$estimate
    mtx_p[i,i] <- ct$p.value
  }
  if (!quiet) {
    cat("1 ")
  }
  
  # Compute off-diagonals
  for (i in 2:struc_num) {
    for (j in 1:(i-1)) {
      ct <- cor.test(v[,i], v[,j], cor.method=method)
      mtx[i,j] <- mtx[j,i] <- ct$estimate
      mtx_p[i,j] <- mtx_p[j,i] <- ct$p.value
    }
    if (!quiet) {
      cat(paste(i, " ", sep=" ")) 
    }
  }
  
  
  # Correct for multiple comparisons
  p_vals <- p.adjust(mtx_p[lower.tri(mtx_p, diag=FALSE)], method = p.adjust.method)
  mtx_p[lower.tri(mtx_p, diag=FALSE)] <- p_vals
  mtx_p[upper.tri(mtx_p, diag=FALSE)] <- t(mtx_p)[upper.tri(mtx_p, diag=FALSE)]
  if (!quiet) {
    cat("Done.\n")
  }
  
  # Return
  out <- list(cor=mtx, p=mtx_p)
  return(out)
}
