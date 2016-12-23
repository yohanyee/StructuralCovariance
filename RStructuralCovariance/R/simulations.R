simulate_cormatrix_null_distribution <- function(v, method="pearson", quiet=FALSE) {
  vinfo <- get_volume_info(v)
  cmat <- matrix(NA, nrow=vinfo$structures$num, ncol=vinfo$structures$num, dimnames = list(vinfo$structures$names, vinfo$structures$names))
  
  #TODO
}