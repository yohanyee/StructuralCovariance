#' @export
get_volume_info <- function(v) {
  subjects <- list(num=dim(v)[1], names=rownames(v))
  structures <- list(num=dim(v)[2], names=colnames(v))
  out <- list(subjects=subjects, structures=structures)
  return(out)
}

#' @export
resample_volumes <- function(v, size=NULL, replace=TRUE, prob=NULL) {
  vinfo <- get_volume_info(v)
  size <- (if (is.null(size)) vinfo$subjects$num else size)
  size <- (if (size > vinfo$subjects$num | size < 3) vinfo$subjects$num else size)
  indices <- sample(1:size, size=size, replace=replace, prob=prob)
  v <- v[indices,]
  return(v)
}

#' @export
normalize_volumes <- function(v) {
  return(v/rowSums(v))
}

#' @export
correlate_volumes <- function(v, method="pearson", alternative="two.sided", quiet=FALSE, ...) {
  vinfo <- get_volume_info(v)
  cmat <- pmat <- matrix(NA, nrow=vinfo$structures$num, ncol=vinfo$structures$num, dimnames = list(vinfo$structures$names, vinfo$structures$names))
  
  # Compute diagonals
  if (!quiet) {cat("1 ")}
  diag(cmat) <- 1
  diag(pmat) <- 0

  # Compute off-diagonals
  for (i in 2:vinfo$structures$num) {
    if (!quiet) {cat(paste(i, " ", sep=""))}
    for (j in 1:(i-1)) {
      ct <- cor.test(v[,i], v[,j], method=method, alternative=alternative, ...)
      cmat[i,j] <- cmat[j,i] <- ct$estimate
      pmat[i,j] <- pmat[j,i] <- ct$p.value
    }
  }
  cat("Done.\n")
  
  # Coerce to RSC matrix objects
  cmat <- as.RSC.cmatrix(cmat)
  pmat <- as.RSC.pmatrix(pmat)
  
  cormatrix <- as.RSC.cormatrix(list(correlation=cmat, p.value=pmat))
  return(cormatrix)
}

#' @export
p.adjust <- function(p, method="fdr", n=length(p)) {
  UseMethod("p.adjust", p)
}

#' @export
p.adjust.default <- function(p, method="fdr", n=length(p)) {
  return(stats::p.adjust(p, method, n))
}

#' @export
p.adjust.matrix <- function(x, method="fdr") {
  if (method=="none") { return(x) }
  if (!all(x==t(x))) { stop("Matrix must be symmetric") }
  p.values.adjusted <- stats::p.adjust(x[lower.tri(x, diag=FALSE)], method = method)
  x[lower.tri(x, diag=FALSE)] <- p.values.adjusted
  x[upper.tri(x, diag=FALSE)] <- t(x)[upper.tri(x, diag=FALSE)]
  return(x)
}

#' @export
p.adjust.RSC.pmatrix <- function(x, method="fdr") {
  x <- as.RSC.pmatrix(p.adjust.matrix(x, method=method), p.adjust.method=method)
  return(x)
}

#' @export
p.adjust.RSC.cormatrix <- function(x, method="fdr") {
  x$p.value <- as.RSC.pmatrix(p.adjust.matrix(x$p.value, method=method), p.adjust.method=method)
  return(x)
}

#' @export
zero_percentile <- function(x) {
  zp <- ecdf(x)(0)
  return(zp)
}

#map_functions <- function(x, function_vector, parallel=FALSE, par.cores=4) {
#  
#}