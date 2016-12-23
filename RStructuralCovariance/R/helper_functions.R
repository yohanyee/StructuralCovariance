get_volume_info <- function(v) {
  subjects <- list(num=dim(v)[1], names=rownames(v))
  structures <- list(num=dim(v)[2], names=colnames(v))
  out <- list(subjects=subjects, structures=structures)
  return(out)
}

resample_volumes <- function(v, size=NULL, replace=TRUE, prob=NULL) {
  vinfo <- get_volume_info(v)
  size <- (if (is.null(size)) vinfo$subjects$num else size)
  size <- (if (size > vinfo$subjects$num | size < 3) vinfo$subjects$num else size)
  indices <- sample(1:size, size=size, replace=replace, prob=prob)
  v <- v[indices,]
  return(v)
}

normalize_volumes <- function(v) {
  return(v/rowSums(v))
}

correlate_volumes <- function(v, method="pearson", quiet=FALSE) {
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
      ct <- cor.test(v[,i], v[,j], cor.method=method)
      cmat[i,j] <- cmat[j,i] <- ct$estimate
      pmat[i,j] <- pmat[j,i] <- ct$p.value
    }
  }
  cat("Done.\n")
  
  # Coerce to RSC matrix objects
  cmat <- as.RSC.cmatrix(cmat)
  pmat <- as.RSC.pmatrix(pmat)
  
  cormatrix <- as.RSC.cormatrix(list(correlation=cmat, p.value=pmat))
  return(out)
}

p.adjust.matrix <- function(x, method="fdr") {
  if (method=="none") { return(x) }
  if (!all(x==t(x))) { stop("Matrix must be symmetric") }
  p.values.adjusted <- p.adjust(x[lower.tri(x, diag=FALSE)], method = method)
  x[lower.tri(x, diag=FALSE)] <- p.values.adjusted
  x[upper.tri(x, diag=FALSE)] <- t(x)[upper.tri(x, diag=FALSE)]
  return(x)
}

p.adjust.RSC.pmatrix <- function(x, method="fdr") {
  x <- as.RSC.pmatrix(p.adjust.matrix(x, method=method), p.adjust.method=method)
  return(x)
}

p.adjust.RSC.cormatrix <- function(x, method="fdr") {
  x$p.value <- as.RSC.pmatrix(p.adjust.matrix(x$p.value, method=method), p.adjust.method=method)
  return(x)
}