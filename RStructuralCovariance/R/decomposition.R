#o <- hclust(dist(cov2cor(X)), method="average")$order

#corrplot(cov2cor(X)[o,o], method="color", tl.col="black", tl.cex=0.2, col=col2rev(100))


pca_reconstruct <- function(X, n=5) {
  mean_vec <- colMeans(X)
  pca <- prcomp(X, retx = TRUE, scale=FALSE)
  X_recon <- scale(pca$x[,1:n] %*% t(pca$rotation[,1:n]), center=-mean_vec, scale=FALSE)
  return(as.matrix(Matrix::nearPD(X_recon)$mat))
}

#X_recon <- cov2cor(pca_reconstruct(X, n=20))
#X_recon[X_recon < -1] <- -1
#X_recon[X_recon > 1] <- 1

#corrplot(cov2cor(X)[o,o], method="color", tl.col="black", tl.cex=0.2, col=col2rev(100))
#corrplot(cov2cor(X_recon)[o,o], method="color", tl.col="black", tl.cex=0.2, col=col2rev(100))
#corrplot(X_recon[o,o], method="color", tl.col="black", tl.cex=0.2, col=col2rev(100))


#mvrnorm(n=10, mu = colMeans(X), Sigma=lowrank_cov$mat)

#X_cor <- X_recon/max(abs(X_recon))
#pca_reconstruct(X, n=20)
