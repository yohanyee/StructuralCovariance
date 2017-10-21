library(data.tree)

#' @export
hanatFromAnatMatrix <- function(anatMatrix, norm.function=NULL, hclust.method="complete", verbose=TRUE, ...) {
  # Set up constants and association matrix
  num_strucs <- dim(anatMatrix)[2]
  assoc_mtx <- matrix(nrow=num_strucs, ncol=num_strucs, dimnames = list(colnames(anatMatrix), colnames(anatMatrix)))
  diag(assoc_mtx) <- 1
  
  # Pre-normalize volumes if required
  if (is.null(norm.function)) {
    normAnatMatrix <- anatMatrix
  } else {
    normAnatMatrix <- do.call(norm.function, args = list(X=anatMatrix, ...))
  }
  
  # Set up text progress bar
  if (verbose) {
    cat("Computing association matrix.\n")
    pb <- txtProgressBar(max=(num_strucs*(num_strucs-1)/2), style=3)
    prog <- 1
  }
  
  # Compute association matrix
  for (i in 2:num_strucs) {
    for (j in 1:(i-1)) {
      if (is.null(norm.function)) {
        norm_vol_union <- normAnatMatrix[,i] + normAnatMatrix[,j]
      } else {
        norm_vol_union <- do.call(norm.function, args = list(X=(anatMatrix[,i] + anatMatrix[,j]), ...))
      }
      assoc_mtx[i, j] <- assoc_mtx[j, i] <- min(cor(norm_vol_union, normAnatMatrix[,i]), cor(norm_vol_union, normAnatMatrix[,j]))
      if (verbose) {
        setTxtProgressBar(pb, prog)
        prog <- prog + 1
      }
    }
  }
  
  if (verbose) {
    close(pb)
    cat("Done.\n")
  }
  
  # Compute dendrogram and tree
  anat_dendro <- as.dendrogram(hclust(dist(assoc_mtx), method=hclust.method))
  hanat_tree <- as.Node(anat_dendro, name = "brain")
  
  # Add volumes to hierarchy
  # RMINC::addVolumestoHierarchy() uses anatIDs attribute of input volumes. 
  # Since anatCombined objects don't have this, I've adapted that code to work here by matching column names instead
  
  # First, add volumes to leaf nodes
  hanat_tree$Do(function(x) {
    if (isLeaf(x)) {
      x$volumes <- anatMatrix[,x$name]
      x$meanVolume <- mean(x$volumes)
    }
  })
  
  # Then, traverse tree and aggregate values
  hanat_tree$Do(function(x) {
    x$volumes <- Aggregate(x, "volumes", rowSums)
  }, traversal = "post-order", filterFun = isNotLeaf
  )
  
  hanat_tree$Do(function(x) {
    x$meanVolume <- Aggregate(x, "meanVolume", sum)
  }, traversal = "post-order"
  )
  
  # Lastly, if norm.function is set, then calculate normalized volumes
  if (!is.null(norm.function)) {
    
    # Set up text progress bar
    if (verbose) {
      cat("Adding normalized volumes to nodes.\n")
      pb <- txtProgressBar(max=hanat_tree$totalCount, style=3)
      prog <- 1
    }
    
    hanat_tree$Do(function(x) {
      x$normVolumes <- normAnatMatrix[,x$name]
      if (verbose) {
        setTxtProgressBar(pb, prog)
        prog <- prog + 1
      }
    }, filterFun = isLeaf
    )
    
    hanat_tree$Do(function(x) {
      x$normVolumes <- do.call(norm.function, args = list(X=x$volumes, ...))
      if (verbose) {
        setTxtProgressBar(pb, prog)
        prog <- prog + 1
      }
    }, filterFun = isNotLeaf
    )
    
    if (verbose) {
      close(pb)
      cat("Done.\n")
    }
  }
  
  return(hanat_tree)
}