library(data.tree)

#' Construct and populate an anatomical hierarchy from volume data
#' 
#' Given volume data as a rectangular array, construct an anatomical hierarchy by grouping structures that share volume patterns.
#' More specifically, a pairwise association matrix is constructed for each pair of structures, where the association is defined as the minimum of the correlation between each of the two structures and the sum of both structures.
#' In other words, structures are grouped in a way that parent structures (with volume the sum of its children) correlate the most with its children structure. 
#' The association matrix is then hierarchically clustered, and a data.tree is constructed from the resulting dendrogram. 
#' Volumes can (and should) be normalized when computing associations.
#'  
#' @param anatMatrix A rectangular array of volume data, with structure names as column names.
#' @param norm.function A function to normalize the volume data.
#' @param train.indices an optional set of indices corresponding to the rows of \code{anatMatrix} indicating the data on which 1) the anatomical hierarchy will be determined, and 2) if \code{norm.function} is set to \code{normByModel}, the data that should be used in the model fitting. If \code{NULL}, then all rows are used in the hierarchy construction.
#' @param hclust.method Method to be passed to \code{hclust()}.
#' @param verbose Be verbose.
#' @param ... Further arguments to be passed to the function given by \code{norm.function}.
#' @return An anatomical hierarchy with volume attributes (\code{volumes}, \code{meanVolume}, and if applicable, \code{normVolumes})
#' 
#' @export
hanatFromAnatMatrix <- function(anatMatrix, norm.function=NULL, train.indices=NULL, hclust.method="complete", verbose=TRUE, ...) {
  # Set up constants and association matrix
  num_strucs <- dim(anatMatrix)[2]
  assoc_mtx <- matrix(nrow=num_strucs, ncol=num_strucs, dimnames = list(colnames(anatMatrix), colnames(anatMatrix)))
  diag(assoc_mtx) <- 1
  if (is.null(train.indices)) {
    train.indices <- all.indices
  }
  if (!all(train.indices %in% all.indices)) {
    stop("train_indices are not valid indices")
  }
  
  # Pre-normalize volumes if required
  if (is.null(norm.function)) {
    normAnatMatrix <- anatMatrix
  } else {
    if (quote(norm.function)=="normByModel") {
      normAnatMatrix <- do.call(norm.function, args = list(X=anatMatrix, train.indices=train.indices, ...))
    } else {
      normAnatMatrix <- do.call(norm.function, args = list(X=anatMatrix, ...))
    }
    
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
        if (quote(norm.function)=="normByModel") {
          norm_vol_union <- do.call(norm.function, args = list(X=(anatMatrix[,i] + anatMatrix[,j]), train.indices=train.indices, ...))
        } else {
          norm_vol_union <- do.call(norm.function, args = list(X=(anatMatrix[,i] + anatMatrix[,j]), ...))
        }
      }
      assoc_mtx[i, j] <- assoc_mtx[j, i] <- min(cor(norm_vol_union[train.indices], normAnatMatrix[train.indices,i]), cor(norm_vol_union[train.indices], normAnatMatrix[train.indices,j]))
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
      if (quote(norm.function)=="normByModel") {
        x$normVolumes <- do.call(norm.function, args = list(X=x$volumes, train.indices=train.indices, ...))
      } else {
        x$normVolumes <- do.call(norm.function, args = list(X=x$volumes, ...))
      }
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
  
  # Name all nodes with (somewhat) more meaningful names
  hanat_tree$Do(function(x) {
    child_name <- names(which.max(lapply(x$children, "[[", "meanVolume")))
    if (startsWith(child_name, "Level ")) {
      new_name <- gsub(paste0("Level ", x$level + 1, ": "), paste0("Level ", x$level, ": "), child_name)
    } else {
      new_name <- paste0("Level ", x$level, ": ", child_name)
    }
    x$name <- new_name
  }, filterFun = isNotLeaf, traversal = "post-order"
  )
  
  return(hanat_tree)
}