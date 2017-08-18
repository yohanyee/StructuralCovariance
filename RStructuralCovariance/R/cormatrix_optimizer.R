# A method of successively normalizing both rows and columns of a matrix, a la Brad Efron and further described by Olshen et al.
# On convergence: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2868388/ (note that the proof contained within is flawed according to the 2013 paper published by the same authors)
# The type option takes the same arguments as in norm()
#' @export
successive_norm <- function(X, type="F", tol=1e-8, verbose=TRUE, na.set=0) {
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

# Given a named matrix, get a list of unique index (i,j) pairs to avoid recomputing symmetric elements
#' @export
get_unique_indexing <- function(mtx, diagonal_index=-1, set_repeats=NULL) {
  
  # Set row and column names, and determine which names are shared or unique between rows and columns
  rnames <- rownames(mtx)
  cnames <- colnames(mtx)
  different_row_columns <- c(setdiff(rnames, cnames), setdiff(cnames, rnames)) # Names not shared by columns and rows
  common_row_columns <- setdiff(c(rnames, cnames), different_row_columns) # Names shared by columns and rows
  num_common_row_columns <- length(common_row_columns)
  
  # Construct template index matrix
  index_mtx <- matrix(nrow=dim(mtx)[1], ncol=dim(mtx)[2], dimnames = list(rnames, cnames))
  
  # Start indexing 
  l <- 1
  
  # First, work on block of matrix that has common row and column names
  if (num_common_row_columns >= 1) {
    
    # In this block, set the diagonal to index specified
    if (is.null(diagonal_index)) {
      for (o in 1:num_common_row_columns) {
        index_mtx[common_row_columns[o], common_row_columns[o]] <- l
        l <- l + 1
      }
    } else {
      for (o in 1:num_common_row_columns) {
        index_mtx[common_row_columns[o], common_row_columns[o]] <- diagonal_index
      }
    }
    
    # Index the off-diagonals
    if (num_common_row_columns >= 2) {
      if (is.null(set_repeats)) {
        for (m in 2:num_common_row_columns) {
          for (n in 1:(m-1)) {
            index_mtx[common_row_columns[m], common_row_columns[n]] <- l
            index_mtx[common_row_columns[n], common_row_columns[m]] <- l
            l <- l + 1
          }
        }      
      } else {
        # Set repeated elements to NA
        for (m in 2:num_common_row_columns) {
          for (n in 1:(m-1)) {
            index_mtx[common_row_columns[m], common_row_columns[n]] <- l
            index_mtx[common_row_columns[n], common_row_columns[m]] <- set_repeats
            l <- l + 1
          }
        }
      }
    }
  }
  
  # Index the remaining elements
  num_remaining_elements <- length(index_mtx[is.na(index_mtx)])
  index_mtx[is.na(index_mtx)] <- l:(num_remaining_elements+l-1)
  
  # Provide indexing data in form of melted dataframe
  index_df <- reshape2::melt(index_mtx, value.name="index")
  index_df <- index_df[order(index_df$index),]
  index_df$unique <- !duplicated(index_df$index)
  index_df$unique[index_df$index <= 0] <- FALSE
  colnames(index_df)[1:2] <- c("rownames", "colnames")
  
  # Determine number of unique elements to be computed
  unique_elements = length(which(index_df$unique))

  # Return data
  return(list(lookup_table=index_mtx, indexing_dataframe=index_df, computable_elements=unique_elements))
}

# Given a named matrix, fill in symmetric pairs (including diagonals if necessary)
#' @export
fill_symmetric_elements <- function(mtx, diagonal_index=-1, neg_one_value=1) {
  
  # Get indexing information
  indexing <- get_unique_indexing(mtx, diagonal_index=diagonal_index)
  
  # Replace -1 index with replacement value if specified
  if (!is.null(neg_one_value)) {
    mtx[which(indexing$lookup_table==-1)] <- neg_one_value
  }
  
  # Determine remaining NA values left to replace
  num_na <- length(which(is.na(mtx)))
  
  # If NA values exist, then loop over all elements replacing NA values until either end of matrix reached or all NAs replace
  if ((num_na) >= 1) {
    
    # Start counting number of NAs replaced
    finished_na <- 0
    
    # Loop over matrix
    for (m in 1:dim(mtx)[1]) {
      for (n in 1:dim(mtx)[2]) {
        
        # If value is an NA, then start replacement procedure
        if (is.na(mtx[m, n])) {
          
          # Determine what index the NA is located at, and find other values at the same index
          l <- indexing$lookup_table[m,n]
          matched_elements <- mtx[indexing$lookup_table==l]
          matched_elements_index_filled <- which(!is.na(matched_elements))
          
          # Replace only if there is a non NA value at elements at a given index to replace with
          if (length(matched_elements_index_filled) >= 1) {
            
            # If multiple replacement values exist and are not the same, throw a warning
            if (length(matched_elements_index_filled) >= 2) {
              filled_elements <- matched_elements[matched_elements_index_filled]
              if (max(filled_elements) != min(filled_elements)) {
                warning(paste("Multiple values found to replace matrix element [", rownames(mtx)[m], ", ", colnames(mtx)[n], "]: ", paste(filled_elements, collapse=" "), ". Using first value.", sep=""))
              }
            }
            
            # Replace NA(s) with first legal replacement value, and increment counter with number of NAs replaced
            fill_value <- matched_elements[matched_elements_index_filled[1]]
            mtx[indexing$lookup_table==l] <- fill_value
            finished_na <- finished_na + length(which(is.na(matched_elements)))
          }
        }
        
        #If total number of initial NAs replaced, then return out of function
        if (finished_na==num_na) {
          return(mtx)
        }
      }
    }
  }
  
  # At the end of the loop, return out of function
  return(mtx)  
}

# Construct a template correlation matrix from column names
#' @export
construct_template_matrix <- function(strucs_source, strucs_target) {
  strucs_source <- as.character(strucs_source)
  strucs_target <- as.character(strucs_target)
  mtx <- matrix(NA, nrow=length(strucs_source), ncol=length(strucs_target), dimnames = list(strucs_source, strucs_target))
  return(mtx)
}

# Construct an empty matrix like another matrix
#' @export
construct_like_matrix <- function(like_mtx) {
  strucs_source <- as.character(rownames(like_mtx))
  strucs_target <- as.character(colnames(like_mtx))
  mtx <- matrix(NA, nrow=length(strucs_source), ncol=length(strucs_target), dimnames = list(strucs_source, strucs_target))
  return(mtx)
}

# Populate a matrix by calling a function at every element in an efficient manner
# Here, efficient means that only unique elements are computed (e.g., symmetric elements are not)
# function_string must be associated with a function that takes at least two named arguments: rowname, colname
#' @export
apply_on_matrix <- function(mtx, function_string, diagonal_index=-1, neg_one_value=1, indexing=NULL, ...) {
  
  # Throw a warning if there already are values in the matrix
  prefilled_values <- length(which(!is.na(mtx)))
  if (prefilled_values >= 1) {
    warning(paste("There already are", prefilled_values, "non NA values in the matrix"))
  }
  
  # Get indexing information
  if (is.null(indexing)) {
    indexing <- get_unique_indexing(mtx, diagonal_index=diagonal_index)
  }
  index_df <- subset(indexing$indexing_dataframe, indexing$indexing_dataframe$index != -1)
  index_df_compute <- subset(index_df, index_df$unique==TRUE)
  index_df_copy <- subset(index_df, index_df$unique==FALSE)
  
  # Replace -1 index with replacement value if specified
  if (!is.null(neg_one_value)) {
    mtx[which(indexing$lookup_table==-1)] <- neg_one_value
  }
  
  # Populate matrix at unique elements
  for (li in 1:dim(index_df_compute)[1]) {
    rowname <- as.character(index_df_compute$rownames[li])
    colname <- as.character(index_df_compute$colnames[li])
    args <- list(rowname=rowname, colname=colname, ...)
    mtx[rowname, colname] <- do.call(function_string, args)
  }
  
  # Copy non-unique elements (if they exist)
  if (dim(index_df_copy)[1] >=1 ) {
    for (lc in 1:dim(index_df_copy)[1]) {
      
      # Get row and column names of element at which to copy
      l <- index_df_copy$index[lc]
      rowcopy <- as.character(index_df_copy$rownames[lc])
      colcopy <- as.character(index_df_copy$colnames[lc])
      
      # Get row and column names of element to copy
      li <- which(index_df_compute$index==l)
      rowname <- as.character(index_df_compute$rownames[li])
      colname <- as.character(index_df_compute$colnames[li])
      
      # Copy
      mtx[rowcopy, colcopy] <- mtx[rowname, colname]
    }
  }
  # Return computed matrix
  return(mtx)
}

# Construct a matrix by mapping a function over each element
#' @export
construct_matrix <- function(df, strucs_source, strucs_target, function_string, ...) {
  mtx <- construct_template_matrix(strucs_source, strucs_target) 
  mtx <- apply_on_matrix(mtx, function_string, df=df, ...)
  return(mtx)
}

# Conventional cor function adapted for efficient correlation matrix computation
#' @export
corr <- function(df, rowname, colname) {
  return(cor(df[, rowname], df[, colname]))
}

# Dot product adapted for efficient pairwise matrix computation
#' @export
dot <- function(df, rowname, colname) {
  return(sum(df[,rowname]*df[,colname]))
}

# Initialize correlation pre-computation terms for matrix optimization
# These are objects that hold the precomputed terms to help with efficient correlation computation
# Object values should be updated with every iteration (as optimal row is dropped)
# Source and target terms: Matrix (2 x {s,t}) with first row containing sum, second is sum of squared term
# Interaction terms: array (s x t) of multiplication terms between source and target values
# n, integer of remaining rows
# dropped: vector of integers that show the dropped rows in order
#' @export
cormatrix_optimizer_precompute_terms <- function(df, strucs_source, strucs_target, indexing_for_interactions=NULL) {
  T_source <- df[,strucs_source]
  T_target <- df[,strucs_target]
  T_source_terms <- construct_template_matrix(c("sum", "sumsq"), colnames(T_source))
  T_source_terms[1,] <- colSums(T_source)
  T_source_terms[2,] <- colSums(T_source^2)
  T_target_terms <- construct_template_matrix(c("sum", "sumsq"), colnames(T_target))
  T_target_terms[1,] <- colSums(T_target)
  T_target_terms[2,] <- colSums(T_target^2)
  if (is.null(indexing_for_interactions)) {
    T_interaction_terms <- construct_matrix(df, strucs_source=strucs_source, strucs_target=strucs_target, function_string="dot", diagonal_index=NULL) 
  } else {
    T_interaction_terms <- construct_matrix(df, strucs_source=strucs_source, strucs_target=strucs_target, function_string="dot", indexing=indexing_for_interactions) 
  }
  out <- list(df=df, strucs_source=strucs_source, strucs_target=strucs_target, 
              source=T_source_terms, target=T_target_terms, interaction=T_interaction_terms, 
              n=dim(df)[1], n_original=dim(df)[1], dropped_rows=c())
  return(out)
}

# Optimized Pearson correlation that uses precomputed terms to efficiently calculate the correlation coefficient after dropping a row
#' @export
cormatrix_optimizer_drop_row_compute_element <- function(terms, drop_row, rowname, colname) {
  sdrop <- terms$df[drop_row,rowname]
  tdrop <- terms$df[drop_row,colname]
  sterms <- terms$source[,rowname]
  s <- as.numeric(sterms[1])
  s2 <- as.numeric(sterms[2])
  tterms <- terms$target[,colname]
  t <- as.numeric(tterms[1])
  t2 <- as.numeric(tterms[2])
  n <- terms$n
  st <- as.numeric(terms$interaction[rowname, colname])
  r <- (((n-1)*(st - (sdrop*tdrop)))-((s-sdrop)*(t-tdrop)))/(sqrt(((n-1)*(s2-(sdrop)^2))-(s-sdrop)^2)*sqrt(((n-1)*(t2-(tdrop)^2))-(t-tdrop)^2))
  return(r)
}

# Optimized Pearson correlation that uses precomputed terms to efficiently calculate the correlation coefficient 
#' @export
cormatrix_optimizer_compute_element <- function(terms, rowname, colname) {
  sterms <- terms$source[,rowname]
  s <- as.numeric(sterms[1])
  s2 <- as.numeric(sterms[2])
  tterms <- terms$target[,colname]
  t <- as.numeric(tterms[1])
  t2 <- as.numeric(tterms[2])
  n <- terms$n
  st <- as.numeric(terms$interaction[rowname, colname])
  r <- (n*st-s*t)/(sqrt(n*s2-s^2)*sqrt(n*t2-t^2))
  return(r)
}

# Efficiently compute the correlation matrix using precomputed terms
#' @export
cormatrix_optimizer_compute_matrix <- function(terms, indexing_for_computation=NULL) {
  mtx <- construct_template_matrix(terms$strucs_source, terms$strucs_target)
  mtx <- apply_on_matrix(mtx, "cormatrix_optimizer_compute_element", terms=terms, indexing=indexing_for_computation)
  return(mtx)
}

# Efficiently compute the correlation matrix after dropping a row using precomputed terms
#' @export
cormatrix_optimizer_drop_row_and_compute_matrix <- function(terms, drop_row, indexing_for_computation=NULL) {
  mtx <- construct_template_matrix(terms$strucs_source, terms$strucs_target)
  mtx <- apply_on_matrix(mtx, "cormatrix_optimizer_drop_row_compute_element", terms=terms, drop_row=drop_row, indexing=indexing_for_computation)
  return(mtx)
}

# Compare (correlate) two matrices by vectorizing unique elements
#' @export
cormatrix_optimizer_correlate <- function(mtx, basemtx, indexing_for_comparison=NULL) {
  if (is.null(indexing_for_comparison)) {
    indexing <- get_unique_indexing(mtx, diagonal_index = -1, set_repeats = -1)
  } else {
    indexing <- indexing_for_comparison
  }
  indices <- which(indexing$lookup_table >= 1)
  r <- cor(mtx[indices], basemtx[indices])
  return(r)
}

# Compare (correlate) two matrices by vectorizing unique elements, with one of the matrices efficiently computed after dropping a row
#' @export
cormatrix_optimizer_drop_row_and_correlate <- function(terms, drop_row, basemtx, indexing_for_computation=NULL, indexing_for_comparison=NULL) {
  mtx <- cormatrix_optimizer_drop_row_and_compute_matrix(terms, drop_row, indexing_for_computation)
  r <- cormatrix_optimizer_correlate(mtx, basemtx, indexing_for_comparison=indexing_for_comparison)
  return(r)
}

# Update precomputed terms after dropping row(s)
#' @export
cormatrix_optimizer_update_terms <- function(terms, drop_rows, indexing_for_interactions=NULL) {
  strucs_source <- terms$strucs_source
  strucs_target <- terms$strucs_target
  for (dr in as.character(drop_rows)) {
    sdrop <- terms$df[dr, strucs_source]
    tdrop <- terms$df[dr, strucs_target]
    
    # Change source values
    terms$source[1,] <- terms$source[1,] - sdrop
    terms$source[2,] <- terms$source[2,] - (sdrop^2)
    
    # Change target values
    terms$target[1,] <- terms$target[1,] - tdrop
    terms$target[2,] <- terms$target[2,] - (tdrop^2)
    
    # Change interaction values
    if (is.null(indexing_for_interactions)) {
      dr_interaction <- construct_matrix(df=t(terms$df[dr,]), strucs_source = strucs_source, strucs_target=strucs_target, function_string = "dot", diagonal_index=NULL)
    } else {
      dr_interaction <- construct_matrix(df=t(terms$df[dr,]), strucs_source = strucs_source, strucs_target=strucs_target, function_string = "dot", indexing=indexing_for_interactions)
    }
    terms$interaction <- terms$interaction - dr_interaction
    
    # Change n
    terms$n <- terms$n - 1
    
    # Add to dropped terms
    terms$dropped_rows <- c(terms$dropped_rows, dr)
  }
  return(terms)
}

# Determine which row(s) to drop
#' @export
cormatrix_optimizer_passthru <- function(terms, basemtx,
                                         batch_size=1, cor_objective=1, 
                                         probabilistic=FALSE, probabilistic_weight_exponent=2, 
                                         parallel=4, cluster_initialize=TRUE, cluster_shutdown=TRUE,
                                         indexing_for_computation=NULL, indexing_for_comparison=NULL, indexing_for_interactions=NULL) {
  avail_rows <- rownames(terms$df)[!(rownames(terms$df) %in% terms$dropped_rows)]
  num_rows <- length(avail_rows)
  objective_df <- data.frame(row=character(num_rows), cor=numeric(num_rows), dist_to_objective=numeric(num_rows))
  objective_df$row <- as.character(objective_df$row)
  if (parallel==1) {
    i <- 1
    for (ar in avail_rows) {
      r <- cormatrix_optimizer_drop_row_and_correlate(terms, drop_row = ar, basemtx = basemtx, indexing_for_computation=indexing_for_computation, indexing_for_comparison=indexing_for_comparison)
      objective_df$row[i] <- as.character(ar)
      objective_df$cor[i] <- r
      objective_df$dist_to_objective[i] <- abs(r - cor_objective)
      i <- i + 1
      cat(".")
    }
    cat("Done.\n")
  } else {
    parallel <- min(parallel, avail_rows)
    if (cluster_initialize) {
      library(doParallel)
      cl <- makeCluster(parallel, outfile="", homogeneous=TRUE, rscript_args="--vanilla")
      registerDoParallel(cl) 
    }
    iterations <- length(avail_rows)
    prog <- txtProgressBar(min=1, max=iterations, style=3)
    export_functions <- c("cormatrix_optimizer_drop_row_compute_element", "cormatrix_optimizer_drop_row_and_correlate", 
                          "cormatrix_optimizer_drop_row_and_compute_matrix", "cormatrix_optimizer_correlate",
                          "construct_template_matrix", "apply_on_matrix", "cormatrix_optimizer_compute_element", "get_unique_indexing")
    export_variables <- c("avail_rows", "terms", "basemtx", "cor_objective", "indexing_for_computation", "indexing_for_comparison", "indexing_for_interactions")
    objective_df <- foreach(i=1:iterations, .export = c(export_functions), .combine=rbind, .verbose = FALSE) %dopar% {
      ar <- avail_rows[i]
      r <- cormatrix_optimizer_drop_row_and_correlate(terms, drop_row = ar, basemtx = basemtx, indexing_for_computation=indexing_for_computation, indexing_for_comparison=indexing_for_comparison)
      setTxtProgressBar(prog, i) 
      data.frame(row=as.character(ar), cor=r, dist_to_objective=abs(r - cor_objective))
    }
    close(prog)
    if (cluster_shutdown) {
      registerDoSEQ()
    }
  }

  # Determine which rows to drop
  if (probabilistic==TRUE) {
    nearness <- 1/((objective_df$dist_to_objective)^(probabilistic_weight_exponent))
    weights <- nearness/sum(nearness)
    rows_to_drop <- objective_df$row[sample(1:length(weights), size = batch_size, replace = FALSE, prob = weights)]
  } else {
    rows_to_drop <- objective_df$row[order(objective_df$dist_to_objective)[1:batch_size]]
  }
  rows_to_drop <- as.character(rows_to_drop)
  
  # Update terms
  newterms <- cormatrix_optimizer_update_terms(terms = terms, drop_rows = rows_to_drop, indexing_for_interactions=indexing_for_interactions)
  
  # Compute new correlation with basemtx and return that
  newmtx <- cormatrix_optimizer_compute_matrix(terms = newterms, indexing_for_computation=indexing_for_computation)
  r <- cormatrix_optimizer_correlate(newmtx, basemtx, indexing_for_comparison=indexing_for_comparison)
  computation_info <- list(probabilistic=probabilistic, probabilistic_weight_exponent=probabilistic_weight_exponent, parallel=parallel)
  out <- list(terms=newterms, r=r, dropped_rows=rows_to_drop, computation_info=computation_info)
  return(out)
}

# Repeatedly drop rows to determine an enriched set of rows that optimize the correlation between two matrices
#' @export
cormatrix_optimizer_optimize <- function(X, Y, strucs_source, strucs_target, batch_sizes=NULL, batch_definitions=NULL, precompute_indexing=TRUE, cor_objective=1, min_rows=3, tol=1e-6, max_passes=-1, logfile="optimization_output.log", outfile="optimization_output.RData", baserowfile="optimization_base.txt", enrichedrowfile="optimization_enriched.txt", is_Allen_gene=FALSE, ...) {
  
  cat("\n##################\n")
  cat(paste("# INITIALIZATION #\n"))
  cat("##################\n\n")
  
  # Start timing
  start.time <- Sys.time()
  cat("\n")
  cat(paste("* Initializing optimization process\n"))
  cat(paste("* Number of prunable rows:", dim(X)[1], "\n"))
  cat(paste("* Starting time:", start.time, "\n"))
  cat("\n")
  
  # Precompute base matrix (to which the X data's correlation matrix is optimized, the indexing terms, and the correlation terms)
  cat("\n")
  cat(paste("* Computing base matrix to optimize to\n"))
  cat("\n")
  basemtx <- construct_matrix(Y, strucs_source, strucs_target, function_string = "corr")
  if (precompute_indexing) {
    cat("\n")
    cat(paste("* Precomputing matrix indices\n"))
    cat("\n")
    indexing_for_computation <- get_unique_indexing(basemtx, diagonal_index = -1, set_repeats = NULL)
    indexing_for_comparison <- get_unique_indexing(basemtx, diagonal_index = -1, set_repeats = -1)
    indexing_for_interactions <- get_unique_indexing(basemtx, diagonal_index = NULL, set_repeats = NULL)
  } else {
    indexing_for_computation <- NULL
    indexing_for_comparison <- NULL
    indexing_for_interactions <- NULL
  }
  cat("\n")
  cat(paste("* Precomputing correlation terms\n"))
  cat("\n")
  terms <- cormatrix_optimizer_precompute_terms(X, strucs_source=strucs_source, strucs_target=strucs_target, indexing_for_interactions=indexing_for_interactions)
  input_terms <- terms
  
  # Compute base data comparison before dropping rows
  cat("\n")
  cat(paste("* Computing starting matrix\n"))
  startmtx <- cormatrix_optimizer_compute_matrix(terms=terms, indexing_for_computation=indexing_for_computation)
  r <- cormatrix_optimizer_correlate(startmtx, basemtx, indexing_for_comparison=indexing_for_comparison)
  cat(paste("* Initial r-value:", r, "\n"))
  cat("\n")
  
  # Base gene set
  cat("\n")
  cat(paste("* Determining base gene set\n"))
  cat("\n")
  base_set <- rownames(terms$df)
  
  # Collect initial data
  cat("\n")
  cat(paste("* Logging initial data\n"))
  cat("\n")
  current.time <- Sys.time()
  optimizer_df <- data.frame(pass=numeric(terms$n_original+1), num_rows_left=numeric(terms$n_original+1), r=numeric(terms$n_original+1), dropped_rows=character(terms$n_original+1), current_time=character(terms$n_original+1), elapsed_time=character(terms$n_original+1))
  optimizer_df$pass[1] <- 0
  optimizer_df$num_rows_left[1] <- terms$n
  optimizer_df$r[1] <- r
  optimizer_df$dropped_rows[1] <- ""
  optimizer_df$dropped_rows <- as.character(optimizer_df$dropped_rows)
  optimizer_df$current_time[1] <- current.time
  optimizer_df$current_time <- as.character(optimizer_df$current_time)
  optimizer_df$elapsed_time[1] <- current.time - start.time
  optimizer_df$elapsed_time <- as.character(optimizer_df$elapsed_time)
  
  # Log data if required
  if (!is.null(logfile)) {
    cat("\n")
    cat(paste("* Writing initial data to file:", logfile, "\n"))
    cat("\n")
    logrow <- optimizer_df[1,]
    if (file.exists(logfile)) {
      logcopy <- paste(logfile, "backup", paste(sample(c(letters, LETTERS, 0:9), 10), collapse = ""), sep="_")
      file.copy(logfile, logcopy)
      warning(paste("File:", logfile, "exists! Saving it as a backup to", logcopy))
    }
    write.table(logrow, file = logfile, append = FALSE, quote=4, sep = ",", row.names = FALSE, col.names = TRUE)
  }
  
  # Begin optimization
  cat("\n################\n")
  cat(paste("# OPTIMIZATION #\n"))
  cat("################\n\n")
  terminated <- FALSE
  objective_reached <- FALSE
  termination_reason <- NA
  pass <- 1
  while ((terms$n > min_rows) & (abs(r-cor_objective)>=tol) & ((pass <= max_passes) | max_passes < 0)) {
    # Print
    current.time <- Sys.time()
    cat("\n")
    cat(paste("* Iteration:", pass, "of", num_passes, "\n"))
    cat(paste("* Number of rows remaining:", terms$n, "\n"))
    cat(paste("* R-value:", r, "\n"))
    cat(paste("* Current time:", current.time, "\n"))
    cat(paste("* Time elapsed:", current.time - start.time, "\n"))
    cat(paste("* Estimated time of completion:", , "\n"))
    cat(paste("* Progress:", , "\n"))
    cat("\n")
    
    # Set batch size
    if (!is.null(batch_definitions)) {
      if (!is.null(batch_sizes)) {
        warning("Only one of batch_definitions and batch_sizes should be set. Using batch_definitions")
      } 
      this_batch_size <- batch_definitions$batch_size[which(batch_definitions$below > terms$n & batch_definitions$aboveeq <= terms$n)]
      if (length(this_batch_size) > 1) {
        this_batch_size <- this_batch_size[1]
        warning(paste("Multiple batch sizes fit this definition. Using batch size:", this_batch_size))
      }
    } else {
      if (!is.null(batch_sizes)) {
        this_batch_size <- ifelse(is.na(batch_sizes[pass]), batch_sizes[length(batch_sizes)], batch_sizes[pass])
      } else {
        this_batch_size <- 1
      }
    }
    
    # Return if batch size will result in going below min_rows
    if ((terms$n - this_batch_size) < min_rows) {
      terminated <- TRUE
      termination_reason <- "ROW_LIMIT_REACHED"
      break
    }
    
    # Optimize
    if (pass==1) {
      cluster_start <- TRUE
      cluster_stop <- FALSE
    } else {
      cluster_start <- FALSE
      cluster_stop <- FALSE
    }
    #optimize <- cormatrix_optimizer_passthru(terms, basemtx, batch_size=this_batch_size, cor_objective=cor_objective,
    #                                        indexing_for_computation=indexing_for_computation, indexing_for_comparison=indexing_for_comparison, indexing_for_interactions=indexing_for_interactions, 
    #                                        cluster_initialize=cluster_start, cluster_shutdown=cluster_stop, ...)
    optimize <- cormatrix_optimizer_passthru(terms, basemtx, batch_size=this_batch_size, cor_objective=cor_objective, indexing_for_computation=indexing_for_computation, indexing_for_comparison=indexing_for_comparison, indexing_for_interactions=indexing_for_interactions, parallel=6, cluster_initialize = cluster_start, cluster_shutdown = cluster_stop)
    terms <- optimize$terms
    
    # Collect data
    current.time <- Sys.time()
    optimizer_df$pass[pass+1] <- pass
    optimizer_df$num_rows_left[pass+1] <- terms$n
    optimizer_df$r[pass+1] <- optimize$r
    optimizer_df$dropped_rows[pass+1] <- paste(optimize$dropped_rows, collapse=";")
    optimizer_df$current_time[pass+1] <- current.time
    optimizer_df$elapsed_time[pass+1] <- current.time - start.time
    
    # Log data if required
    if (!is.null(logfile)) {
      logrow <- optimizer_df[pass+1,]
      write.table(logrow, file = logfile, append = TRUE, quote=4, sep = ",", row.names = FALSE, col.names = FALSE)
    }
    
    # Increase pass
    pass <- pass + 1
    r <- optimize$r
  }
  if (optimize$computation_info$parallel > 1) {
    registerDoSEQ()
  }
  
  # Record termination
  terminated <- TRUE
  if (terms$n <= min_rows) {
    termination_reason <- "ROW_LIMIT_REACHED"
  } else if (abs(r-cor_objective)<tol) {
    termination_reason <- "OBJECTIVE_REACHED"
    objective_reached <- TRUE
  } else if (pass > max_passes) {
    termination_reason <- "ITERATION_LIMIT_REACHED"
  } else {
    termination_reason <- "UNKNOWN"
  }
  cat("\n")
  cat(paste("* Optimization terminated\n"))
  cat(paste("* Termination due to:", termination_reason, "\n"))
  cat(paste("* Cor objective reached:", objective_reached, "\n"))
  cat("\n")

  cat("\n################\n")
  cat(paste("# FINISHING UP #\n"))
  cat("################\n\n")
    
  # Compute peak enrichment data
  cat("\n")
  cat(paste("* Computing peak enrichment data and matrix\n"))
  cat("\n")
  peak_row <- which.min(abs(optimizer_df$r-cor_objective))
  peak_r <- optimizer_df$r[peak_row]
  removed_set <- unlist(strsplit(optimizer_df$dropped_rows[2:peak_row], ";"))
  peak_enriched_set <- setdiff(rownames(terms$df), removed_set)
  peakmtx <- construct_matrix(terms$df[(which(rownames(terms$df) %in% peak_enriched_set)),], strucs_source, strucs_target, function_string = "corr")
  
  # If Allen Institute data, then relabel rownames to reflect unique genes
  if (is_Allen_gene) {
    cat("\n")
    cat(paste("* Converting Allen Institute rownames to gene lists\n"))
    cat("\n")
    base_set <- unique(sapply(strsplit(base_set, "_sid"), "[[", 1))
    peak_enriched_set <- unique(sapply(strsplit(peak_enriched_set, "_sid"), "[[", 1))
  }
  
  # End timing
  end.time <- Sys.time()
  
  # Construct output object
  cat("\n")
  cat(paste("* Constructing output\n"))
  cat("\n")
  out <- list(inputs=list(X=X, Y=Y, strucs_source=strucs_source, strucs_target=strucs_target, 
                          batch_sizes=batch_sizes, batch_definitions=batch_definitions, 
                          precompute_indexing=precompute_indexing, cor_objective=cor_objective, 
                          probabilistic=optimize$computation_info$probabilistic, probabilistic_weight_exponent=optimize$computation_info$probabilistic_weight_exponent, 
                          parallel=optimize$computation_info$parallel, 
                          min_rows=min_rows, tol=tol, max_passes=max_passes, logfile=logfile, outfile=outfile),
              terms=list(start=input_terms, final=terms),
              matrices=list(base=basemtx, start=startmtx, peak=peakmtx),
              optimization=list(data=optimizer_df[1:pass,],  passes=(pass-1), terminated=terminated, termination_reason=termination_reason, objective_reached=objective_reached),
              enrichment=list(base_set=base_set, peak_r=peak_r, peak_enriched_set=peak_enriched_set),
              debug=list(sysinfo=as.list(Sys.info()), timing=list(start=start.time, end=end.time, walltime=(end.time-start.time)), packages=search())
              )
  
  # Save object if desired
  if (!is.null(outfile)) {
    cat("\n")
    cat(paste("* Saving output to:", outfile,"\n"))
    cat("\n")
    save("out", file = outfile)
  }
  
  if (!is.null(baserowfile)) {
    cat("\n")
    cat(paste("* Saving base set to:", baserowfile,"\n"))
    cat("\n")
    write.table(base_set, file = baserowfile, row.names = FALSE, col.names = FALSE, quote=FALSE)
  }
  
  if (!is.null(enrichedrowfile)) {
    cat("\n")
    cat(paste("* Saving enriched set to:", enrichedrowfile,"\n"))
    cat("\n")
    write.table(peak_enriched_set, file = enrichedrowfile, row.names = FALSE, col.names = FALSE, quote=FALSE)
  }
  
  cat("\n########\n")
  cat(paste("# DONE #\n"))
  cat("########\n\n")
  
  #Return
  return(out)
}