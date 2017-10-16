# Regress out effects of whole brain volume, sex, and coil
#' @export
regress_out_confounds <- function(x, train_indices, scans) {
  df <- data.frame(absvol=x, study=scans$Study_Name, genotype=scans$Genotype, sex=scans$Mouse_Sex, coil=scans$Scan_Coil, background=scans$Background, wbv=vols_wbv, distcorr=scans$Distortion_Correction_File)
  
  df_train <- df[train_indices,]
  model_train <- lm(absvol ~ wbv * sex + coil, data=df_train)
  residuals_train <- df_train$absvol - predict(model_train, newdata=df_train)
  
  test_indices <- setdiff(1:dim(scans)[1], train_indices)
  df_test <- df[test_indices,]
  residuals_test <- df_test$absvol - predict(model_train, newdata=df_test)
  
  x_out <- numeric(dim(scans)[1])
  x_out[train_indices] <- residuals_train
  x_out[test_indices] <- residuals_test
  return(x_out)
}

#' @export
construct_raw_cormatrix_at_level <- function(tree, level, attribute="vols_reg") {
  this_level <- level
  level_strucs <- names(tree$Get(attribute = "name", filterFun = function(x) x$level == this_level))
  num_level_strucs <- length(level_strucs)
  cormatrix <- matrix(nrow=num_level_strucs, ncol=num_level_strucs, dimnames = list(level_strucs, level_strucs))
  for (i in 1:dim(cormatrix)[1]) {
    for (j in 1:dim(cormatrix)[2]) {
      struc_i <- rownames(cormatrix)[i]
      struc_j <- colnames(cormatrix)[j]
      cormatrix[i, j] <- cor(GetAttribute(FindNode(tree, struc_i), attribute = attribute),
                             GetAttribute(FindNode(tree, struc_j), attribute = attribute))
    }
    cat(".")
  }
  cat("Done.\n")
  return(cormatrix)
}

#' @export
construct_raw_cormatrices <- function(tree) {
  max_height <- GetAttribute(tree, "height")
  raw_level_cormatrices <- list()
  for (l in 1:max_height) {
    print(l)
    raw_level_cormatrices[[l]] <- construct_raw_cormatrix_at_level(tree, level = l)
  }
  return(raw_level_cormatrices)
}

#' @export
get_posterior_data <- function(tree, struc_i, struc_j, thres_cor=0.3, num_d=10, precision_per_level=100, p2r=p2r_table) {
  si_ancestors <- rev(names(FindNode(tree, struc_i)$Get('name', traversal = "ancestor")))
  sj_ancestors <- rev(names(FindNode(tree, struc_j)$Get('name', traversal = "ancestor")))
  max_height <- max(length(si_ancestors), length(sj_ancestors))
  
  init_cor <- cor(GetAttribute(FindNode(tree, si_ancestors[1]), "vols_reg"), GetAttribute(FindNode(tree, sj_ancestors[1]), "vols_reg"))
  updated_params <- fit_beta(init_cor, thres_cor, num_samples=num_d, precision=precision)
  
  for (iter in 2:max_height) {
    si_struc <- si_ancestors[min(iter, length(si_ancestors))]
    sj_struc <- sj_ancestors[min(iter, length(sj_ancestors))]
    #print(paste0("Working on: ", si_struc, " <> ", sj_struc))
    
    measured_cor <- cor(GetAttribute(FindNode(tree, si_struc), "vols_reg"), GetAttribute(FindNode(tree, sj_struc), "vols_reg"))
    updated_params <- update_beta_prior(a=updated_params[1], b=updated_params[2], 
                                        rho = measured_cor, threshold = thres_cor, num_samples = num_d, precision = precision_per_level)
  }
  initial_cor <- measured_cor
  map_proba <- (updated_params[1] - 1) / (updated_params[1] + updated_params[2] - 2)
  posterior_cor <- beta_to_rho(a=updated_params[1], b=updated_params[2], p2r_table=p2r, precision = precision_per_level*max_height, retval = "median")
  return(list(prob=map_proba, r_initial=initial_cor, r_posterior=posterior_cor))
}

# posterior distances
posterior_sqerr <- function()

#' @export
get_cormatrices <- function(tree, strucs=NULL, subjects=NULL, thres_cor=0.3, num_d=10, precision_per_level=100, p2r=p2r_table, ordering=NULL) {
  if (is.null(subjects)) {
    subjects <- 1:length(tree$vols_reg)
  }
  if (is.null(strucs)) {
    strucs <- names(tree$Get("name"))
  }
  treeCopy <- Clone(tree)
  treeCopy$Do(function(node) node$vols_reg <- node$vols_reg[subjects])
  
  pb <- txtProgressBar(min=0, max=((length(strucs))*(length(strucs)+1)/2), style = 3)
  k <- 1
  raw_cormtx <- posterior_cormtx <- posterior_proba <- matrix(nrow=length(strucs), ncol=length(strucs), dimnames = list(strucs, strucs)) 
  for (i in 1:length(strucs)) {
    for (j in 1:i) {
      posterior_data <- get_posterior_data(treeCopy, strucs[i], strucs[j], thres_cor = thres_cor, num_d=num_d, precision_per_level = precision_per_level, p2r = p2r)
      raw_cormtx[i, j] <- raw_cormtx[j, i] <- posterior_data$r_initial
      posterior_cormtx[i, j] <- posterior_cormtx[j, i] <- posterior_data$r_posterior
      posterior_proba[i, j] <- posterior_proba[j, i] <- posterior_data$prob
      setTxtProgressBar(pb, k)
      k <- k + 1
    }
  }
  posterior_proba[posterior_proba > 1] <- 1
  if (is.null(ordering)) {
    d <- dist(posterior_cormtx)
    fit <- hclust(d, method="average")
    o <- fit$order
  } else {
    o <- ordering
  }
  return(list(raw=raw_cormtx, posterior=posterior_cormtx, proba=posterior_proba, ordering=o))
}