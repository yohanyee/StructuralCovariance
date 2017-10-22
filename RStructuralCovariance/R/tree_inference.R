# TODO: test on two structures that share the same parent
#' Correlation coefficient inference via Bayesian updating on hierarchical anatomy
#' 
#' @description Function that computes the posterior probability that the correlation between two structures' properties are greater than some threshold, along with a posterior correlation coefficient, given a tree structure of anatomy. 
#' See the details section below for more information.
#' 
#' @details Correlation matrices constructed from small samples (e.g. n=10) drawn from the same population are generally not very robust, in the sense that across multiple draws, each correlation coefficient varies quite a bit.
#' This is a statistical issue: the distribution of correlation coefficients, when sampling from a population with known correlation, is very wide when the sample size is small. 
#' For example, when sampling the correlation from 10 observations repeatedly drawn from a population correlation of 0.5, 95\% of the distribution p(r) lies between -0.15 and 0.86. 
#' When the number of observations increases to 100, the interval becomes much narrower, around 0.33 to 0.63. 
#' 
#' Typically biological studies are not fortunate to have sample sizes that large. Brain structures are constrained in their correlation patterns however (e.g. two structures that develop from a parent structure are both likely to have similar correlations in volumes to the parent structure); 
#' and exploiting such constraints might provide better estimates of correlations between structural properties than the raw sample correlation coefficient.
#' This function uses a Bayesian framework and exploits the hierarchical structure of anatomy in order to better estimate the correlation coefficient between two given structures. 
#' Briefly, the probability that the correlation between two structures is greater than a threshold value is inferred by recursively computing correlations between their parent structures, and updating this probability (modeled as a beta distribution).
#' A mapping from probability to correlation is also achieved via the Fisher transform to obtain a posterior correlation.
#' 
#' @param hanat Hierarchical anatomy tree giving relationship between structures and their parents.
#' @param struc_i Structure name (must be contained as \code{name} attribute in some node of \code{hanat}).
#' @param struc_j Structure name (must be contained as \code{name} attribute in some node of \code{hanat}).
#' @param indices Optional set of indices that indicate which volumes/observations to use in computing the posterior data. This is useful if you store all volumes across different experimental groups in the tree, and want to separately compute the posterior data in different groups.
#' @param thres_cor Optional, threshold correlation from which to calculate the probability. If \code{NULL}, thres_cor is automatically set as the raw correlation coefficient (all observations, independent of indices chosen).
#' @param precision_cor Optional parameter to tune the width of the distribution of the correlation coefficient p(r). By default this is equal to the number of indices (if provided) or volumes/observations (if indices are not provided).
#' @param precision_sampling Optional parameter to tune the quality of the estimated distributions, which are calculated by sampling a number of times given by this parameter.
#' @param volattr Attribute name in the anatomy tree \code{hanat} that contains the volume data from which correlations are computed.
#' @return a list of posterior data including the probability that the posterior correlation is above the input threshold, and (if \code{p2r.table} is provided), a posterior correlation.
#' 
#' @export
posterior_inference <- function(hanat, struc_i, struc_j, indices=NULL, thres_cor=NULL, precision_cor=NULL, precision_sampling=100, volattr="normVolumes") {
  
  # Set inputs
  inputs <- list(struc_i=struc_i, struc_j=struc_j, indices=indices, thres_cor=thres_cor, precision_cor=precision_cor, precision_sampling=precision_sampling, volattr=volattr)
  
  # Check that structures are different
  if (struc_i==struc_j) {
    stop("Structures must be different")
  }
  
  # Set thres_cor if required
  if (is.null(thres_cor)) {
    thres_cor <- cor(GetAttribute(FindNode(hanat, struc_i), volattr), GetAttribute(FindNode(hanat, struc_j), volattr))
  }
  
  # Set indices
  all.indices <- 1:length(GetAttribute(hanat, volattr))
  if (is.null(indices)) {
    indices <- all.indices
  }
  if (!all(indices %in% all.indices)) {
    stop("indices are not valid")
  }
  
  # Set precision_cor if NULL
  if (is.null(precision_cor)) {
    precision_cor <- length(indices)
  }
  
  # Get list of ancestors for each input structure
  si_ancestors <- rev(names(FindNode(hanat, struc_i)$Get('name', traversal = "ancestor")))
  sj_ancestors <- rev(names(FindNode(hanat, struc_j)$Get('name', traversal = "ancestor")))
  level_smaller <- min(length(si_ancestors), length(sj_ancestors))
  level_larger <- max(length(si_ancestors), length(sj_ancestors))
  
  # Remove overlapping ancestors
  level_merge <- 1
  while(si_ancestors[1]==sj_ancestors[1]) {
    si_ancestors <- si_ancestors[-1]
    sj_ancestors <- sj_ancestors[-1]
    level_merge <- level_merge + 1
  }
  
  # Define function for bootstrapping measured correlation
  boot_cor <- function(hanat, si, sj, volattr, indices, n=1000, apply.fun="median") {
    x <- GetAttribute(FindNode(hanat, si), volattr)[indices]
    y <- GetAttribute(FindNode(hanat, sj), volattr)[indices]
    dfs <- data.frame(x=x, y=y)
    n_obs <- dim(dfs)[1]
    distribution <- replicate(n, cor(dfs[sample(1:n_obs, size = n_obs, replace=TRUE),])[1,2])
    return(do.call(apply.fun, args=list(distribution)))
  }
  
  # Initialize data
  init_cor <- boot_cor(hanat, si_ancestors[1], sj_ancestors[1], volattr, indices)
  init_params <- fit_beta(init_cor, thres_cor, precision_cor = precision_cor, precision_sampling = precision_sampling)
  update_cors_list <- list()
  updated_params_list <- list()
  update_cors_list$init_cor <- init_cor
  updated_params_list$init_params <- init_params
  
  # Set updated params
  updated_params <- init_params
  num_updates <- 0
  
  # Update
  num_iters <- (level_larger - level_merge + 1)
  if (num_iters >= 2) {
    for (iter in 2:num_iters) {
      si_struc <- si_ancestors[min(iter, length(si_ancestors))]
      sj_struc <- sj_ancestors[min(iter, length(sj_ancestors))]
      measured_cor <- boot_cor(hanat, si_struc, sj_struc, volattr, indices)
      update_cors_list[[paste("update", num_updates + 1, sep="_")]] <- measured_cor
      updated_params <- update_beta_prior(a=updated_params[1], b=updated_params[2], rho = measured_cor, threshold = thres_cor, precision_cor = precision_cor, precision_sampling = precision_sampling)
      updated_params_list[[paste("update", num_updates + 1, sep="_")]] <- updated_params
      num_updates <- num_updates + 1
    }
  }
  
  # Outputs
  out <- list()

  out_prior <- list()
  out_prior$r_raw <- cor(GetAttribute(FindNode(hanat, struc_i), volattr)[indices], GetAttribute(FindNode(hanat, struc_j), volattr)[indices])
  out_prior$r_bootstrapped <- measured_cor
  out_prior$proba <- (init_params[1] - 1) / (init_params[1] + init_params[2] - 2)
  out_prior$params <- init_params
  
  out_posterior <- list()
  out_posterior$proba <- (updated_params[1] - 1) / (updated_params[1] + updated_params[2] - 2)
  out_posterior$params <- updated_params
  out_posterior$r <- p2r(out_posterior$proba, cor_thres = thres_cor, n=length(indices))
  
  out_updates <- list()
  out_updates$update_cors_list <- update_cors_list
  out_updates$updated_params_list <- updated_params_list
  out_updates$num_updates <- num_updates
  
  out_options <- list()
  out_options$n_indices <- length(indices)
  out_options$precision_cor <- precision_cor
  out_options$thres_cor <- thres_cor
  
  out$inputs <- inputs
  out$levels <- list(larger=level_larger, smaller=level_smaller, merge=level_merge)
  out$prior <- out_prior
  out$posterior <- out_posterior
  out$updates <- out_updates
  out$options <- out_options
  
  return(out)
}


# Helper functions for updating the beta distribution parameters ----

# Modeling the probability that for a given population correlation rho (d data points), the sample correlation lies above the threshold
# This function determines the beta parameters for that model
fit_beta <- function(rho, threshold, precision_cor=10, precision_sampling=100) {
  cor_dist <- sample_from_cor(rho=rho, d=precision_cor, samples=precision_sampling)
  a <- length(which(cor_dist > threshold))
  b <- precision_sampling - a
  out <- c(a, b)
  names(out) <- c("a", "b")
  return(out)
}

# Helper function to update beta distribution parameters
update_beta <- function(a, b, successes, failures) {
  a_new <- a + successes
  b_new <- b + failures
  out <- c(a_new, b_new)
  names(out) <- c("a", "b")
  return(out)
}

# Helper functions for working with the beta distribution in the context of correlations ----

# Helper function to sample correlation coefficients given a population correlation
sample_from_cor <- function(rho, d=10, samples=100) {
  return(replicate(n=samples, cor(MASS::mvrnorm(d, mu = c(0,0), Sigma = matrix(c(1,rho,rho,1), ncol = 2), empirical = FALSE))[2,1]))
}

# Update the beta parameters given a new measured correlation
update_beta_prior <- function(a, b, rho, threshold, precision_cor=10, precision_sampling=100) {
  cor_dist <- sample_from_cor(rho=rho, d=precision_cor, samples=precision_sampling)
  successes <- length(which(cor_dist > threshold))
  failures <- precision_sampling - successes
  new_params <- update_beta(a=a, b=b, successes=successes, failures=failures)
  return(new_params)
}

# Convert a probability (that r > cor_thres) to a correlation coefficient via the Fisher transform
p2r <- function(prob, cor_thres, n) {
  return(tanh(qnorm(prob, mean=atanh(cor_thres), sd = sqrt(n)/sqrt(n-3))))
}

# Old functions ----

# Transform probability of connection above threshold, to a correlation coefficient
# beta_to_rho <- function(a, b, p2r.table, precision_sampling=1000, retval="median") {
#   probability_values <- rbeta(precision_sampling, a, b)
#   rho_values <- numeric(precision_sampling)
#   for (i in 1:precision_sampling) {
#     probability <- probability_values[i]
#     rho_values[i] <- p2r.table$rho[which.min(abs(p2r.table$prob - probability))]
#   }
#   if (retval=="distribution") {
#     return(rho_values)
#   } else {
#     return(do.call(retval, list(x=rho_values)))
#   }
# }

# Given threshold rho, generate a transformation table mapping probability of being above threshold to correlation coefficient
# Is this just tan(prob) / arctanh(prob), with a scaling factor to account for num_samples??
# TODO: this can be parallelized
# probability_to_rho_table <- function(threshold=0., stepsize=0.01, precision_cor=10, precision_sampling=10000, show_progress=TRUE) {
#   rhos <- seq(from=-1, to=1, by=stepsize)
#   num_rows <- length(rhos)
#   p2r.table <- data.frame(rho=numeric(num_rows), prob=numeric(num_rows))
#   if (show_progress) {
#     pb <- txtProgressBar(min=0, max=num_rows, style = 3)
#   }
#   for (i in 1:length(rhos)) {
#     rho <- rhos[i]
#     probability <- length(which(sample_from_cor(rho, d = precision_cor, samples=precision_sampling) > threshold))/precision_sampling
#     p2r.table$rho[i] <- rho
#     p2r.table$prob[i] <- probability
#     if (show_progress) {
#       setTxtProgressBar(pb, i)
#     }
#   }
#   return(p2r.table)
# }
# 
