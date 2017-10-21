library(MASS)

# Helper function to update beta distribution parameters
update_beta <- function(a, b, successes, failures) {
  a_new <- a + successes
  b_new <- b + failures
  out <- c(a_new, b_new)
  names(out) <- c("a", "b")
  return(out)
}

# Helper function to sample correlation coefficients given a population correlation
sample_from_cor <- function(rho, d=10, samples=100) {
  return(replicate(n=samples, cor(mvrnorm(d, mu = c(0,0), Sigma = matrix(c(1,rho,rho,1), ncol = 2), empirical = FALSE))[2,1]))
}

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

# Update the beta parameters given a new measured correlation
update_beta_prior <- function(a, b, rho, threshold, precision_cor=10, precision_sampling=100) {
  cor_dist <- sample_from_cor(rho=rho, d=precision_cor, samples=precision_sampling)
  successes <- length(which(cor_dist > threshold))
  failures <- precision_sampling - successes
  new_params <- update_beta(a=a, b=b, successes=successes, failures=failures)
  return(new_params)
}

# Transform probability of connection above threshold, to a correlation coefficient
# THIS IS WHAT SHOULD BE DONE
beta_to_rho <- function(a, b, p2r.table, precision_sampling=1000, retval="median") {
  probability_values <- rbeta(precision_sampling, a, b)
  rho_values <- numeric(precision_sampling)
  for (i in 1:precision_sampling) {
    probability <- probability_values[i]
    rho_values[i] <- p2r.table$rho[which.min(abs(p2r.table$prob - probability))]
  }
  if (retval=="distribution") {
    return(rho_values)
  } else {
    return(do.call(retval, list(x=rho_values)))
  }
}

# Given threshold rho, generate a transformation table mapping probability of being above threshold to correlation coefficient
# Is this just tan(prob) / arctanh(prob), with a scaling factor to account for num_samples??
# TODO: this can be parallelized
probability_to_rho_table <- function(threshold=0., stepsize=0.01, precision_cor=10, precision_sampling=10000, show_progress=TRUE) {
  rhos <- seq(from=-1, to=1, by=stepsize)
  num_rows <- length(rhos)
  p2r.table <- data.frame(rho=numeric(num_rows), prob=numeric(num_rows))
  if (show_progress) {
    pb <- txtProgressBar(min=0, max=num_rows, style = 3)
  }
  for (i in 1:length(rhos)) {
    rho <- rhos[i]
    probability <- length(which(sample_from_cor(rho, d = precision_cor, samples=precision_sampling) > threshold))/precision_sampling
    p2r.table$rho[i] <- rho
    p2r.table$prob[i] <- probability
    if (show_progress) {
      setTxtProgressBar(pb, i)
    }
  }
  return(p2r.table)
}

# TODO: test on two structures that share the same parent
posterior_inference <- function(hanat, struc_i, struc_j, indices=NULL, p2r.table=NULL, thres_cor=0., precision_cor=10, precision_sampling=100, volattr="normVolumes") {
  
  # Check that structures are different
  if (struc_i==struc_j) {
    stop("Structures must be different")
  }
  
  # Set indices
  all.indices <- 1:length(GetAttribute(hanat, volattr))
  if (is.null(indices)) {
    indices <- all.indices
  }
  if (!all(indices %in% all.indices)) {
    stop("indices are not valid")
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
  out$r_measured_raw <- cor(GetAttribute(FindNode(hanat, struc_i), volattr)[indices], GetAttribute(FindNode(hanat, struc_j), volattr)[indices])
  out$r_measured_bootstrapped <- measured_cor
  out$posterior_proba <- (updated_params[1] - 1) / (updated_params[1] + updated_params[2] - 2)
  out$posterior_params <- updated_params
  out$update_cors_list <- update_cors_list
  out$updated_params_list <- updated_params_list
  out$num_updates <- num_updates
  out$levels <- list(larger=level_larger, smaller=level_smaller, merge=level_merge)
  
  if (is.null(p2r.table)) {
    out$r_posterior <- NA
  } else {
    out$r_posterior <- beta_to_rho(a=updated_params[1], b=updated_params[2], p2r.table=p2r.table, precision_sampling = precision_sampling*10, retval = "median")
  }
  return(out)
}
