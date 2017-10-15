require(MASS)

# Helper function to update beta distribution parameters
update_beta <- function(a, b, successes, failures) {
  a_new <- a + successes
  b_new <- b + failures
  return(c(a_new, b_new))
}

# Helper function to sample correlation coefficients given a population correlation
sample_from_cor <- function(rho, d=10, size=100) {
  return(replicate(n=size, cor(mvrnorm(d, mu = c(0,0), Sigma = matrix(c(1,rho,rho,1), ncol = 2), empirical = FALSE))[2,1]))
}

# Modeling the probability that for a given population correlation rho (d data points), the sample correlation lies above the threshold
# This function determines the beta parameters for that model
fit_beta <- function(rho, threshold, num_samples=10, precision=100) {
  cor_dist <- sample_from_cor(rho=rho, d=num_samples, size=precision)
  a <- length(which(cor_dist > threshold))
  b <- precision - a
  return(c(a, b))
}

# Update the beta parameters given a new measured correlation
update_beta_prior <- function(a, b, rho, threshold, num_samples=10, precision=100) {
  cor_dist <- sample_from_cor(rho=rho, d=num_samples, size=precision)
  successes <- length(which(cor_dist > threshold))
  failures <- precision - successes
  new_params <- update_beta(a=a, b=b, successes=successes, failures=failures)
  return(new_params)
}

# Given threshold rho, generate a transformation table mapping probability of being above threshold to correlation coefficient
# Is this just tan(prob) / arctanh(prob), with a scaling factor to account for num_samples??
# TODO: this can be parallelized
probability_to_rho_table <- function(threshold, stepsize=0.01, num_samples=10, precision=10000, show_progress=TRUE) {
  rhos <- seq(from=-1, to=1, by=stepsize)
  num_rows <- length(rhos)
  p2r_table <- data.frame(rho=numeric(num_rows), prob=numeric(num_rows))
  if (show_progress) {
    pb <- txtProgressBar(min=0, max=num_rows, style = 3)
  }
  for (i in 1:length(rhos)) {
    rho <- rhos[i]
    sample_from_cor(rho, d = num_samples, size=precision)
    probability <- length(which(sample_from_cor(rho, d = num_samples, size=precision) > threshold))/precision
    p2r_table$rho[i] <- rho
    p2r_table$prob[i] <- probability
    if (show_progress) {
      setTxtProgressBar(pb, i)
    }
  }
  return(p2r_table)
}

# Transform probability of connection above threshold, to a correlation coefficient
# THIS IS WHAT SHOULD BE DONE
beta_to_rho <- function(a, b, p2r_table, precision=100, retval="median") {
  probability_values <- rbeta(precision, a, b)
  rho_values <- numeric(precision)
  for (i in 1:precision) {
    probability <- probability_values[i]
    rho_values[i] <- p2r_table$rho[which.min(abs(p2r_table$prob - probability))]
  }
  if (retval=="distribution") {
    return(rho_values)
  } else {
    return(do.call(retval, list(x=rho_values)))
  }
}

# # setup
# thres_cor <- 0.3
# num_mice <- 10
# precision <- 100
# 
# # Precompute table for determining correlation from probability
# # Note that the precision here is manually set, because precision in this context is different from the other uses
# # Here, we'd like to be as precise as possible in mapping probability to correlations, as opposed to precision=a+b as a measure of uncertainty of the beta distribution
# p2r_table <- probability_to_rho_table(threshold=thres_cor, stepsize = 0.01, num_samples = num_mice, precision = 10000)
# 
# # prior
# init_cor <- 1
# prior_params <- fit_beta(init_cor, thres_cor, num_samples=num_mice, precision=precision)
# hist(rbeta(10000, prior_params[1], prior_params[2]), 100)
# prior_cor <- beta_to_rho(a=prior_params[1], b=prior_params[2], p2r_table=p2r_table, precision = precision, retval = "median")
# print(prior_cor)
# 
# # measured cor
# measured_cor <- 0.01
# updated_params <- update_beta_prior(a=prior_params[1], b=prior_params[2], rho = measured_cor, threshold = thres_cor, num_samples = num_mice, precision = precision)
# hist(rbeta(10000, updated_params[1], updated_params[2]), 100)
# posterior_cor <- beta_to_rho(a=updated_params[1], b=updated_params[2], p2r_table=p2r_table, precision = precision, retval = "median")
# print(posterior_cor)
# 
# # another measured cor
# measured_cor <- 0.01
# prior_params <- updated_params
# updated_params <- update_beta_prior(a=prior_params[1], b=prior_params[2], rho = measured_cor, threshold = thres_cor, num_samples = num_mice, precision = precision)
# hist(rbeta(10000, updated_params[1], updated_params[2]), 100)
# posterior_cor <- beta_to_rho(a=updated_params[1], b=updated_params[2], p2r_table=p2r_table, precision = precision, retval = "median")
# print(posterior_cor)
# 
# # another measured cor
# measured_cor <- 0.01
# prior_params <- updated_params
# updated_params <- update_beta_prior(a=prior_params[1], b=prior_params[2], rho = measured_cor, threshold = thres_cor, num_samples = num_mice, precision = precision)
# hist(rbeta(10000, updated_params[1], updated_params[2]), 100)
# posterior_cor <- beta_to_rho(a=updated_params[1], b=updated_params[2], p2r_table=p2r_table, precision = precision, retval = "median")
# print(posterior_cor)
# 
# # another measured cor
# measured_cor <- 0.5
# prior_params <- updated_params
# updated_params <- update_beta_prior(a=prior_params[1], b=prior_params[2], rho = measured_cor, threshold = thres_cor, num_samples = num_mice, precision = precision)
# hist(rbeta(10000, updated_params[1], updated_params[2]), 100)
# posterior_cor <- beta_to_rho(a=updated_params[1], b=updated_params[2], p2r_table=p2r_table, precision = precision, retval = "median")
# print(posterior_cor)
