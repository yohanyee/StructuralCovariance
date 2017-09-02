invwishart_updating <- function(prior_covmtx, prior_n, update_covmtx, update_n, confidence_intervals=NULL, ci_sample_size=10000) {
  S <- prior_n*prior_covmtx
  v <- dim(S)[1]
  
  out <- list()
  out[["posterior_mean"]] <- (update_covmtx*update_n + S)/(update_n - 1)
  out[["posterior_mode"]] <- (update_covmtx*update_n + S)/(2*v + update_n + 1)
  
  if (is.numeric(confidence_intervals)) {
    posterior_samples <- array(dim = c(dim(S)[1], dim(S)[2], ci_sample_size), dimnames = list(rownames(S), colnames(S), paste("sample", 1:ci_sample_size, sep="_")))
    for (i in 1:ci_sample_size) {
      posterior_samples[,,i] <- LaplacesDemon::rinvwishart(nu = (update_n + v), S = (update_covmtx*update_n + S))
    }
    
    out[["posterior_samples"]] <- posterior_samples
    out[["posterior_confidence_intervals"]] <- list()
    for (cint in confidence_intervals) {
      out[["posterior_confidence_intervals"]][[paste(100*cint, "%", sep="")]] <- apply(posterior_samples, c(1,2), "quantile", cint)
    }
  }
  
  return(out)
}

update_cormatrix_data <- function(base_vols, new_vols, confidence_intervals=NULL, ci_sample_size=10000) {
  prior_covmtx <- cov(base_vols)
  prior_n <- dim(base_vols)[1]
  update_covmtx <- cov(new_vols)
  update_n <- dim(new_vols)[1]
  
  bayes_updated <- invwishart_updating(prior_covmtx=prior_covmtx, prior_n=prior_n, update_covmtx=update_covmtx, update_n=update_n, confidence_intervals=confidence_intervals, ci_sample_size=ci_sample_size)
  ...
  TODO
  ...
}