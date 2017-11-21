library(igraph)

clusters_initialize_random <- function(n, k) {
  membership <- sample(1:k, size=n, replace=TRUE)
  return(membership)
}
  
clusters_flip <- function(membership, size=1) {
  indices_to_flip <- sample(1:length(membership), size=size, replace=FALSE)
  unique_classes <- unique(membership)
  for (i in indices_to_flip) {
    new_class <- sample(setdiff(unique_classes, membership[i]), size=1)
    membership[i] <- new_class
  }
  return(membership)
}

energy_matrix_distance <- function(vols, membership, lambda=1) {
  k <- length(unique(membership))
  p <- dim(vols)[2]
  num_features <- p*(p-1)/2
  cormats_squeezed <- matrix(nrow=k, ncol=num_features)
  for (cl in 1:k) {
    cormat <- cor(vols[which(membership==cl),])
    corvec <- cormat[upper.tri(cormat)]
    cormats_squeezed[cl,] <- corvec
  }
  mean_dist <- mean(dist(cormats_squeezed))
  energy <- 1/mean_dist
  size_penalty <- exp(-min(table(membership))/10)
  return((energy + lambda*size_penalty))
}

energy_integration_segregation <- function(vols, membership, threshold=0.6) {
  k <- length(unique(membership))
  df <- data.frame(integration=numeric(k), segregation=numeric(k))
  for (cl in 1:k) {
    cormat <- cor(vols[which(membership==cl),])
    gph <- graph.adjacency(cormat > threshold)
    df$integration[cl] <- transitivity(gph)
    df$segregation[cl] <- modularity(cluster_walktrap(gph))
  }
  return(det(0.001/cov(df)))
}

cluster_individuals <- function(vols, energy_function, cluster_size=9, batch_size=10, temperature=1, max_iter=100000, tol=1e-03, r=100, ...) {
  n <- dim(vols)[1]
  k <- cluster_size
  inverse_temperature <- 1/temperature
  membership <- clusters_initialize_random(n, k)
  this_energy <- do.call(energy_function, args = list(vols=vols, membership=membership,...))
  
  print(paste("Iteration:", "0", "| Energy:", this_energy, "| Slope:", NA))

  energy_vec <- numeric(max_iter+1)
  energy_vec[1] <- this_energy
  
  iter <- 1
  continue <- TRUE
  while (continue) {
    while (TRUE) {
      new_membership <- clusters_flip(membership, size=batch_size)
      new_membership_num_classes <- length(unique(new_membership))
      if (new_membership_num_classes==k) {
        break
      }
    }
    new_energy <- do.call(energy_function, args = list(vols=vols, membership=new_membership,...))
    if (new_energy < this_energy) {
      membership <- new_membership
      this_energy <- new_energy
    } else {
      prob_of_flip <- exp(-inverse_temperature*(new_energy - this_energy))
      flip <- sample(c(TRUE, FALSE), size=1, prob = c(prob_of_flip, 1-prob_of_flip))
      if (flip) {
        membership <- new_membership
        this_energy <- new_energy
      }
    }
    
    energy_vec[iter+1] <- this_energy
    
    if ((iter %% r)==0) {
      slope <- lm(energy_vec[(iter-r+1):(iter)] ~ c(1:r))$coefficients[2]
      print(paste("Iteration:", iter, "| Energy:", this_energy, "| Slope:", slope))
      if (slope >= tol) {
        continue <- FALSE
        print("Stopping: no energy change.")
      }
    }
    if (iter>=max_iter) {
      continue <- FALSE
      print("Stopping: reached maximum iterations.")
    }
    iter <- iter + 1
  }
  iter_vec <- 0:(iter-1)
  energy_vec <- energy_vec[1:iter]
  
  return(list(membership=membership, iteration_energies=data.frame(iteration=iter_vec, energy=energy_vec)))
}

res <- cluster_individuals(vols, "energy_matrix_distance", cluster_size = 20, batch_size = 100, tol=1e-03, temperature = 100, max_iter = 1000)

res <- cluster_individuals(vols, "energy_integration_segregation", cluster_size = 20, batch_size = 100, tol=1e-03, temperature = 1000, threshold=0.7, max_iter = 10000, r=10)
a <- res$iteration_energies
qplot(a$iteration, a$energy)

membership <- res$membership

scans[which(membership==1),c("Mouse_ID", "Study_Name", "Genotype")]

sf <- scans[,c("Study_Name", "Is_Wildtype")]
sf$membership <- membership
sf$Is_Wildtype <- factor(sf$Is_Wildtype)

geno_table <- table(sf$Is_Wildtype, sf$membership)
study_table <- table(sf[which(sf$Is_Wildtype=="MUT"),]$Study_Name, sf[which(sf$Is_Wildtype=="MUT"),]$membership)

study_table

chisq.test(geno_table)
chisq.test(study_table)
