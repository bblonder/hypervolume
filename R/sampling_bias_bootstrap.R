sampling_bias_bootstrap <- function(name, hv, n = 10, points_per_resample = 'sample_size', cores = 1, verbose = TRUE, mu = NULL, sigma = NULL, cols_to_bias = 1:ncol(hv@Data), weight_func = NULL) {
  # Check if cluster registered to doparallel backend exists
  exists_cluster = TRUE
  if(cores > 1 & getDoParWorkers() == 1) {
    # If no cluster is registered, create a new one based on use input
    cl = makeCluster(cores)
    clusterEvalQ(cl, {
      library(hypervolume)
      library(mvtnorm)
    })
    registerDoParallel(cl)
    exists_cluster = FALSE
  }
  
  # Create folder to store bootstrapped hypervolumes
  dir.create(file.path('./Objects', name))
  if(verbose) {
    pb = progress_bar$new(total = n)
  }
  
  # Apply weights to data before bootstrapping
  foreach(i = 1:n, .combine = c) %dopar% {
    if(is.null(weight_func)) {
      if(length(mu) == 1) {
        weights = dnorm(hv@Data[,cols_to_bias], mean = mu, sd = sqrt(sigma))
      } else {
        weights = dmvnorm(hv@Data[,cols_to_bias], mean = mu, sigma = diag(sigma))
      }
    } else {
      weights = weight_func(hv@Data[,cols_to_bias])
    }
    if(points_per_resample == 'sample_size') {
      points = apply(rmultinom(nrow(hv@Data), 1, weights) == 1, 2, which)
    } else {
      points = apply(rmultinom(points_per_resample, 1, weights) == 1, 2, which)
    }
    sample_dat = hv@Data[points,]
    h = copy_param_hypervolume(hv, sample_dat, name = paste("resample", as.character(i)))
    path = paste0(h@Name, '.rds')
    saveRDS(h, file.path('./Objects', name, path))
    if(verbose) {
      pb$tick()
    }
  }
  
  # If a cluster was created for this specific function call, close cluster and register sequential backend
  if(!exists_cluster) {
    stopCluster(cl)
    registerDoSEQ()
  }
  return(file.path(getwd(), 'Objects', name))
}