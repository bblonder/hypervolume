weighted_bootstrap <- function(name, hv, n = 10, points_per_resample = 'sample_size', cores = 1, verbose = TRUE, to_file = TRUE, mu = NULL, sigma = NULL, cols_to_weigh = 1:ncol(hv@Data), weights = NULL) {
  # Check if cluster registered to doparallel backend exists
  exists_cluster = TRUE
  if(cores > 1 & getDoParWorkers() == 1) {
    # If no cluster is registered, create a new one based on user input
    cl = makeCluster(cores)
    clusterEvalQ(cl, {
      library(hypervolume)
      library(mvtnorm)
    })
    registerDoParallel(cl)
    exists_cluster = FALSE
  }
  if(cores >1 & getDoParWorkers() > 1) {
    print("Exsisting cluster registered to doParallel. The existing cluster will be used instead of user specified number of cores.")
  }
  
  on.exit({
    # If a cluster was created for this specific function call, close cluster and register sequential backend
    if(!exists_cluster) {
      stopCluster(cl)
      registerDoSEQ()
    }
  })
  
  if(to_file){
    # Create folder to store bootstrapped hypervolumes
    dir.create(file.path('./Objects', name))
  } else {
    hv_list = new("HypervolumeList")
  }
  if(verbose) {
    pb = progress_bar$new(total = n)
  }
  
  # Calculate weights from data before bootstrapping
  list = foreach(i = 1:n, .combine = c) %dopar% {
    if(is.null(weights)) {
      if(length(mu) == 1) {
        weights = dnorm(hv@Data[,cols_to_weigh], mean = mu, sd = sqrt(sigma))
      } else {
        weights = dmvnorm(hv@Data[,cols_to_weigh], mean = mu, sigma = diag(sigma))
      }
    } else {
      weights = weights
    }
    if(points_per_resample == 'sample_size') {
      points = sample(1:nrow(hv@Data), size = nrow(hv@Data), replace = TRUE, prob = weights)
    } else {
      points = sample(1:nrow(hv@Data), size = points_per_resample, replace = TRUE, prob = weights)
    }
    sample_dat = hv@Data[points,]
    h = copy_param_hypervolume(hv, sample_dat, name = paste("resample", as.character(i)))
    if(to_file) {
      path = paste0(h@Name, '.rds')
      saveRDS(h, file.path('./Objects', name, path))
    } 
    if(verbose) {
      pb$tick()
    }
    if(!to_file) {
      h
    } 
  }
  
  if(to_file) {
    return(file.path(getwd(), 'Objects', name))
  } else {
    hv_list@HVList = c(list)
    return(hv_list)
  }
}