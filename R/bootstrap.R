bootstrap <- function(name, hv, n = 10, points_per_resample = 'sample_size', cores = 1, verbose = TRUE) {
  # Check if cluster registered to doparallel backend exists
  exists_cluster = TRUE
  if(cores > 1 & getDoParWorkers() == 1) {
    # If no cluster is registered, create a new one based on use input
    cl = makeCluster(cores)
    clusterEvalQ(cl, {
      library(hypervolume)
    })
    registerDoParallel(cl)
    exists_cluster = FALSE
  }
  
  # Create folder to store bootstrapped hypervolumes
  dir.create(file.path('./Objects', name))
  if(verbose) {
    pb = progress_bar$new(total = n)
  }
  
  # Construct n hypervolumes from points_per_sample points sampled with replacement from original data
  foreach(i = 1:n) %dopar% {
    if(points_per_resample == 'sample_size') {
      sample_dat = hv@Data[sample(1:nrow(hv@Data), nrow(hv@Data), replace = TRUE),]
    } else {
      sample_dat = hv@Data[sample(1:nrow(hv@Data), points_per_resample, replace = TRUE),]
    }
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
  
  # Absolute path to hypervolume objects
  return(file.path(getwd(), 'Objects', name))
}