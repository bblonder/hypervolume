k_split <- function(name, hv, k = 5, cores = 1, verbose = TRUE) {
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
  
  # Leaves out 1/k of the data points for each resample. Each set of 1/k points is disjoint.
  npoints = nrow(hv@Data)
  shuffled = hv@Data[sample(1:npoints, npoints, replace = FALSE),]
  if(verbose) {
    pb = progress_bar$new(total = k*3)
  }
  foreach(i = 1:k, .combine = c) %dopar% {
    range = (floor((i-1)*npoints/k) + 1):floor(i*npoints/k)
    train_dat = shuffled[-1 * range,]
    h = copy_param_hypervolume(hv, data = train_dat, name = paste("split", as.character(i)))
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
