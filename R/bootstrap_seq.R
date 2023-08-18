bootstrap_seq <- function(name, hv, n = 10, seq = 3:nrow(hv@Data), cores = 1, verbose = TRUE) {
  # Check if cluster registered to doparallel backend exists
  exists_cluster = TRUE
  if(cores > 1 & getDoParWorkers() == 1) {
    # If no cluster is registered, create a new one based on user input
    cl = makeCluster(cores)
    clusterEvalQ(cl, {
      library(hypervolume)
      library(foreach)
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
  
  # Create folder to store bootstrapped hypervolumes
  dir.create(file.path('./Objects', name))
  
  # Construct a hypervolume from samples of increasingly many points and bootstrap new hypervolume
  foreach(i = seq) %dopar% {
    subdir = paste('sample size', as.character(i))
    dir.create(file.path('./Objects', name, subdir))
    bootstrap(file.path(name, subdir), hv, n, i, verbose = verbose)
  }
  
  
  return(file.path(getwd(), 'Objects', name))
}