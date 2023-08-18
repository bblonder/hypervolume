hypervolume_overlap_confidence <- function(path1, path2, CI = .95, cores = 1) {
  exists_cluster = TRUE
  if(cores > 1 & getDoParWorkers() == 1) {
    # If no cluster is registered, create a new one based on user input
    cl = makeCluster(cores)
    clusterEvalQ(cl, {
      library(hypervolume)
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
  
  if(list.files(path1)[1] != "resample 1.rds" | list.files(path2)[1] != "resample 1.rds") {
    if(!exists_cluster) {
      stopCluster(cl)
      registerDoSEQ()
    }
    stop("Invalid input")
  }
  distribution = foreach(i = list.files(path1), .combine = rbind) %:%
    foreach(j = list.files(path2), .combine = rbind) %dopar% {
      h1 = readRDS(file.path(path1, i))
      h2 = readRDS(file.path(path2, j))
      hypervolume_overlap_statistics(hypervolume_set(h1, h2, check.memory = FALSE))
    }
  results = list(
    "jaccard" = quantile(distribution[,"jaccard"], c(.5 - CI/2, .5 + CI/2)),
    "sorensen" = quantile(distribution[,"sorensen"], c(.5 - CI/2, .5 + CI/2)),
    "frac_unique_1" = quantile(distribution[,"frac_unique_1"], c(.5 - CI/2, .5 + CI/2)),
    "frac_unique_2" = quantile(distribution[,"frac_unique_2"], c(.5 - CI/2, .5 + CI/2)),
    "distribution" = distribution
  )
  
  
  return(results)
}
