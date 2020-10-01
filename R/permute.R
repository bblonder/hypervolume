#' Permutations of hypervolumes from two input hypervolumes
#' 
#' @param name Name of directory inside ./Objects for storing hypervolume objects
#' @param hv1 A hypervolume
#' @param hv2 A hypervolume
#' @param n Number of permutations to take
#' @param cores size of cluster to register to doParallel backend
#' @param verbose A boolean, prints a progress bar if run sequentially
#' @return absolute path to a directory of permuted hypervolumes
#' @example 
#' data("quercus")
#' qsample1 = quercus[sample(1:nrow(quercus), 150),]
#' qsample2 = quercus[sample(1:nrow(quercus), 150),]
#' hv1 = hypervolume(qsample1[,2:3])
#' hv2 = hypervolume(qsample2[,2:3])
#' perm_path = permute("Quercus_perm_150", hv1, hv2, n = 200, cores = 20)
#' 
#' perm_null = foreach(i = list.files(perm_path), .combine = rbind) %dopar% {
#' h1 = readRDS(file.path(perm_path, i, "hv1.rds"))
#' h2 = readRDS(file.path(perm_path, i, "hv2.rds"))
#' hypervolume_overlap_statistics(hypervolume_set(h1, h2, check.memory = FALSE))
#' }

permute <- function(name, hv1, hv2, n = 50, cores = 1, verbose = TRUE) {
  # Check if cluster registered to doparallel backend exists
  exists_cluster = TRUE
  if(cores > 1 & getDoParWorkers() == 1) {
    # If no cluster is registered, create a new one based on use input
    cl = makeCluster(cores)
    clusterEvalQ(cl, {
      library(hypervolume)
      source('Utils.R')
    })
    registerDoParallel(cl)
    exists_cluster = FALSE
  }
  
  # Create folder to store permuted hypervolumes
  dir.create(file.path('./Objects', name))
  if(verbose) {
    pb = progress_bar$new(total = n)
  }
  
  # Create n pairs of hypervolumes by permuting combined data of hv1 and hv2
  foreach(i = 1:n, .combine = c) %dopar% {
    combined_data = rbind(hv1@Data, hv2@Data)
    
    # take a sample of points to include in first permuted hypervolume.
    perm_idx = sample(1:nrow(combined_data), nrow(hv1@Data))
    
    # Use sampled indices to construct first hypervolume and the rest of the data to construct the second hypervolume
    h1 = copy_param_hypervolume(hv1, combined_data[perm_idx,], name = "hv1")
    h2 = copy_param_hypervolume(hv2, combined_data[-1 * perm_idx,], name = "hv2")
    
    subdir = file.path("./Objects", name, paste0("permutation", as.character(i)))
    # Save files to subdirectories of ./Object/<name> called permutation <i>
    dir.create(subdir)
    name1 = paste0(h1@Name, '.rds')
    name2 = paste0(h2@Name, '.rds')
    saveRDS(h1, file.path('./Objects', name, paste0("permutation", as.character(i)), name1))
    saveRDS(h2, file.path('./Objects', name, paste0("permutation", as.character(i)), name2))
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