bootstrap_n_hypervolumes <- function(name, hv, n = 10, points_per_resample = 'sample_size', cores = 1, verbose = TRUE, seed = NULL) {
  
  # check if hv is an object of class Hypervolume
  # if yes, transform it to a list
  if(inherits(hv, "Hypervolume")){
    hv <- hypervolume_join(hv)
  }
  
  # get random state from global environment
  old_state <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  on.exit({
    # restore random state
    if (!is.null(old_state)) {
      assign(".Random.seed", old_state, envir = .GlobalEnv, inherits = FALSE)
    }
  }, add = TRUE)
  
  
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
  
  on.exit({
    # If a cluster was created for this specific function call, close cluster and register sequential backend
    if(!exists_cluster) {
      stopCluster(cl)
      registerDoSEQ()
    }
  }, add = TRUE)
  
  
  # Create folder to store bootstrapped hypervolumes
  dir.create(file.path('./Objects', name), recursive = TRUE, showWarnings = FALSE)

  # check if names are equal, if yes, add a number to make them unique
  hv_names <- unlist(lapply(hv@HVList, function(x) x@Name))
  if(length(unique(hv_names)) != length(hv@HVList)){
    hv_names <- paste(hv_names, 1:length(hv@HVList), sep = "_")
  }
  
  # write the original order of hypervolumes, needed to avoid poor alignment when using
  # hypervolume_n_occupancy
  write.table(hv_names, file = file.path('./Objects', name, "log.txt"), row.names = FALSE)
  
  if(verbose){
    Fun <- function(...) progress_bar_foreach(iterator = n, fun = function (a, ...)c(a, list(...)), clear = FALSE) 
  } else {
    Fun <- function(...) function (a, ...)c(a, list(...))
  }
  
  
  for(j in 1:length(hv@HVList)){
    # Construct n hypervolumes from points_per_sample points sampled with replacement from original data
    
    path_j <- file.path('./Objects', name, hv_names[j])
    name_j <- dir.create(path_j, showWarnings = FALSE)
    hv_j <- hv@HVList[[j]]
    
    
    if(verbose){
      cat("\n Permuting hypervolume", hv_names[j], "...\n")
    }

    foreach(i = 1:n, .combine = Fun()) %dopar% {
      if(points_per_resample == 'sample_size') {
        sample_dat = hv_j@Data[sample(1:nrow(hv_j@Data), nrow(hv_j@Data), replace = TRUE),]
      } else {
        sample_dat = hv_j@Data[sample(1:nrow(hv_j@Data), points_per_resample, replace = TRUE),]
      }
      h = suppressMessages(copy_param_hypervolume(hv_j, sample_dat, name = paste("resample", as.character(i)))) 
      path = paste0(h@Name, '.rds')
      saveRDS(h, file.path(path_j, path))
    }
  }


  # Absolute path to hypervolume objects
  return(file.path(getwd(), 'Objects', name))
}
