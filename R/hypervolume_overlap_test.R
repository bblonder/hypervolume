hypervolume_overlap_test <- function(hv1, hv2, path, alternative = "one-sided", bins = 100, cores = 1) {
  if(alternative != "one-sided" & alternative != "two-sided") {
    stop("invalid alternative hypothesis")
  }
  observed = hypervolume_overlap_statistics(hypervolume_set(hv1, hv2, check.memory = FALSE))
  
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
  
  if(length(path) == 1) {
    if(length(list.files(path)) == 0) {
      if(!exists_cluster) {
        stopCluster(cl)
        registerDoSEQ()
      }
      stop("Invalid input path")
    }
    if(list.files(path)[1] == "permutation1") {
      distribution = foreach(i = list.files(path), .combine = rbind) %dopar% {
        h1 = readRDS(file.path(path, i, "hv1.rds"))
        h2 = readRDS(file.path(path, i, "hv2.rds"))
        hypervolume_overlap_statistics(hypervolume_set(h1, h2, check.memory = FALSE))
      }
    } else if(list.files(path)[1] == "resample 1.rds") {
      files = list.files(path)
      half = floor(length(files)/2)
      distribution = foreach(i = 1:half, .combine = rbind) %:%
        foreach(j = (half+1):length(files), .combine = rbind) %dopar% {
          h1 = readRDS(file.path(path, files[i]))
          h2 = readRDS(file.path(path, files[j]))
          hypervolume_overlap_statistics(hypervolume_set(h1, h2, check.memory = FALSE))
        }
    } else {
      if(!exists_cluster) {
        stopCluster(cl)
        registerDoSEQ()
      }
      stop("Invalid input path")
    }
  } else if(length(path) == 2) {
    if(list.files(path[1])[1] == list.files(path[2])[1] & list.files(path[1])[1] == "resample 1.rds") {
      distribution = foreach(i = list.files(path[1]), .combine = rbind) %:%
        foreach(j = list.files(path[2]), .combine = rbind) %dopar% {
          h1 = readRDS(file.path(path[1], i))
          h2 = readRDS(file.path(path[2], j))
          hypervolume_overlap_statistics(hypervolume_set(h1, h2, check.memory = FALSE))
        }
    } else {
      if(!exists_cluster) {
        stopCluster(cl)
        registerDoSEQ()
      }
      stop("invalid input paths")
    }
  } else {
    if(!exists_cluster) {
      stopCluster(cl)
      registerDoSEQ()
    }
    stop("invalid input paths")
  }
  
  # Build list of useful results
  if(alternative == "one-sided") {
    p_values = list(
      "jaccard" = mean(distribution[,"jaccard"] <= observed["jaccard"]),
      "sorensen" = mean(distribution[,"sorensen"] <= observed["sorensen"]),
      "frac_unique_1" = mean(distribution[,"frac_unique_1"] >= observed["frac_unique_1"]),
      "frac_unique_2" = mean(distribution[,"frac_unique_2"] >= observed["frac_unique_2"])
    )
  } else if(alternative == "two-sided") {
    p_values = list(
      "jaccard" = mean(abs(distribution[,"jaccard"] - mean(distribution[,"jaccard"])) >= observed["jaccard"]),
      "sorensen" = mean(abs(distribution[,"sorensen"] - mean(distribution[,"sorensen"])) >= observed["sorensen"]),
      "frac_unique_1" = mean(abs(distribution[,"frac_unique_1"] - mean(distribution[,"frac_unique_1"])) >= observed["frac_unique_1"]),
      "frac_unique_2" = mean(abs(distribution[,"frac_unique_2"] - mean(distribution[,"frac_unique_2"])) >= observed["frac_unique_2"])
    )
  }
  plots = list(
    "jaccard" = ggplot(data.frame(distribution), aes(x = jaccard, y = ..density..)) + 
      geom_histogram(bins = bins, alpha = .7) + 
      geom_vline(aes(xintercept = observed["jaccard"]), color = "red") + 
      ggtitle("Distribution of jaccard index"),
    "sorensen" = ggplot(data.frame(distribution), aes(x = sorensen, y = ..density..)) + 
      geom_histogram(bins = bins, alpha = .7) + 
      geom_vline(aes(xintercept = observed["sorensen"]), color = "red") + 
      ggtitle("Distribution of sorensen index"),
    "frac_unique_1" = ggplot(data.frame(distribution), aes(x = frac_unique_1, y = ..density..)) + 
      geom_histogram(bins = bins, alpha = .7) + 
      geom_vline(aes(xintercept = observed["frac_unique_1"]), color = "red") + 
      ggtitle("Distribution of fraction of Hypervolume 1 that is unique"),
    "frac_unique_2" = ggplot(data.frame(distribution), aes(x = frac_unique_2, y = ..density..)) + 
      geom_histogram(bins = bins, alpha = .7) + 
      geom_vline(aes(xintercept = observed["frac_unique_2"]), color = "red") + 
      ggtitle("Distribution of fraction of Hypervolume 2 that is unique")
  )
  result = list(p_values = p_values, plots = plots, distribution = distribution)
  
  if(!exists_cluster) {
    stopCluster(cl)
    registerDoSEQ()
  }
  
  return(result)
}