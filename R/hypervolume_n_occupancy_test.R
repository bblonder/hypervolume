hypervolume_n_occupancy_test <- function(observed, path, alternative = "two_sided", CI = 0.95, cores = 1) {
  
  # intitialize multi core calculations
  
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
  
  # check if the first element of the first pairwise combination is called permutation1.rds
  # n <- list all the files in the first pa

  if(list.files(file.path(path, list.files(path)[1]))[1] != "permutation1.rds" ) {
    if(!exists_cluster) {
      stopCluster(cl)
      registerDoSEQ()
    }
    stop("Invalid input")
  }


  # name of the groups of the observed HypervolumeList
  group_names <- unlist(lapply(observed@HVList, function(x)x@Name))
  
  # combination of these names
  observed_combn <- combn(group_names, 2)
  
  # number of the pairwise combination inferred from the parent directory
  groups <- length(list.files(path))
  
  # check if the number of columns of observed_combn is the same as the one in groups
  if(ncol(observed_combn) != groups){
    stop("The number of pairwise combination does not match between observed and path.")
  }
  
  # store the names of the permuted files. Names are the same for each pairwise combination
  n <- list.files(file.path(path, list.files(path)[1]))
  
  # initialize an empty list for storing the pairwise differences between
  # ValueAtRandomPoints of the observed argument
  observed_res <- vector(ncol(observed_combn), mode = "list")
  
  # calculate pairwise differences
  for(i in 1:ncol(observed_combn)){
    group_1 <- which(group_names == observed_combn[1, i])
    group_2 <- which(group_names == observed_combn[2, i])
    result <- lapply(observed[[c(group_1, group_2)]]@HVList, function(x) x@ValueAtRandomPoints)
    result <- do.call(cbind, result)
    observed_res[[i]] <- result[ , 1] - result [ , 2]
  }
  
  

  # calculate the difference in the resampled hypervolumes
  # perform a nested loop with foreach and store the results
  # this loop will return a list of length equal to the number of pairwise combinations
  # each list contains a marix, where the rose are the RandomPoints and columns represent
  # the difference between groups
  
  simu_res <- foreach(i = list.files(path)) %:% foreach(j = n, .combine = cbind) %dopar% {
    h1 = readRDS(file.path(path, i, j))
    result <- lapply(h1@HVList, function(x) x@ValueAtRandomPoints)
    result <- do.call(cbind, result)
    simulated_res <- result[ , 1] - result [ , 2]
  }
  
  
  # initialize a list to store the results
  result_list <- vector(groups, mode = "list")
  
  # number of permutations
  n <- length(list.files(file.path(path, list.files(path)[1])))
  
  # calculate the probability that observed differences are significantly different
  # from simulated differences
  # retain only ValueAtRandomPoints for which the difference is significant
  # se the other ValueAtRandomPoints to 0
  
  for(i in 1:groups){
    # simulated results for the i-th combination
    distribution <- simu_res[[i]]
    
    # observed result for the i-th combination
    obs <- observed_res[[i]]
    
    # calculates probabilities: can be more, less or two_sided
    # credits to vegan oecosimu
    if(alternative == "less"){
      result <- sweep(distribution, 1, obs, ">=")
      result <- apply(result, 1, sum)
      p <- pmin(1, (result  + 1)/(n + 1))
    }
    
    if(alternative == "more"){
      result <- sweep(distribution, 1, obs, "<=")
      result <- apply(result, 1, sum)
      p <- pmin(1, (result  + 1)/(n + 1))
    }
    
    if(alternative == "more_less"){
      result_less <- sweep(distribution, 1, obs, ">=")
      result_more <- sweep(distribution, 1, obs, "<=")
      result_less <- apply(result_less, 1, sum)
      result_more <- apply(result_more, 1, sum)
      p_less <- pmin(1, (result_less  + 1)/(n + 1))
      p_more <- pmin(1, (result_more  + 1)/(n + 1))
      p <- apply(data.frame(p_less, p_more), 1, max)
    }
    
    
    if(alternative == "two_sided"){
      result <- sweep(abs(distribution), 1, abs(obs), ">=")
      result <- apply(result, 1, sum)
      p <- pmin(1, (result  + 1)/(n + 1))
    }
    
    # retain only those probabilities that are greater than CI
    p_obs <- p >=  CI
    
    # retain only ValueAtRandomPoints for which the difference is significant
    if(sum(p_obs) == 0){
      result_list[[i]] <- rep(0, length(obs))
    } else {
      obs[! p_obs] <- 0
      result_list[[i]] <- obs
    }
  }
  

  # labels of the pairwise combinations
  label_res <- unlist(apply(observed_combn, 2, function(x) paste(x[1], x[2], sep = "__")))
  
  # cbind the results of pairwise combinations
  result_combn <- do.call(cbind, result_list)
  
  # assign colnames
  colnames(result_combn) <- label_res

  
  # create an empty list to store the hypervolume results 
  hv_list_res <- vector(ncol(observed_combn), mode = "list")
  
  # create an empty hypervolume as a model for storing results
  empty_hypervolume <- new("Hypervolume")
  
  # get mindensity from observed, it is the same ofor all the hypervolumes
  mindensity <- observed@HVList[[1]]@PointDensity
  
  
  vol_list <- unlist(lapply(observed@HVList, function(x) x@Volume))
  res <- lapply(observed@HVList, function(x) x@ValueAtRandomPoints)
  res <- do.call(cbind, res)
  res[res > 0] <- 1
  
  # volume of the significant fraction for th groups under comparison
  intersection_weights <- sweep(res, 1, apply(res, 1, sum), "/")
  vol_tot <- sum(apply(intersection_weights, 2, function(x) mean(x[x > 0], na.rm = TRUE)) * vol_list, na.rm = TRUE)

  
  # dimnames
  cn <- dimnames(observed@HVList[[1]]@RandomPoints)[[2]]
  dn <- list(NULL, cn)
  
  # create the hypervolumes to return to the user
  for(i in 1:ncol(observed_combn)){
    
    # position of the groups under comparison in observed_combn
    group_1 <- which(group_names == observed_combn[1, i])
    group_2 <- which(group_names == observed_combn[2, i])
    
    # vol_list <- c(observed@HVList[group_names == observed_combn[1, i]][[1]]@Volume, 
    #               observed@HVList[group_names == observed_combn[2, i]][[1]]@Volume)
    # 

    result <- lapply(observed[[c(group_1, group_2)]]@HVList, function(x) x@ValueAtRandomPoints)
    result <- do.call(cbind, result)
    # result <- result[colSums(result) != 0, ]
    result[result != 0] <- 1
    
    
    res_vol <- sum(result_combn[, i] != 0) / nrow(res) * vol_tot
    
    # get Data from observed, merge data for the two groups under comparison
    Data <- lapply(observed[[c(group_1, group_2)]]@HVList, function(x) x@Data)
    Data <- unique(do.call(rbind, Data))
    
    # volume of the significant fraction for th groups under comparison
    # res_vol <- sum(rowSums(result_combn[, i, drop = FALSE]) != 0)
    
    # ValueAtRandomPoints, obtained from result_combn
    empty_hypervolume@ValueAtRandomPoints <- result_combn[, i]
    
    # the name of the method is now n_occupancy_test
    empty_hypervolume@Method <- "n_occupancy_test"
    
    # label of the two groups under comparison
    empty_hypervolume@Name <- label_res[i]
    
    # assign merged Data
    empty_hypervolume@Data <- Data
    
    # number of dimensions
    empty_hypervolume@Dimensionality <- ncol(observed@HVList[[1]]@RandomPoints)
    
    # empty parameters
    empty_hypervolume@Parameters = list()
    
    # The coordinates of RandomPoints are the same across all the hypervolumes
    empty_hypervolume@RandomPoints = observed@HVList[[1]]@RandomPoints
    
    # assign the volume
    empty_hypervolume@Volume <- res_vol
    
    # meandensity is the same across all the hypervolumes 
    empty_hypervolume@PointDensity <- sum(result_combn[, i] != 0) / res_vol
    
    # set dimnames
    dimnames(empty_hypervolume@RandomPoints) = dn
    dimnames(empty_hypervolume@Data) = dn
    
    hv_list_res[[i]] <- empty_hypervolume
    
    
  }
  
  
  # return the results, if only 1 comparison is present return an object of class Hypervolume
  # otherwise returns a HypervolumeList
  if(ncol(observed_combn) == 1){
    result <- empty_hypervolume
  } else {
    hv_list_res <- new("HypervolumeList", HVList = hv_list_res)
    result <- hv_list_res
  }
  
  if(!exists_cluster) {
    stopCluster(cl)
    registerDoSEQ()
  }
  
  return(result)

}
