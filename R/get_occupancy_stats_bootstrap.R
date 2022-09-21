get_occupancy_stats_bootstrap <- function(path, FUN, remove_zeroes = TRUE, method = "pairwise", res_type = "summary", verbose = TRUE, cores = 1){
  
  # check if the path exists
  if(!dir.exists(path)){
    stop("This path does not exist.")
  }
  
  if(identical(res_type, "raw") & identical(method, "pairwise")){
    stop("The res_type raw is not available when method is equal to pairwise")
  }
  
  
  # initialize multi-core calculations
  
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
  
  
  file_list <- list.files(path, pattern = ".rds")
  iterations <- length(file_list)
  
  if(verbose){
    cat("Begin stats calculation...")
    cat("\n")
  }
  
  col_n <- sapply(readRDS(file.path(path, file_list[1]))@HVList, function(x) x@Name) 
  
  if(verbose){
    result <- suppressWarnings(foreach(i = 1:iterations, .combine = progress_bar_foreach(iterations, fun = rbind)) %dopar%
      {
        h1 <- readRDS(file.path(path, file_list[i]))
        res <- sapply(h1@HVList, function(x) FUN(x@ValueAtRandomPoints[x@ValueAtRandomPoints > 0]))
        return(res)
      }) 
    
    
  } else {
    result <- suppressWarnings(foreach(i = 1:iterations, .combine = rbind) %dopar%
      {
        h1 <- readRDS(file.path(path, file_list[i]))
        res <- sapply(h1@HVList, function(x) FUN(x@ValueAtRandomPoints[x@ValueAtRandomPoints > 0]))
        return(res)
      }) 
  }

  if (verbose){
    cat("Stats calculation completed...\n")
  }

  
  rownames(result) <- NULL
  colnames(result) <- col_n
  
  if(identical(method, "raw")){
    final_res <- result
  }
  
  
  
  if(identical(method, "all")){
    
    res <- result
    
    if(identical(res_type, "raw")){
      final_res <- res
      rownames(final_res) <- NULL
    }
    
    if(identical(res_type, "summary")){
      res <- data.frame(melt(data.table(res), measure.vars = colnames(res)))
      final_res <- do.call(data.frame, aggregate(value ~ variable, res,
                                                 FUN = function(x) c(mean = mean(x), sd = sd(x), min = min(x),
                                                                     quantile_2.5 = quantile(x, 0.025),
                                                                     median = quantile(x, 0.5),
                                                                     quantile_97.5 = quantile(x, 0.975),
                                                                     max = max(x),
                                                                     skewness = skewness(x),
                                                                     kurtosis = kurtosis(x))))
      colnames(final_res) <- c("hypervolume", "mean", "sd", "min", "quantile_2.5", "median", "quantile_97.5", "max",
                               "skewness", "kurtosis")
    }
  }
  
  
  
  if(identical(method, "pairwise")){
    
    res <- result
    
    file_combn <- combn(1:ncol(res), 2)
    char_combn <- combn(colnames(res), 2, FUN = function(x) paste(x[1], x[2], sep = " - "))
    res_pairwise <- matrix(NA, ncol = 9, nrow = ncol(file_combn))
    colnames(res_pairwise) <- c("mean", "sd", "min", "quantile_2.5", "median", "quantile_97.5", "max",
                                "skewness", "kurtosis")
    
    
    
    for(i in 1:ncol(file_combn)){
      res_temp <- apply(res[, file_combn[, i]], 1, function(x) x[1] - x[2])
      res_pairwise[i, 1] <-  mean(res_temp)
      res_pairwise[i, 2] <-  sd(res_temp)
      res_pairwise[i, 3] <-  min(res_temp)
      res_pairwise[i, 4] <-  quantile(res_temp, c(0.025))
      res_pairwise[i, 5] <-  quantile(res_temp, c(0.5))
      res_pairwise[i, 6] <-  quantile(res_temp, c(0.975))
      res_pairwise[i, 7] <-  max(res_temp)
      res_pairwise[i, 8] <-  skewness(res_temp)
      res_pairwise[i, 9] <-  kurtosis(res_temp)
      
    }
    final_res <- data.frame(comparison = char_combn, res_pairwise)
  }
  
  final_res
}
