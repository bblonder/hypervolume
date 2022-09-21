get_occupancy_intersection_bootstrap <- function(path, method = "n_wise", res_type = "summary", m = 2, relative = FALSE,  tol = 1e-10){
  
  if(! any(identical(res_type, "raw") | identical(res_type, "summary"))){
    stop("res_type must be raw or summary")
  }
  
  if(! any(identical(method, "all") | identical(method, "n_wise"))){
    stop("method must be all or n_wise")
  }
  
  
  file_ls <- list.files(path, pattern = ".rds")
  
  if(identical(method, "all")){
    
    res <- lapply(file_ls, function(z) occupancy_to_intersection(readRDS(file.path(path, z)), method = "all", tol = tol))
    
    if(relative){
      res <- lapply(res, get_volume)
      res_union <- lapply(file_ls, function(x)get_volume(occupancy_to_union(readRDS(file.path(path, x)), tol = tol)))
      res <- do.call(rbind, res)
      res_union <- do.call(c, res_union)
      res <- sweep(res, 1, res_union, "/")
    } else {
      res <- lapply(res, get_volume)
      res <- do.call(rbind, res)
    }

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
      final_res[, "hypervolume"] <- as.character(final_res[, "hypervolume"])
    }
  }
  

  if(identical(method, "n_wise")){
    res <- lapply(file_ls, function(z) occupancy_to_intersection(readRDS(file.path(path, z)), method = "n_wise", m = m, tol = tol))
    
    if(relative){
      res <- lapply(res, get_volume)
      res_union <- lapply(file_ls, function(x)get_volume(occupancy_to_union(readRDS(file.path(path, x)), method = "n_wise", m = m, tol = tol)))
      res <- do.call(rbind, res)
      res_union <- do.call(rbind, res_union)
      res <- res / res_union
    } else {
      res <- lapply(res, get_volume)
      res <- do.call(rbind, res)
    }
    
    
    
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
      colnames(final_res) <- c("comparison", "mean", "sd", "min", "quantile_2.5", "median", "quantile_97.5", "max",
                               "skewness", "kurtosis")
      final_res[, "comparison"] <- gsub("_", " - ", final_res[, "comparison"])
    }
  }
  final_res
}
