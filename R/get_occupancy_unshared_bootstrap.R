get_occupancy_unshared_bootstrap <- function(path, method = "pairwise", res_type = "summary", relative = FALSE,  tol = 1e-10){
  
  if(! any(identical(res_type, "raw") | identical(res_type, "summary"))){
    stop("res_type must be raw or summary")
  }
  
  if(! any(identical(method, "all") | identical(method, "pairwise"))){
    stop("method must be all or pairwise")
  }
  
  file_ls <- list.files(path, pattern = ".rds")
  
  if(identical(method, "all")){
    
    res <- lapply(file_ls, function(z) occupancy_to_unshared(readRDS(file.path(path, z)), method = "all", tol = tol))
    
    if(relative){
      res <- lapply(res, get_volume)
      res_union <- lapply(file_ls, function(x)get_volume(occupancy_to_union(readRDS(file.path(path, x)), method = "all", tol = tol)))
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
  
  if(identical(method, "pairwise")){
    res <- lapply(file_ls, function(z) occupancy_to_unshared(readRDS(file.path(path, z)), method = "all", tol = tol))
    res <- do.call(rbind, lapply(res, function(x) get_volume(x)))
    
    combn_hyper <-  readRDS(file.path(path, file_ls[1]))
    combn_names <- sapply(combn_hyper@HVList, function(x) x@Name)
    combn_names <- combn(combn_names, 2)
    
    combn_res <- matrix(NA, ncol = ncol(combn_names), nrow = length(file_ls))
    colnames(combn_res) <- apply(combn_names, 2, function(x) paste(x, collapse = "_"))
    
    for(i in 1:ncol(combn_names)){
      combn_res[, i] <- res[, combn_names[1,i]] - res[, combn_names[2,i]]
    }
    
    res <- combn_res
    
    if(relative){
      res_union <- lapply(file_ls, function(z) get_volume(occupancy_to_union(readRDS(file.path(path, z)), method = "n_wise", m = 2, tol = tol)))
      res_union <- do.call(rbind, res_union)
      res_union <- res_union[,match(colnames(res_union), colnames(res))]
      res <- res/res_union
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
