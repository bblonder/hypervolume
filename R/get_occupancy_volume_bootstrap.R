get_occupancy_volume_bootstrap <- function(path, method = "all", res_type = "raw", relative = FALSE, tol = 1e-10){

  if(! any(identical(method, "all") | identical(method, "pairwise"))){
    stop("method must be all or pairwise")
  }
  
  if(! any(identical(res_type, "raw") | identical(res_type, "summary"))){
    stop("res_type must be raw or pairwise")
  }

  file_ls <- list.files(path, pattern = ".rds")
  
  if(identical(method, "all")){
    if(identical(res_type, "raw")){
      if(relative){
        res <- lapply(file_ls, function(z) unlist(get_relative_volume(readRDS(file.path(path, z)))))
      } else {
        res <- lapply(file_ls, function(z) unlist(get_volume(readRDS(file.path(path, z)))))
      }
      res <- do.call(rbind, res)
      final_res <- res
    }
    
    if(identical(res_type, "summary")){
      if(relative){
        res <- lapply(file_ls, function(z) unlist(get_relative_volume(readRDS(file.path(path, z)))))
      } else {
        res <- lapply(file_ls, function(z) unlist(get_volume(readRDS(file.path(path, z)))))
      }
      res <- do.call(rbind, res)
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
    if(identical(res_type, "raw")){
      stop("res_type raw is currently not available for method pairwise. Set res_type to summary instead.")
    }
    
    if(identical(res_type, "summary")){
      if(relative){
        res <- lapply(file_ls, function(z) unlist(get_relative_volume(readRDS(file.path(path, z)))))
      } else {
        res <- lapply(file_ls, function(z) unlist(get_volume(readRDS(file.path(path, z)))))
      }
      res <- do.call(rbind, res)
      
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
    
  }
  

  final_res
}
