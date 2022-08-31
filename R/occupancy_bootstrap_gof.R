occupancy_bootstrap_gof <- function(path, FUN){
  
  # get the log
  occupancy_log <- read.table(file.path(path, "log_occupancy.txt"), header = TRUE)
  unique_iterations <- unique(occupancy_log[, "n"])
  
  
  if(! is.function(FUN)){
    if(identical(FUN, "mae")){
      res_fun <- lapply(1:length(unique_iterations), function(x) mean(abs(occupancy_log[occupancy_log$n == x, "input"] - occupancy_log[occupancy_log$n == x, "re_computed"])))
      res_fun <- unlist(res_fun)
      res <- data.frame(metric = "mae", 
                 mean = mean(res_fun),
                 sd = sd(res_fun),
                 min = min(res_fun),
                 perc_2.5 = quantile(res_fun, 0.025),
                 perc_25 = quantile(res_fun, 0.25),
                 median = median(res_fun),
                 perc_75 = quantile(res_fun, 0.75),
                 perc_97.5 = quantile(res_fun, 0.975),
                 max = max(res_fun))
      rownames(res) <- NULL
    }
    
    if(identical(FUN, "rmse")){
      res_fun <- lapply(1:length(unique_iterations), function(x) sqrt(mean((occupancy_log[occupancy_log$n == x, "input"] - occupancy_log[occupancy_log$n == x, "re_computed"])^2)))
      res_fun <- unlist(res_fun)
      res <- data.frame(metric = "rmse", 
                        mean = mean(res_fun),
                        sd = sd(res_fun),
                        min = min(res_fun),
                        perc_2.5 = quantile(res_fun, 0.025),
                        perc_25 = quantile(res_fun, 0.25),
                        median = median(res_fun),
                        perc_75 = quantile(res_fun, 0.75),
                        perc_97.5 = quantile(res_fun, 0.975),
                        max = max(res_fun))
      rownames(res) <- NULL
    }
  }
  
  if(is.function(FUN)){
    res_fun <- lapply(1:length(unique_iterations), function(x) FUN(occupancy_log[occupancy_log$n == x, "input"] - occupancy_log[occupancy_log$n == x, "re_computed"]))
    res_fun <- unlist(res_fun)
    fun_name <- as.character(substitute(FUN))
    if(length(fun_name) > 1){
      fun_name = "custom"
    }
    res <- data.frame(metric = fun_name, 
                      mean = mean(res_fun),
                      sd = sd(res_fun),
                      min = min(res_fun),
                      perc_2.5 = quantile(res_fun, 0.025),
                      perc_25 = quantile(res_fun, 0.25),
                      median = median(res_fun),
                      perc_75 = quantile(res_fun, 0.75),
                      perc_97.5 = quantile(res_fun, 0.975),
                      max = max(res_fun))
    rownames(res) <- NULL
  }
  
  res
  
}


