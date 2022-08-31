hypervolume_n_occupancy_bootstrap <- function(path, name = NULL, classification = NULL, method = "subsample", FUN = mean, num.points.max = NULL, verbose = TRUE,
                                              distance.factor = 1, check.hyperplane = FALSE, box_density = 5000, thin = FALSE, quant.thin = 0.5, seed = NULL){
  
  # list the number of hypervolumes for which we have bootstrap
  file_ls <- list.files(path)
  
  # extract the log file, with the information about the order of hypervolumes
  log_ls <- read.table(file.path(path, file_ls[grepl(".txt", file_ls)]), header = TRUE)
  file_ls <- file_ls[! grepl(".txt", file_ls)]
  
  match_ls <- match(log_ls[,1], file_ls)
  
  # count the number of permutations for each hypervolume
  n <- length(list.files(file.path(path,file_ls[1])))
  
  
  dir.create(file.path('./Objects', name), recursive = TRUE, showWarnings = FALSE)
  suppressWarnings(write.table(data.frame(n = numeric(), input = numeric(), re_computed = numeric(), ratio = numeric()),
              file = file.path('./Objects', name, "log_occupancy.txt"), 
              row.names = FALSE,
              append = FALSE)) 
  
  for(i in 1:n){
    hyper_perm <- lapply(file_ls, function(z) readRDS(file.path(path, z, paste("resample ", i, ".rds", sep = ""))) )
    hyper_perm <- hypervolume_join(hyper_perm[match_ls])
    hyper_perm_occupancy <- n_occupancy_function(hv_list = hyper_perm, classification = classification, method = method,
                                                 FUN = FUN, num.points.max = num.points.max, verbose = verbose,
                                                 distance.factor = distance.factor, check.hyperplane = check.hyperplane, 
                                                 box_density = box_density, print_log = TRUE, append = TRUE, iteration = i,
                                                 thin = thin, quant.thin = quant.thin,
                                                 seed = seed, path = file.path('./Objects', name, "log_occupancy.txt"))
    path_i = paste0("permutation ", i, '.rds')
    saveRDS(hyper_perm_occupancy, file.path('./Objects', name, path_i))
  }
  
  return(file.path('./Objects', name))
  
}