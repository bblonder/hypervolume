hypervolume_n_occupancy <- function(hv_list, classification = NULL, method = "subsample", FUN = mean, num.points.max = NULL, verbose = TRUE,
                                    distance.factor = 1, check.hyperplane = FALSE, box_density = 5000, thin = FALSE, 
                                    quant.thin = 0.5, seed = NULL, print_log = FALSE){
  
  n_occupancy_function(hv_list = hv_list, classification = classification, method = method,
                       FUN = FUN, num.points.max = num.points.max, verbose = verbose,
                       distance.factor = distance.factor, check.hyperplane = check.hyperplane, 
                       box_density = box_density, print_log = print_log, thin = thin, 
                       quant.thin = quant.thin, seed = seed, append = FALSE)
}