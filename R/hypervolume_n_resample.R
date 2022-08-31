hypervolume_n_resample <- function(name, hv_list, n = 10, points_per_resample = 'sample_size', cores = 1, verbose = TRUE, seed = NULL){
  # Check for valid inputs
  if (n <= 0) {
    stop("Invalid value for n")
  } else if (points_per_resample != "sample size" & points_per_resample <= 0) {
    stop("Invalid value for points_per_resample")
  } else if (cores < 1) {
    stop("cores must be greater than or equal to 1")
  } else if (inherits(hv_list, "Hypervolume")) {
    stop("For one hypervolume use hypervolume_resample")
  }
  
  bootstrap_n_hypervolumes(name = name, hv = hv_list, n = n, points_per_resample = points_per_resample, cores = cores, verbose = verbose, seed = seed)
  
}