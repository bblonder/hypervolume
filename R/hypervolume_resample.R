hypervolume_resample <- function(name, hv, method, n = 10, points_per_resample = 'sample_size', seq = 3:nrow(hv@Data), k = 5, cores = 1, verbose = TRUE, to_file = TRUE, mu = NULL, sigma = NULL, cols_to_weigh = 1:ncol(hv@Data), weight_func = NULL) {
  # Check for valid inputs
  if (n <= 0) {
    stop("Invalid value for n")
  } else if (points_per_resample != "sample size" & points_per_resample <= 0) {
    stop("Invalid value for points_per_resample")
  } else if (seq[1] <= 0) {
    stop("Invalid input for seq")
  } else if (method == 'weighted bootstrap' & is.null(weight_func) & (length(mu) != length(sigma) | length(mu) != length(cols_to_weigh))) {
    stop("mu, sigma, and cols_to_weigh must have same length")
  } else if (cores < 1) {
    stop("cores must be greater than or equal to 1")
  }
  
  # Create Objects folder in current working directory if one doesn't already exist, then call functions
  dir.create('./Objects', showWarnings = FALSE)
  if (method == 'bootstrap') {
    return(bootstrap(name, hv, n, points_per_resample, cores, verbose, to_file))
  } else if (method == 'bootstrap seq') {
    return(bootstrap_seq(name, hv, n, seq, cores, verbose))
  } else if (method == 'weighted bootstrap') {
    return(weighted_bootstrap(name, hv, n, points_per_resample, cores, verbose, to_file, mu, sigma, cols_to_weigh, weight_func))
  } else if (method == "k_split") {
    return(k_split(name, hv, k, cores, verbose))
  }
}