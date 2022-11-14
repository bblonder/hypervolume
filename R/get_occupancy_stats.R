get_occupancy_stats <- function(hv, FUN, remove_zeroes = TRUE){

  # check if hv is an object of class Hypervolume
  # if yes, transform it to a list
  if(inherits(hv, "Hypervolume")){
    hv <- hypervolume_join(hv)
  }
  
  # check the method, this function is intended to work
  # with the occupancy routine only
  check_method <- sapply(hv@HVList, function(hv) hv@Method)
  
  if(any(! check_method %in% c("n_occupancy", "n_occupancy_permute", "n_occupancy_test"))){
    stop("get_occupancy_stats works for methods n_occupancy, n_occupancy_permute, n_occupany_test")
  }
  
  # extract the names of hypervolumes
  h_names <- sapply(hv@HVList, function(hv) hv@Name)

  # apply the function over the hypervolumes' ValueAtRandomPoints
  # remove zeroes will remove points not included in the hypervolume under evaluation
  if(remove_zeroes){
    res <- lapply(hv@HVList, function(hv) FUN(hv@ValueAtRandomPoints[hv@ValueAtRandomPoints != 0]))
  } else {
    res <- lapply(hv@HVList, function(hv) FUN(hv@ValueAtRandomPoints))
  }
  
  # format results and return them
  res <- t(do.call(rbind, res))
  colnames(res) <- h_names

  if(nrow(res) == 1){
    res <- as.vector(res)
    names(res) <- h_names
  }

  res
}
