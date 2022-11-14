occupancy_filter <- function(hv, operator = NULL, filter = NULL, tol = 1e-10){

  if(inherits(hv, "Hypervolume")){
    if(! any(hv@Method %in% c("n_occupancy", "n_occupancy_test", "occupancy_to_unshared")) ){
      stop("Construction method must be one of n_occupancy or n_occupancy_test")
    }
    i_vol <- get_volume(hv)
    i_vol_length <- sum(hv@ValueAtRandomPoints != 0)
    occ_temp <- hv@ValueAtRandomPoints
    occ_temp[eval(call(operator, occ_temp, filter))] <- 0 
    hv@ValueAtRandomPoints <- occ_temp
    hv@Volume <- i_vol / i_vol_length * sum(occ_temp != 0)
  }
  
  if(inherits(hv, "HypervolumeList")){
    hyper_methods <- sapply(hv@HVList, function(x) x@Method)
    if(any(! hyper_methods %in% c("n_occupancy", "n_occupancy_test", "occupancy_to_unshared"))){
      stop("Construction method must be one of n_occupancy or n_occupancy_test")
    }
    
    # get volume of the hyper list and try to reconstruct the volume of the union
    # the reconstruction is made for each hypervolume in the hypervolume list
    vols <- get_volume(hv)
    all_length <- lapply(hv@HVList, function(x)x@ValueAtRandomPoints != 0)
    all_length <- do.call(cbind, all_length)
    all_length <- apply(all_length, 1, sum)
    all_length <- sum(all_length != 0)
    i_length <- sapply(hv@HVList, function(z) sum(z@ValueAtRandomPoints != 0))
    all_vol <- vols*all_length/i_length
    
    # if reconstructed volume differs above a certain threshold stop the function
    
    if(any(dist(all_vol) > tol)){
      stop("The reconstruction of the total volume failed. Try to: \n
    \t 1. Increase the tolerance if results satisfy your needs. \n 
    \t 2. Check the accuracy of the results of hypervolume_n_occupancy. \n 
    Contact the maintener if none of the above applies to you.")
    }
    
    for( i in 1:length(hv@HVList)){
      occ_temp <- hv@HVList[[i]]@ValueAtRandomPoints
      occ_temp[eval(call(operator, occ_temp, filter))] <- 0 
      hv@HVList[[i]]@ValueAtRandomPoints <- occ_temp
      hv@HVList[[i]]@Volume <-  all_vol[[1]] / all_length * sum(occ_temp != 0)
    }
  }

  hv
  
}