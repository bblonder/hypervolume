get_relative_volume <- function(hv_list, tol = 1e-10){
  
  # check if the object is an HypervolumeList
  # calculation makes no sense for a single hypervolume
  
  if(! inherits(hv_list, "HypervolumeList")){
    stop("hv_list must be an HypervolumeList")
  }

  if(length(hv_list@HVList) < 2){
    stop("At least two hypervolumes are needed for a relative volume")
  }
  
  # check if hypervolumes method were built within the occupancy routine
  hyper_methods <- sapply(hv_list@HVList, function(x) x@Method)

  if(any(! hyper_methods %in% c("n_occupancy", "n_occupancy_test", "n_occupancy_permute",
                                "occupancy_to_union", "occupancy_to_ushared", "occupancy_to_intersection"))){
    stop("Construction method must be one of n_occupancy, n_occupancy_test, n_occupancy_permute,
         occupancy_to_union, occupancy_to_unshared or occupancy_to_intersection")
  }


  # get volume of the hyper list and try to reconstruct the volume of the union
  # the reconstruction is made for each hypervolume in the hypervolume list
  vols <- get_volume(hv_list)
  #all_length <- length(hv_list@HVList[[1]]@ValueAtRandomPoints)
  combine_vrp <- lapply(hv_list@HVList, function(x) x@ValueAtRandomPoints)
  combine_vrp <- do.call(cbind, combine_vrp)
  zeroes_check <- apply(combine_vrp, 1, function(x) sum(abs(x)))
  all_length <- sum(zeroes_check != 0)
  i_length <- unlist(lapply(hv_list@HVList, function(z) sum(z@ValueAtRandomPoints != 0)))
  all_vol <- vols*all_length/i_length

  # if reconstructed volume differs above a certain threshold stop the function
  
  if(any(dist(all_vol) > tol)){
    stop("The reconstruction of the total volume failed. Try to: \n
    \t 1. Increase the tolerance if results satisfy your needs. \n 
    \t 2. Check the accuracy of the results of hypervolume_n_occupancy. \n 
    Contact the maintener if none of the above applies to you.")
  }

  # # get the relative volume
  # if(identical(method, "entire")){
  #   res <- vols / all_vol[[1]]
  # }

  # # get the relative volume of
  # if(identical(method, "unique_fraction")){
  #   all_occupancy <- do.call(cbind, lapply(hv_list@HVList, function(z) z@ValueAtRandomPoints))
  #   all_occupancy[all_occupancy != 0] <- 1
  #   res <- rep(NA, length(vols))
  # 
  #   for(i in 1:length(vols)){
  #     res_temp <- all_occupancy[, i] - apply(all_occupancy[, -i, drop = FALSE], 1, function(z) sum(any(z > 0)))
  #     res_temp <- sum(res_temp[res_temp > 0])*all_vol[[1]]/all_length
  #     res[i] <- res_temp / all_vol[[1]]
  #   }
  # 
  #   names(res) <- names(vols)
  # 
  # }

  res <- vols / all_vol[[1]]
  res
}
