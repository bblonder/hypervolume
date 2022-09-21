occupancy_to_intersection <- function(hv_list, method = "all", m = 2, tol = 1e-10){
  
  # check method
  if(! any(identical(method, "all") | identical(method, "n_wise"))){
    stop("method can be all or n_wise")
  }
  
  # unshared volume makes no sense for single hypervolumes
  if(inherits(hv_list, "Hypervolume")){
    stop("A HypervolumeList object is needed")
  }
  

  if(inherits(hv_list, "HypervolumeList")){
    
    # check construction method, it needs to be n_occupancy
    hyper_methods <- sapply(hv_list@HVList, function(x) x@Method)
    
    if(any(! hyper_methods %in% c("n_occupancy", "n_occupancy_test", "occupancy_to_union", "occupancy_to_unshared"))){
      stop("Construction method must be n_occupancy, n_occupancy_test, occupancy_to_union or occupancy_to_unshared")
    }
    
    if(identical(method, "all")){
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
      
      # extract data
      all_data <- lapply(hv_list@HVList, function(x) x@Data)
      all_data <- do.call(rbind, all_data)
      all_data <- unique(all_data)
      rownames(all_data) <- NULL
      
      # if reconstructed volume differs above a certain threshold stop the function
      
      if(any(dist(all_vol) > tol)){
        stop("The reconstruction of the total volume failed. Try to: \n
    \t 1. Increase the tolerance if results satisfy your needs. \n 
    \t 2. Check the accuracy of the results of hypervolume_n_occupancy. \n 
    Contact the maintener if none of the above applies to you.")
      }
      
      all_occupancy <- do.call(cbind, lapply(hv_list@HVList, function(z) z@ValueAtRandomPoints))
      all_occupancy[all_occupancy != 0] <- 1
      
      res <- apply(all_occupancy, 1, sum)   
      res[res < ncol(all_occupancy)] <- 0
      res[res == ncol(all_occupancy)] <- 1
      # create a temporary hypervolume to store the results of the loop
      temp_hyper <-  hv_list@HVList[[1]]
      temp_hyper@ValueAtRandomPoints <- res
      temp_hyper@Volume <- all_vol[1] * sum(res) / length(res)
      temp_hyper@Method <- "occupancy_to_intersection"
      temp_hyper@Name <- "intersection"
      temp_hyper@Data <- all_data
      temp_hyper@Parameters <- list()
      pt_d <- sum(res)/temp_hyper@Volume
      names(pt_d) <- NULL
      temp_hyper@PointDensity <- pt_d
      hv_list <- temp_hyper
      
    }
    
    if(identical(method, "n_wise")){
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
      
      # Extract occupancy values
      all_occupancy <- do.call(cbind, lapply(hv_list@HVList, function(z) z@ValueAtRandomPoints))
      all_occupancy[all_occupancy != 0] <- 1
      
      # extract the names of individual hypervolumes
      hyper_names <- sapply(hv_list@HVList, function(x) x@Name)
      
      # sanity checks on m
      if(m > ncol(all_occupancy)){
        stop("m cannot be greater than the number of hypervolume in the HypervolumeList")
      }
      
      if(m < 2){
        stop("m cannot be lower than 2")
      }
      
      
      # create combinations and the name combination
      
      combn_hyper <- combn(length(hv_list@HVList), m)
      combn_names <- combn(hyper_names, m)
      
      
      
      if(m > 2){
        comb_names <- apply(combn_names, 2, function(x) paste(abbreviate(x, minlength = 3), collapse = "_"))
      } else{
        comb_names <- apply(combn_names, 2, function(x) paste(x, collapse = "_") )
      }
      
      
      # create a temporary hypervolume to store the results of the loop
      temp_hyper <-  hv_list@HVList[[1]]
      
      # create a list to store the results of the loop
      hv_temp <- list()
      
      for(i in 1:ncol(combn_names)){
        temp_combn <- all_occupancy[, combn_hyper[, i], drop = FALSE]
        res <- apply(temp_combn, 1, sum)   
        res[res != m] <- 0
        res[res == m] <- 1
        
        combn_data <- lapply(combn_hyper[, i], function(x) hv_list@HVList[[x]]@Data)
        combn_data <- do.call(rbind, combn_data)
        combn_data <- unique(combn_data)
        rownames(combn_data) <- NULL
        
        
        temp_hyper@ValueAtRandomPoints <- res
        temp_hyper@Volume <- all_vol[1] * sum(res) / length(res)
        temp_hyper@Method <- "occupancy_to_intersection"
        temp_hyper@Name <- comb_names[i]
        temp_hyper@Data <- combn_data
        temp_hyper@Parameters <- list()
        pt_d <- sum(res)/temp_hyper@Volume
        names(pt_d) <- NULL
        temp_hyper@PointDensity <- pt_d
        hv_temp[[i]] <- temp_hyper
      }
      
      hv_list <- hypervolume_join(hv_temp)
      if(length(hv_list@HVList) == 1){
        hv_list <- hv_list@HVList[[1]]
      }
      
    }
  }
  
  hv_list
}