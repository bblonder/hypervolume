occupancy_to_unshared <- function(hv_list, method = "all", tol = 1e-10){
  
  # check method
  if(! any(identical(method, "all") | identical(method, "pairwise"))){
    stop("method can be all or pairwise")
  }
  
  # unshared volume makes no sense for single hypervolumes
  if(inherits(hv_list, "Hypervolume")){
    stop("A HypervolumeList object is needed")
  }
  

  if(inherits(hv_list, "HypervolumeList")){
    
    # check construction method, it needs to be n_occupancy
    hyper_methods <- sapply(hv_list@HVList, function(x) x@Method)
    
    if(any(! hyper_methods %in% c("n_occupancy", "n_occupancy_test", "occupancy_to_union", "occupancy_to_intersection"))){
      stop("Construction method must be n_occupancy, n_occupancy_test, occupancy_to_union or occupancy_to_intersection")
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
      
      # names
      all_names <- sapply(hv_list@HVList, function(x) x@Name)
      
      # if reconstructed volume differs above a certain threshold stop the function
      
      
      if(any(dist(all_vol) > tol)){
        stop("The reconstruction of the total volume failed. Try to: \n
    \t 1. Increase the tolerance if results satisfy your needs. \n 
    \t 2. Check the accuracy of the results of hypervolume_n_occupancy. \n 
    Contact the maintener if none of the above applies to you.")
      }
      
      all_occupancy <- do.call(cbind, lapply(hv_list@HVList, function(z) z@ValueAtRandomPoints))
      all_occupancy[all_occupancy != 0] <- 1
      
      
      temp_hyper <- hv_list[[1]]
      hv_temp <- list()
      for(i in 1:length(hv_list@HVList)){
        res_temp <- all_occupancy[, i] - apply(all_occupancy[, -i, drop = FALSE], 1, function(z) sum(any(z > 0)))
        res_temp[res_temp < 1] <- 0
        
        temp_hyper@ValueAtRandomPoints <- res_temp
        temp_hyper@Volume <-  all_vol[[1]] / all_length * sum(res_temp)
        temp_hyper@Method <- "occupancy_to_unshared"
        temp_hyper@Name <- all_names[i]
        temp_hyper@Data <- hv_list[[i]]@Data
        temp_hyper@Parameters <- list()
        pt_d <- sum(abs(res_temp))/temp_hyper@Volume
        names(pt_d) <- NULL
        temp_hyper@PointDensity <- pt_d
        hv_temp[[i]] <- temp_hyper
      }
      
      hv_list <- hypervolume_join(hv_temp)
      if(length(hv_list@HVList) == 1){
        hv_list <- hv_list@HVList[[1]]
      }
      
    }
    
    if(identical(method, "pairwise")){
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
      
      # create combinations and the name combination
      combn_hyper <- combn(length(hv_list@HVList), 2)
      combn_names <- combn(hyper_names, 2)
      
      # create a temporary hypervolume to store the results of the loop
      temp_hyper <-  hv_list@HVList[[1]]
      
      # create a list to store the results of the loop
      hv_temp <- list()
      
      for(i in 1:ncol(combn_names)){
        a <- combn_hyper[1, i] 
        b <- combn_hyper[2, i]
        res_temp <- all_occupancy[, a] - all_occupancy[, b, drop = FALSE]
        
        combn_data <- lapply(combn_hyper[, i], function(x) hv_list@HVList[[x]]@Data)
        combn_data <- do.call(rbind, combn_data)
        combn_data <- unique(combn_data)
        rownames(combn_data) <- NULL
        
        temp_hyper@ValueAtRandomPoints <- res_temp[,1]
        temp_hyper@Volume <-  all_vol[[1]] / all_length * sum(res_temp != 0)
        temp_hyper@Name <- paste(combn_names[, i], collapse = "_")
        temp_hyper@Method <- "occupancy_to_unshared"
        temp_hyper@Data <- combn_data
        temp_hyper@Parameters <- list()
        pt_d <- sum(abs(res_temp))/temp_hyper@Volume
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