get_centroid_weighted <- function(hv)
{
  if (inherits(hv,"Hypervolume"))
  {
    # with the new Method n_overlap and n_overlap_test it makes sense to calculate 
    # weighted mean
    if(identical(hv@Method, "n_occupancy") | identical(hv@Method, "n_occupancy_test")){
      # get RandomPoints, no need to remove zeroes because they are removed during hypervolume_n_occupancy
      # when classification is NULL
      ran_points <- hv@RandomPoints
      apply(ran_points, 2, function(x) weighted.mean(x, hv@ValueAtRandomPoints))
    } else {
      stop("weighted centroids works for methods occupancy and occupany_test")
    }
    
    
  }
  
  else if (inherits(hv,"HypervolumeList"))
  
    { 
    method_list <- unique(unlist(lapply(hv@HVList, function(x) x@Method)))
    if(identical(method_list, "n_occupancy") | identical(method_list, "n_occupancy_test")){
      
      # get data for each element of the HypervolumeList
      
      # get RandomPoints
      ran_points <- lapply(hv@HVList, function(x) x@RandomPoints)
      
      # get ValueAtRandomPoints
      value_at_random_points <- lapply(hv@HVList, function(x) x@ValueAtRandomPoints)
      
      # get group labels
      group_name <- unlist(lapply(hv@HVList, function(x) x@Name))
      
      # initialize an empty matrix to store the results
      result <- as.data.frame(matrix(NA, ncol = ncol(ran_points[[1]]), nrow = length(group_name)))
      
      # remove zeroes from both RandomPoints and ValueAtRandomPoints and then calculate the weighted mean
      # remember, 0 are kept by hypervolume_n_occupancy because they simplify further calculations
      
      for(i in 1:length(ran_points)){
        # remove coordinates where ValueAtRandomPoints is zero
        temp_ran_points <- ran_points[[i]][value_at_random_points[[i]] > 0,]
        
        # remove ValueAtRandomPoints equal to 0
        temp_values <- value_at_random_points[[i]]
        temp_values <- temp_values[temp_values > 0]
        
        # calculate weighted means for each column of RandomPoints
        result[i, ] <- apply(temp_ran_points, 2, function(x) weighted.mean(x, temp_values))
      }
      
      # assign colnames and rownames to the result matrix
      colnames(result) <- colnames(ran_points[[1]])
      rownames(result) <- group_name
      result
      
    } else {
      stop("weighted centroids works for methods occupancy and occupancy_test")
    }
    
  }
  else
  {
    stop('Wrong class input.')
  }
}