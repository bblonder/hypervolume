hypervolume_n_occupancy <- function(hv_list, classification = NULL, FUN = mean, num.points.max = NULL, verbose = TRUE, distance.factor = 1, check.hyperplane = FALSE){
  
  # check if hv_list is of class HypervolumeList
  if(! class(hv_list) %in% "HypervolumeList"){
    stop("An object of class HypervolumeList is needed.")
  }
  
  if(! is.null(classification)){
    if(length(hv_list@HVList) != length(classification)){
      stop("The length of hv_list must be the same as the length of classification.")
    }
  }

  

  # store some properties of the hypervolumes stored in hv_list:
  # dimensionality of random points, volume, density, names
  np_list <- unlist(lapply(hv_list@HVList, function(x) nrow(x@RandomPoints)))
  vol_list <- unlist(lapply(hv_list@HVList, function(x) x@Volume))
  hv_point_density_list <- np_list / vol_list
  dimhv_list <- unlist(lapply(hv_list@HVList, function(x) ncol(x@RandomPoints)))
  hv_names_list <- unlist(lapply(hv_list@HVList, function(x) x@Name))
  Data <- lapply(hv_list@HVList, function(x) x@Data)
  
  if(any(is.nan(hv_point_density_list))) {
    hv_point_density_list[is.nan(hv_point_density_list)] <- NA
  }
  
  if (length(unique(dimhv_list)) > 1){
    stop('Dimensionality of hypervolumes is not the same.')
  }
  
  
  if (check.hyperplane){
    hv_df <- as.data.frame(hv_list[[i]]@Data)
    
    res <- findLinearCombos(hv_df)
    
    if (res[2] != "NULL") {
      warning("Some data is hyperplanar")
    }
  }

  
  
  
  # This should check if dimhv_list contains only as single value
  # I think it will work in most cases. It will provide strange results when,
  # for example, dimhv_list is a vector of length 1 containing only a 0.
  
  if(length(unique(dimhv_list)) > 1) stop("Dimensionality of hypervolumes is not the same.")
  
  # Extract the method with lapply and check if they are the same
  # I would put a stop here, but feel free to modify if needed
  
  method_check <- unlist(lapply(hv_list@HVList, function(x) x@Method))
  if(length(unique(method_check)) > 1) stop("Hypervolumes building method is not the same.")
  
  
  #calculate max number of points according to dimensionality 
  if (is.null(num.points.max)) {
      num.points.max = ceiling(10^(3 + sqrt(unique(dimhv_list))))
      if (verbose) {
        cat(sprintf("Choosing num.points.max=%.0f (use a larger value for more accuracy.)\n", 
                    num.points.max))
    }
  }

  
  #get the density of the HV with minimum points density
  density_list <- c()
  for (i in 1:length(hv_point_density_list)){
    density_list <- c(density_list, hv_point_density_list[i], num.points.max/hv_list[[i]]@Volume)
    mindensity = min(density_list, na.rm = TRUE)
  }
  
  if (verbose) {
    cat(sprintf("Using minimum density of %f\n", mindensity))
  }
  
  
  #for each HV, keep the number of points that gives the same density to each HV according to their volume
  numpointstokeep_hv_list <- c()
  for (i in 1:length(hv_list@HVList)){
    numpointstokeep_hv <- floor(mindensity*vol_list[i])
    numpointstokeep_hv_list <- c(numpointstokeep_hv_list, numpointstokeep_hv)
  }
  
  
  ####################################################################################################################
  ###   check that no HV is  HV are empty
  ####################################################################################################################
  
  is_any_pointstokeep_null = FALSE
  nopointstokeep <- c()
  
  for (i in 1:length(numpointstokeep_hv_list)){
    if (numpointstokeep_hv_list[i] == 0 | is.null(numpointstokeep_hv_list[i])){
      is_any_pointstokeep_null = TRUE
      nopointstokeep <- c(nopointstokeep, i)
    }
  }
  
  if (is_any_pointstokeep_null == TRUE){
    stop(paste0("hv",nopointstokeep,"has no random points and is empty."))
    
  }
  
  #####################################################################################################################  
  #####################################################################################################################
  
  
  #subsample all hypervolumes to get same points density
  hv_points_ss_list <- list()
  
  for (i in 1:length(hv_list@HVList)){
    # randomly select points with as many points as the numpointstokeep is equal to
    hv_points_ss <- hv_list[[i]]@RandomPoints[sample(1:np_list[i], size = numpointstokeep_hv_list[i]), , drop = FALSE] 
    hv_points_ss_list[[i]] <- hv_points_ss
  }
  
  dim <- dimhv_list[1]
  
  point_density = nrow(hv_points_ss_list[[1]])/hv_list[[1]]@Volume
  
  cutoff_dist = point_density^(-1/dim) * distance.factor
  
  #####################################################################################################################
  #####################################################################################################################
  if (verbose) {
    cat("Beginning ball queries... \n")
  }
  
  
  
  # determine the range to build the convex hull
  hv_points_ss_list <- lapply(hv_list@HVList, function(x) x@RandomPoints)
  total_hv_points_ss <- do.call("rbind", hv_points_ss_list)
  
  
  #compare the set to the resampled HVs, individually 
  
  if (verbose){
    pb <- progress_bar$new(total = length(hv_points_ss_list))
    pb$tick(0)
  }
  
  final_points_intersection_list <- vector(mode = "list", length(hv_points_ss_list))
  
  for (i in 1:length(hv_points_ss_list)){
    
    if (verbose){
      if (!pb$finished){
        pb$update(i/length(hv_points_ss_list))
      }
    }
    
    final_points_intersection_list[[i]] <- evalfspherical(data = hv_points_ss_list[[i]], radius = cutoff_dist, 
                                                                        points =  total_hv_points_ss, verbose = verbose )
    
    # hv_points_in_i = as.data.frame(total_hv_points_ss)[hv_points_in_i_all > 0, 
    #                                                    , drop = FALSE]
    # 
    # final_points_intersection_list[[i]] <- hv_points_in_i
  } 
  
  if (verbose){
    pb$terminate()
  }
  
  res <- do.call("cbind", final_points_intersection_list)
  

  
  total_hv_points_ss <- total_hv_points_ss [rowSums(res) > 0, ]
  
  ### START IMPORTANT !!!!!!!
  # weight each random points as the inverse of their weight. This is intended for assuring that more 
  # dense regions will not be oversampled. It will affect random sampling of points at line 196
  
  weight <- 1 / apply(res, 1, sum)
  
  ### END IMPORTANT
  
  
  rownames(total_hv_points_ss) <- NULL
  colnames(total_hv_points_ss) <- colnames(hv_list@HVList[[1]]@RandomPoints)
  weight <- weight[rowSums(res) > 0]
  res <- res[rowSums(res) > 0, ]
  res[res > 0]<- 1
  rownames(res) <- NULL
  
  
  # resample the points dividing their number by the number of hypervolumes compared
  num_points_to_sample_in_intersection = nrow(res) / length(unique(classification)) #### IMPORTANT, DIVED BY THE NUMBER OF GROUPS
  to_keep <- sample(1:nrow(res), size = num_points_to_sample_in_intersection, prob = weight)
  total_hv_points_ss <- total_hv_points_ss[to_keep, , drop = FALSE]
  final_points_intersection <- res[to_keep, , drop = FALSE]
  
  
  # basically it is the volume of the union
  final_volume_intersection <- nrow(final_points_intersection) / mindensity 
  
  
  # get names
  # get column names assuming first set are correct for all
  cn <- dimnames(hv_list@HVList[[1]]@RandomPoints)[[2]]
  dn <- list(NULL, cn)
  
  
  if(is.null(classification)){
    
    Data <- unique(do.call("rbind", Data))
    #generate a new HV corresponding to the overlap between all the input HVs
    result = new("Hypervolume")
    result@Method = "n_occupancy"
    result@Data = Data
    result@Dimensionality = dim
    result@Volume = final_volume_intersection
    result@PointDensity = mindensity
    result@Parameters = list()
    result@RandomPoints = matrix(total_hv_points_ss, ncol = dim)
    result@ValueAtRandomPoints = apply(final_points_intersection, 1, sum) / ncol(final_points_intersection)
    
    # set dimnames
    dimnames(result@RandomPoints) = dn
    dimnames(result@Data) = dn
    
  }
  
  if(!is.null(classification)){
    
    # split hypervolumes by classification and then apply the aggregation function
    unique_groups <- unique(classification)
    #hv_list_split <- split(hv_list@HVList, classification)
    
    # create an empty list to store the results 
    hv_list_res <- vector(length(unique_groups), mode = "list")
    
    # create an empty hypervolume as a model for storing results
    empty_hypervolume <- new("Hypervolume")
    
    for(i in 1:length(unique_groups)){
      
      # extract final points intersection of the i-th group for further calculation
      result <- final_points_intersection[, classification == unique_groups[i], drop = FALSE]
      
      # merge data points of the hypervolumes of the i-th group
      # unique because some data points can be shared among the hypervolumes
      data_merged <- unique(do.call(rbind, Data[classification == unique_groups[i]]))
      
      # calculate the hypervolume as the number of points or results that are greater than 0 times mindsensity
      res_vol <-  nrow(result[rowSums(result) > 0, ]) / nrow(result) * final_volume_intersection
      
      # apply FUN, mean as the default, to result. This calculates the occupancy aka how many hypervolumes include
      # a given random point
      empty_hypervolume@ValueAtRandomPoints <- apply(result, 1, FUN)
      
      # assign a method      
      empty_hypervolume@Method <- "n_occupancy"
      
      # assign the name of the i-th to the hypervolume
      empty_hypervolume@Name <- unique_groups[i]
      
      # assign data points as did at line 241
      empty_hypervolume@Data <-  data_merged
      empty_hypervolume@Dimensionality <- dim
      empty_hypervolume@Parameters = list()
      
      # Random points are common to each group. Random points whose ValueAtRandomPoints is 0 are not
      # removed. This is for assuring the correct behaviour of furter calculations.
      
      empty_hypervolume@RandomPoints = matrix(as.matrix(as.data.frame(total_hv_points_ss)), 
                                   ncol = dim)
      empty_hypervolume@Volume <- res_vol
      empty_hypervolume@PointDensity = mindensity
      
      # set dimnames
      dimnames(empty_hypervolume@RandomPoints) = dn
      dimnames(empty_hypervolume@Data) = dn
      
      hv_list_res[[i]] <- empty_hypervolume
      

      
      
    }
    
    
    # if classification is NULL an object of class Hypervolume is created
    # otherwise an object of class HypervolumeList is returned
    
    if(length(unique_groups) == 1){
      result <- empty_hypervolume
    } else {
      hv_list_res <- new("HypervolumeList", HVList = hv_list_res)
      result <- hv_list_res
    }
    
  }
  
  result
}