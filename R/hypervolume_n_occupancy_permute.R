hypervolume_n_occupancy_permute <- function(name, hv_list1, hv_list2, verbose = TRUE, n = 9, cores = 1){
  
  
  ########## SOME CHECKS ##########
  
  # check if hv_list2 is of class HypervolumeList
  if(! inherits(hv_list2, "HypervolumeList")){
    stop("An object of class HypervolumeList is needed.")
  }
  
  # check if hv_list2 is of class HypervolumeList
  if(length(hv_list1@HVList) == 1){
    stop("hv_list1 must have length greater than one.")
  }
  
  # check if hv_list1 was build with hypervolume_n_occupancy
  if(hv_list1@HVList[[1]]@Method != "n_occupancy"){
    stop("hv_list1 must be calculated with hypervolume_n_occupancy.")
  }
  
  # # check if the number of elements of hv_list2 match the length of classification
  # if(length(hv_list2@HVList) != length(classification)){
  #   stop("The length of hv_list2 must be the same as the length of classification.")
  # }
  
  
  ##################################################################################################
  # retrieve parameter from hv_list1
  
  parm <- hv_list1[[1]]@Parameters
  
  classification <- parm[["classification"]]
  FUN <- parm[["FUN"]]
  distance.factor <- parm[["distance.factor"]]
  seed <- parm[["seed"]]
  num.points.max <- parm[["num.points.max"]]
  
  # set seed
  # get random state from global environment
  old_state <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  on.exit({
    # restore random state
    if (!is.null(old_state)) {
      assign(".Random.seed", old_state, envir = .GlobalEnv, inherits = FALSE)
    }
  }, add = TRUE)
  
  
  ### unique groups
  unique_groups <- unique(classification)
  
  ### check if user provided a subset of hv_list used in hypervolume_n_occupancy
  ### this could cause a problem because classification is extracted from the occupancy object
  hv_name_list1 <- sapply(hv_list1@HVList, function(x) x@Name)
  if(! setequal(hv_name_list1, classification)){
    stop("Classification in hv_list1 differs from groups classification.")
  }
  
  if(n < 2){
    stop("n must be greater than 1.")
  }
  

  

  ##################################################################################################
  # This is the same code of hypervolume_n_occupancy
  # the only things that is changed is that test points are the same used for building hv_list1
  # and are not selected randomly again.
  
  # for assuring that the function will work, not needed because we use the same RandomPoints of hv_list1
  num.points.max <-  NULL
  
  # store some properties of the hypervolumes stored in hv_list2:
  # dimensionality of random points, volume, density, names
  np_list <- unlist(lapply(hv_list2@HVList, function(x) nrow(x@RandomPoints)))
  vol_list <- unlist(lapply(hv_list2@HVList, function(x) x@Volume))
  hv_point_density_list <- np_list / vol_list
  dimhv_list <- unlist(lapply(hv_list2@HVList, function(x) ncol(x@RandomPoints)))
  hv_names_list <- unlist(lapply(hv_list2@HVList, function(x) x@Name))
  Data <- lapply(hv_list2@HVList, function(x) x@Data)
  
  if(any(is.nan(hv_point_density_list))) {
    hv_point_density_list[is.nan(hv_point_density_list)] <- NA
  }
  
  if (length(unique(dimhv_list)) > 1){
    stop('Dimensionality of hypervolumes is not the same.')
  }
  
  
  # if (check.hyperplane){
  #   hv_df <- as.data.frame(hv_list2[[i]]@Data)
  #   
  #   res <- caret::findLinearCombos(hv_df)
  #   
  #   if (res[2] != "NULL") {
  #     warning("Some data is hyperplanar")
  #   }
  # }
  
  
  
  
  # This should check if dimhv_list2 contains only as single value
  # I think it will work in most cases. It will provide strange results when,
  # for example, dimhv_list2 is a vector of length 1 containing only a 0.
  
  if(length(unique(dimhv_list)) > 1) stop("Dimensionality of hypervolumes is not the same.")
  
  # Extract the method with lapply and check if they are the same
  # I would put a stop here, but feel free to modify if needed
  
  method_check <- unlist(lapply(hv_list2@HVList, function(x) x@Method))
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
    density_list <- c(density_list, hv_point_density_list[i], num.points.max/hv_list2[[i]]@Volume)
    mindensity = min(density_list, na.rm = TRUE)
  }
  
  if (verbose) {
    cat(sprintf("Using minimum density of %f\n", mindensity))
  }
  
  
  #for each HV, keep the number of points that gives the same density to each HV according to their volume
  numpointstokeep_hv_list <- c()
  for (i in 1:length(hv_list2@HVList)){
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

  for (i in 1:length(hv_list2@HVList)){
    # randomly select points with as many points as the numpointstokeep is equal to
    hv_points_ss <- hv_list2[[i]]@RandomPoints[sample(1:np_list[i], size = numpointstokeep_hv_list[i]), , drop = FALSE]
    hv_points_ss_list[[i]] <- hv_points_ss
  }

  dim <- dimhv_list[1]

  point_density = nrow(hv_points_ss_list[[1]])/hv_list2[[1]]@Volume

  cutoff_dist = point_density^(-1/dim) * distance.factor
  
  #####################################################################################################################
  #####################################################################################################################
  if (verbose) {
    cat("Beginning ball queries... \n")
  }
  
  
  
  ############# IMPORTANT ############# 
  ### keep the same random points used in hypervolume_n_occupancy
  total_hv_points_ss <- hv_list1@HVList[[1]]@RandomPoints
  
  
  
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
    
  } 
  
  
  if (verbose){
    pb$terminate()
  }
  
  res <- do.call("cbind", final_points_intersection_list)
  # total_hv_points_ss <- total_hv_points_ss [rowSums(res) > 0, ]
  rownames(total_hv_points_ss) <- NULL
  colnames(total_hv_points_ss) <- colnames(hv_list2@HVList[[1]]@RandomPoints)
  # res <- res[rowSums(res) > 0, ]
  res[res > 0] <- 1
  rownames(res) <- NULL

  vol_list_occ <- unlist(lapply(hv_list1@HVList, function(x) x@Volume))
  res_occ <- lapply(hv_list1@HVList, function(x) x@ValueAtRandomPoints)
  res_occ <- do.call(cbind, res_occ)
  res_occ[res_occ > 0] <- 1
  
  # volume of the significant fraction for th groups under comparison
  intersection_weights <- sweep(res_occ, 1, apply(res_occ, 1, sum), "/")
  final_volume_intersection <- sum(apply(intersection_weights, 2, function(x) mean(x[x > 0], na.rm = TRUE)) * vol_list_occ, na.rm = TRUE)
  final_density <- nrow(res) / final_volume_intersection
  
  
  # basically it is the volume of the union
  # final_volume_intersection <- nrow(res) *  (nrow(final_points_intersection) / sum( final_points_intersection )) / mindensity * ncol(res)
  
  
  # get names
  # get column names assuming first set are correct for all
  cn <- dimnames(hv_list2@HVList[[1]]@RandomPoints)[[2]]
  dn <- list(NULL, cn)
  

  # Check if cluster registered to doparallel backend exists
  exists_cluster = TRUE
  if(cores > 1 & getDoParWorkers() == 1) {
    # If no cluster is registered, create a new one based on use input
    cl = makeCluster(cores)
    clusterEvalQ(cl, {
      library(hypervolume)
    })
    registerDoParallel(cl)
    exists_cluster = FALSE
  }
  
  on.exit({
    # If a cluster was created for this specific function call, close cluster and register sequential backend
    if(!exists_cluster) {
      stopCluster(cl)
      registerDoSEQ()
    }
  }, add = TRUE)
  
  
  # Create folder to store permuted hypervolumes
  dir.create(file.path('./Objects', name), recursive = TRUE, showWarnings = FALSE)

  
  # unique_groups <- unique(classification)
  
  # combinations for pairwise combinations
  pairwise_combn <- combn(unique_groups, 2)
  
  # initialize an empty vector to store the labels of pairwise combinations
  store_labels <- c()
  
  # create directories for each pairwise combination
  for(i in 1:ncol(pairwise_combn)){
    label <- paste(pairwise_combn[1, i], pairwise_combn[2, i], sep = "__")
    store_labels <- c(store_labels, label)
    dir.create(file.path('./Objects', name, label), recursive = TRUE, showWarnings = FALSE)
  }
  
  # initialize an empty hypervolume, used in the following loop
  empty_hypervolume <- new("Hypervolume")
  
  
  # the next loop is for resampling hypervolumes
  
  # basically, for each pairwise combination of the n groups (assigned with classification)
  # we resample the labels and creates a new HypervolumeList of length two
  # elements of the HypervolumeList are the resampled version of hypervolume_n_occupancy
  
  # I nested the foreach loop into a normal loop to avoid to go out of memory
  # with 16 Gb RAM I had problems in performing nested loops with foreach
  # Probably this problem is because of my bad coding ability
  
  if(verbose){
    Fun <- function(...) progress_bar_foreach(iterator = n, fun = function (a, ...)c(a, list(...)), clear = FALSE) 
  } else {
    Fun <- function(...) function (a, ...)c(a, list(...))
  }
  
  
  for(j in 1:length(store_labels)){
    
    if(verbose){
      cat("\nPermuting the pair ", store_labels[j], "...\n", sep = "")
    }
      
      
      foreach(i = 1:n, .combine = Fun()) %dopar% {
        ### extract information that will be useful later
        # create a list of length 2 to host the two hypervolumes of the i-th comparison
        # I think it can be put outside of the loop
        temp_list <- vector(2, mode = "list")
        
        # vector of TRUE FALSE to subset elements of hv_list2 belonging to
        # the two groups under comparison (remember that we are working on pairwise comparisons)
        to_subset <- classification %in% pairwise_combn[, j]
        
        # subset Data of the two groups under comparison
        data_merge <- Data[to_subset]
        
        res_temp <- res[, classification %in% pairwise_combn[, j], drop = FALSE]
        
        # resample the labels of the two groups under comparison
        classification_temp <- sample(classification[to_subset])
        
        ### first hypervolume group
        
        # extract data from resampled hypervolumes of the group 1
        
        # extract data for the first group and keep only the unique values
        data_merge_group <- data_merge[classification_temp == pairwise_combn[1, j]]
        data_merge_group <- unique(do.call(rbind, data_merge_group))
        
        # calculate volume for the first group
        # res_vol <- sum(apply(intersection_weights[, classification_temp == pairwise_combn[1, j], drop = FALSE], 2, function(x) mean(x[x > 0])) * vol_list[classification_temp == pairwise_combn[1, j]])
        
        res_vol <- apply(res_temp [, classification_temp == pairwise_combn[1, j], drop = FALSE], 1, function(x) sum(x > 0))
        res_vol <- sum(res_vol > 0) / nrow(res_temp) * final_volume_intersection
        # res_vol
        
        # assign the ValueAtRandomPoints for the first group
        empty_hypervolume@ValueAtRandomPoints <- apply(res_temp [, classification_temp == pairwise_combn[1, j], drop = FALSE], 1, FUN)
        
        # I have called the Method "n_occupancy_permute"
        empty_hypervolume@Method <- "n_occupancy_permute"
        
        # name of the first group
        empty_hypervolume@Name <- pairwise_combn[1, j]
        
        # set classification
        parm_temp <- parm
        parm_temp[["classification"]] <- classification_temp
        
        p_dens <- apply(res_temp [, classification_temp == pairwise_combn[1, j], drop = FALSE], 1, function(x) sum(x > 0))
        # do assignments for other slots
        empty_hypervolume@Dimensionality <- dim
        empty_hypervolume@Parameters = parm_temp
        empty_hypervolume@RandomPoints = matrix(total_hv_points_ss, ncol = dim)
        empty_hypervolume@Volume <- res_vol
        empty_hypervolume@Data <- data_merge_group
        empty_hypervolume@PointDensity = sum(p_dens > 0) / res_vol
        
        # set dimnames
        dimnames(empty_hypervolume@RandomPoints) = dn
        
        # assign the hypervolume of the first group to the list created at the beginning of the loops
        temp_list[[1]] <- empty_hypervolume
        
        ### second hypervolume group
        # steps are the same as for the first group
        
        data_merge_group <- data_merge[classification_temp == pairwise_combn[2, j]]
        data_merge_group <- unique(do.call(rbind, data_merge_group))
        
        res_vol <- apply(res_temp[, classification_temp == pairwise_combn[2, j], drop = FALSE], 1, function(x) sum(x > 0))
        res_vol <- sum(res_vol > 0) / nrow(res_temp) * final_volume_intersection
        
        p_dens <- apply(res_temp [, classification_temp == pairwise_combn[2, j], drop = FALSE], 1, function(x) sum(x > 0))
        empty_hypervolume@ValueAtRandomPoints <- apply(res_temp[, classification_temp == pairwise_combn[2, j], drop = FALSE], 1, FUN)
        empty_hypervolume@Method <- "n_occupancy_permute"
        empty_hypervolume@Name <- pairwise_combn[2, j]
        empty_hypervolume@Dimensionality <- dim
        empty_hypervolume@Parameters = parm_temp
        empty_hypervolume@RandomPoints = matrix(total_hv_points_ss, ncol = dim)
        empty_hypervolume@Volume <- res_vol
        empty_hypervolume@Data <-  data_merge_group
        empty_hypervolume@PointDensity =  sum(p_dens > 0) / res_vol
        
        # set dimnames
        dimnames(empty_hypervolume@RandomPoints) = dn
        
        # assign the hypervolume of the first group to the list created at the beginning of the loops
        temp_list[[2]] <- empty_hypervolume

        # transform the list in an HypervolumeList
        temp_list <- hypervolume_join(temp_list)
        name1 = paste0("permutation", as.character(i), '.rds')
        saveRDS(temp_list, file.path('./Objects', name, store_labels[j],  name1))
        
      }
    }

  return(file.path(getwd(), 'Objects', name))
}
  