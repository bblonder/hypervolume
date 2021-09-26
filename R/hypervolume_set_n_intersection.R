hypervolume_set_n_intersection <- function (hv_list, num.points.max = NULL, verbose = TRUE,
                                            distance.factor = 1, check.hyperplane = FALSE) {

  ########################################################################################################################
  np_list <- c()
  hv_point_density_list <- c()
  dimhv_list <- c()
  hv_names_list <- c()

  #run through all the hypervolume and for each one store the data in the lists
  for (i in 1:length(hv_list@HVList)){

    hv_names_list <- c(hv_names_list, hv_list[[i]]@Name)
    np <- nrow(hv_list[[i]]@RandomPoints)          #get number of random points
    hv_point_density <- np/hv_list[[i]]@Volume     #calculate point density
    dimhv = ncol(hv_list[[i]]@RandomPoints)        #get dimensionality

    if (is.nan(hv_point_density)) {
      hv_point_density <- NA
    }

    #add the data to the lists
    np_list <- c(np_list, np)
    hv_point_density_list <- c(hv_point_density_list, hv_point_density)
    dimhv_list <- c(dimhv_list, dimhv)

    if (check.hyperplane == TRUE){
      hv_df <- as.data.frame(hv_list[[i]]@Data)

      res <- caret::findLinearCombos(hv_df)

      if (res[2] != "NULL") {
        warning("Some data is hyperplanar")
      }
    }
  }

  #check that each dimensionality is the same as for the 1st HV
  if(length(unique(dimhv_list)) > 1) stop("Dimensionality of hypervolumes is not the same.")

  #check that the building method is unique
  method_check <- unlist(lapply(hv_list@HVList, function(x) x@Method))
  if(length(unique(method_check)) > 1) stop("Hypervolumes building method is not the same.")

  dim = dimhv_list[1]   #dim is the unique dimensionality

  #calculate max number of points according to dimensionality
  if (is.null(num.points.max)) {
    num.points.max = ceiling(10^(3 + sqrt(hv_list[[i]]@Dimensionality)))
    if (verbose) {
      cat(sprintf("Choosing num.points.max=%.0f (use a larger value for more accuracy.)\n",
                  num.points.max))
    }
  }

  #get the density of the HV with minimum points density
  density_list <- c()
  for (i in 1:length(hv_point_density_list)){
    density_list <- c(density_list, hv_point_density_list[i], num.points.max/hv_list[[i]]@Volume)
    mindensity = min(density_list,
                     na.rm = T)
  }

  if (verbose) {
    cat(sprintf("Using minimum density of %f\n", mindensity))
  }

  #for each HV, keep the number of points that gives the same density to each HV according to their volume
  numpointstokeep_hv_list <- c()
  for (i in 1:length(hv_list@HVList)){
    numpointstokeep_hv <- floor(mindensity*hv_list[[i]]@Volume)
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
    hv_points_ss = hv_list[[i]]@RandomPoints[sample(1:np_list[i], size = numpointstokeep_hv_list[i]),
                                             , drop = FALSE] #randomly select points with as many points as the numpointstokeep is equal to
    hv_points_ss_list[[i]] <- hv_points_ss
  }

  point_density = nrow(hv_points_ss_list[[1]])/hv_list[[1]]@Volume

  cutoff_dist = point_density^(-1/dim) * distance.factor

  #####################################################################################################################
  #####################################################################################################################
  if (verbose) {
    cat("Beginning ball queries... \n")
  }

  #get all the points of the HVs in a single set
  total_hv_points_ss <- c()
  for (i in 1:length(hv_points_ss_list)){
    total_hv_points_ss <- rbind (total_hv_points_ss, as.data.frame(hv_points_ss_list[[i]]))
  }
  total_hv_points_ss <-  as.matrix(total_hv_points_ss )


  #intersection_hv_points <-c()
  #volume_intersection_list <- c()
  final_points_intersection_list <- c()

  #compare the set to the resampled HVs, individually

  if (verbose){
    pb <- progress_bar$new(total = length(hv_points_ss_list))
    pb$tick(0)
  }

  for (i in 1:length(hv_points_ss_list)){

    if (verbose){
      if (!pb$finished == TRUE){
        pb$update(i/length(hv_points_ss_list))
      }
    }

    hv_points_in_i_all <- evalfspherical(data = hv_points_ss_list[[i]], radius = cutoff_dist,
                                         points =  total_hv_points_ss, verbose = verbose )

    hv_points_in_i = as.data.frame(total_hv_points_ss)[hv_points_in_i_all > 0,
                                                       , drop = FALSE]

    final_points_intersection_list[[i]] <- hv_points_in_i
  }

  if (verbose){
    pb$terminate()
  }

  #keep only the points that are common for all the dataframes of the list
  final_points_intersection <- unique(purrr::reduce(final_points_intersection_list, dplyr::inner_join, by = names(final_points_intersection_list[[1]])))

  #resample the points dividing their number by the number of hypervolumes compared
  num_points_to_sample_in_intersection = nrow(final_points_intersection)/length(hv_points_ss_list)
  final_points_intersection = final_points_intersection[sample(1:nrow(final_points_intersection),
                                                               size = num_points_to_sample_in_intersection ), , drop = FALSE]

  final_volume_intersection <- nrow(final_points_intersection) / mindensity

  # get names
  # get column names assuming first set are correct for all
  cn <- dimnames(hv_list@HVList[[1]]@RandomPoints)[[2]]
  dn <- list(NULL, cn)

  #generate a new HV corresponding to the overlap between all the input HVs
  result = new("Hypervolume")
  result@Method = "Set operations"
  result@Data = matrix(NaN, nrow = 1, ncol = dim)
  result@Dimensionality = dim
  result@Volume = final_volume_intersection
  result@PointDensity = mindensity
  result@Parameters = list()
  result@RandomPoints = matrix(as.matrix(as.data.frame(final_points_intersection)),
                               ncol = dim)
  result@ValueAtRandomPoints = rep(1, nrow(result@RandomPoints))
  # set dimnames
  dimnames(result@RandomPoints) = dn
  dimnames(result@Data) = dn

  return(result)
}
