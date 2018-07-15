sample_model_rejection <- function(model, range, N.samples, chunk.size=1e3, verbose=TRUE, min.value=0, ...) # range should be the output of e.g. apply(data, 2, range) (2 x n matrix)
{
  # determine dimensionality
  d = ncol(range)
  
  if (is.null(dimnames(range)[[2]]))
  {
    dimnames(range) <- list(NULL,paste("X",1:d,sep=""))
  }
  
  if (verbose==TRUE)
  {
    pb <- progress_bar$new(total = N.samples)
  }
  
  samples = list()
  total_accepted <- 0
  while(total_accepted < N.samples)
  {
    if (verbose==TRUE)
    {
      if (!pb$finished==TRUE)
      {
        pb$update(total_accepted/N.samples)
      }
    }
    
    # generate random points
    random_points <- data.frame(matrix(data=runif(chunk.size*d),ncol=d,nrow=chunk.size,dimnames=list(NULL, dimnames(range)[[2]])))
    
    # rescale points
    for (i in 1:d)
    {
      random_points[,i] <- random_points[,i]*(range[2,i] - range[1,i]) + range[1,i]
    }
    
    # predict values at these points
    samples_this = as.numeric(predict(model, newdata=random_points, ...))
    
    # store the new values
    samples = c(samples, list(cbind(random_points, samples_this)))
    
    # store the number of accepted points
    total_accepted = total_accepted + length(which(samples_this>min.value))
  }
  if (verbose==TRUE)
  {
    pb$terminate()
  }
  
  samples <- as.matrix(rbindlist(samples))

  dimnames(samples) <- list(NULL, c(dimnames(range)[[2]],"value"))
  
  # cutoff the samples to return the correct number of desired points
  index_cutoff = which(cumsum(samples[,ncol(samples)]>min.value) >= N.samples)[[1]]
  samples <- samples[1:index_cutoff,]
  
  return(samples)
}












sample_ellipsoid = function(center, n, scales) {
  k = length(center)
  points = matrix(rnorm(n * k), nrow = n)   # Start with noise
  points = points / sqrt(rowSums(points^2)) # Project to sphere surface
  radii = runif(n, 0)^(1/k) # Each point has a radius, prior to scaling
  for (i in 1:k) {
    points[,i] = points[,i] * radii * scales[i] + center[i]
  }
  return(points)
}

ellipsoid_volume = function(scales) {
  nball_volume(n = length(scales), r = 1) * prod(scales)
}

ellipsoid_inverse_weight = function(samples, centers, scales, verbose) {
  # Returns the number of center-based ellipsoids that contain each sampled point
  
  # scale to unit ball
  for (i in 1:ncol(centers)) {
    samples[ , i] = samples[ , i] / scales[i]
    centers[ , i] = centers[ , i] / scales[i]
  }
  
  # Count the overlap with unit ball
  tree <- kdtree_build(centers,verbose=verbose)
  query <- kdtree_ball_query_multiple(tree, t(samples), 
                                      nrow(samples), ncol(samples), 
                                      r = 1, verb = verbose)
  rm(tree)
  
  return(query)
}






















sample_model_ellipsoid <- function(predict_function=NULL, data, scales, min.value, samples.per.point, chunk.size=1e3, verbose=TRUE, return.full=FALSE)
{
  # Use only complete cases
  data = na.omit(as.matrix(data))
  
  # determine dimensionality
  d = ncol(data)
  
  N.samples <- ceiling(samples.per.point * nrow(data))
  
  if (is.null(dimnames(data)[[2]]))
  {
    dimnames(data) <- list(NULL,paste("X",1:d,sep=""))
  }
  
  if (verbose==TRUE)
  {
    pb <- progress_bar$new(total = N.samples)
    pb$tick(0)
  }
  
  samples = list()
  volume_sampling_extent_all = list()
  total_accepted <- 0
  total_tried <- 0
  
  while(total_accepted < N.samples)
  {
    if (verbose==TRUE)
    {
      if (!pb$finished==TRUE)
      {
        pb$update(total_accepted/N.samples)
      }
    }
    
    ## STEP ONE: Collect samples from ellipsoids around each data point
    full_samples = lapply(1:nrow(data),
                          function(i)
                          {
                            se = sample_ellipsoid(data[i, ], 
                                                  chunk.size, 
                                                  scales = scales)
                            return(data.frame(se))
                          })
    full_samples = as.matrix(rbindlist(full_samples))
    
    
    
    # Discard samples from regions that were over-sampled.
    iweight = ellipsoid_inverse_weight(full_samples, centers = data, 
                                       scales = scales, verbose = verbose)
    # scale weights so that the largest value is 1
    weight = 1/(iweight / min(iweight))
    # This is the average weight of the samples
    mean_weight = sum(1/iweight) / nrow(full_samples)
    # Total volume is the volume of one ellipse, times the number of ellipse,
    # times the proportion of samples that can be retained, times the fraction of points included above threshold
    volume_sampling_extent = ellipsoid_volume(scales) * nrow(data) * mean_weight
    
    # resample sampled points down to uniform density
    included = as.logical(rbinom(length(iweight), size = 1, prob = weight))
    
    # now we have a uniform grid of points around each data point, samples_retained
    samples_retained = full_samples[included, , drop = FALSE]
    
    ### STEP TWO: estimate function
    # predict function value at each point
    predicted_values <- predict_function(samples_retained)

    included_thresholded = ( as.numeric(predicted_values) > min.value )
    
    samples_retained_thresholded = samples_retained[included_thresholded, , drop=FALSE]
    predicted_values_thresholded = predicted_values[included_thresholded]
    
    samples_final_this = cbind(samples_retained_thresholded, predicted_values_thresholded)
    dimnames(samples_final_this) <- list(NULL, c(dimnames(data)[[2]],"value"))
    
    # count progress
    total_tried <- total_tried + nrow(samples_retained) # for eventual volume counting
    total_accepted <- total_accepted + nrow(samples_retained_thresholded) 
    
    # store loop output
    samples <- c(samples, list(data.frame(samples_final_this))) # must be df format for rbindlist
    volume_sampling_extent_all <- c(volume_sampling_extent_all, volume_sampling_extent)
  }
  
  # concatenate results  
  samples <- as.matrix(rbindlist(samples))
  samples <- samples[sample(1:nrow(samples),N.samples),,drop=FALSE]
  # calculate volumes
  volume_sampling_extent_all_mean = mean(unlist(volume_sampling_extent_all),na.rm=T)
  volume = volume_sampling_extent_all_mean * total_accepted / total_tried
  
  if (verbose==TRUE)
  {
    pb$terminate()
  }
  
  if (return.full==TRUE)
  {
    return(list(samples = samples, full_samples=full_samples, volume=volume))
  }
  else
  {
    return(list(samples = samples, volume=volume))
  }
}




sample_model_ellipsoid_dave_broken <- function(predict_function=NULL, data, scales, min.value, samples.per.point, chunk.size=1e3, verbose=TRUE, return.full=FALSE)
{
  # Use only complete cases
  data = na.omit(as.matrix(data))
  
  # determine dimensionality
  d = ncol(data)
  
  N.samples <- ceiling(samples.per.point * nrow(data))
  
  if (is.null(dimnames(data)[[2]]))
  {
    dimnames(data) <- list(NULL,paste("X",1:d,sep=""))
  }
  
  if (verbose==TRUE)
  {
    pb <- progress_bar$new(total = N.samples)
    pb$tick(0)
  }
  
  samples = list()
  volume_sampling_extent_all = list()
  total_accepted <- 0
  total_tried <- 0
  
  while(total_accepted < N.samples)
  {
    if (verbose==TRUE)
    {
      if (!pb$finished==TRUE)
      {
        pb$update(total_accepted/N.samples)
      }
    }
    
    ## STEP ONE: Collect samples from ellipsoids around each data point
    full_samples = lapply(1:nrow(data),
                          function(i)
                          {
                            se = sample_ellipsoid(data[i, ], 
                                                  chunk.size, 
                                                  scales = scales)
                            return(data.frame(se))
                          })
    full_samples = as.matrix(rbindlist(full_samples))
    
    ### STEP TWO: estimate function
    # predict function value at each point
    predicted_values <- predict_function(full_samples)
    
    included_thresholded = (as.numeric(predicted_values) > min.value)
    
    samples_retained_thresholded = full_samples[included_thresholded, , drop=FALSE]
    predicted_values_thresholded = predicted_values[included_thresholded]
    
    # Discard samples from regions that were over-sampled.
    iweight = ellipsoid_inverse_weight(samples_retained_thresholded, centers = data, 
                                       scales = scales, verbose = verbose)
    # scale weights so that the largest value is 1
    weight = 1/(iweight / min(iweight))
    # This is the average weight of the samples
    mean_weight = sum(1/iweight) / nrow(full_samples)
    # Total volume is the volume of one ellipse, times the number of ellipse,
    # times the proportion of samples that can be retained, times the fraction of points included above threshold
    volume_sampling_extent = ellipsoid_volume(scales) * nrow(data) * mean_weight
    
    # resample sampled points down to uniform density
    included = as.logical(rbinom(length(iweight), size = 1, prob = weight))
    
    # now we have a uniform grid of points around each data point
    samples_retained_final = samples_retained_thresholded[included, , drop = FALSE]
    predicted_values_final = predicted_values_thresholded[included]

    samples_final_this = cbind(samples_retained_final, predicted_values_final)
    dimnames(samples_final_this) <- list(NULL, c(dimnames(data)[[2]],"value"))
    
    # count progress
    total_tried <- total_tried + nrow(samples_retained_thresholded) # for eventual volume counting
    total_accepted <- total_accepted + nrow(samples_final_this) 
    
    # store loop output
    samples <- c(samples, list(data.frame(samples_final_this))) # must be df format for rbindlist
    volume_sampling_extent_all <- c(volume_sampling_extent_all, volume_sampling_extent)
  }

  # concatenate results  
  samples <- as.matrix(rbindlist(samples))
  samples <- samples[sample(1:nrow(samples),N.samples),,drop=FALSE]
  # calculate volumes
  volume_sampling_extent_all_mean = mean(unlist(volume_sampling_extent_all),na.rm=T)
  volume = volume_sampling_extent_all_mean * total_accepted / total_tried

  if (verbose==TRUE)
  {
    pb$terminate()
  }
  
  if (return.full==TRUE)
  {
    return(list(samples = samples, full_samples=full_samples, volume=volume))
  }
  else
  {
    return(list(samples = samples, volume=volume))
  }
}



