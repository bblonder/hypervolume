hypervolume_estimate_probability <- function(hv, points, reduction_factor=1, verbose=TRUE, distance_factor=1.0, chunksize=1e4)
{  
  np = nrow(hv@RandomUniformPointsThresholded)
  dimhv = ncol(hv@RandomUniformPointsThresholded)
  dimp = ncol(points)
  
  if (dimp != dimhv)
  {
    stop('Dimensionality of hypervolume and points is not the same.')
  }
  
  if (reduction_factor <= 0 | reduction_factor > 1)
  {
    stop('Reduction factor is not in (0,1].')
  }
  
  # now pick a uniformly random subset of these points
  # assuming that the set of points is already uniformly random
  numpointstokeep_hv = floor(np * reduction_factor)        
  if (reduction_factor < 1)
  {
    hv_points_ss = hv@RandomUniformPointsThresholded[sample(1:np,size=numpointstokeep_hv),]
  }
  else
  {
    hv_points_ss = hv@RandomUniformPointsThresholded
  }
  
  if (verbose==TRUE)
  {
    cat(sprintf('Retaining %d points for %d inclusion tests.\n', numpointstokeep_hv, nrow(points)))
  }
  
  # determine the reduced hypervolume's point density
  point_density = nrow(hv_points_ss) / hv@Volume
  
  # calculate characteristic distances
  cutoff_dist = point_density^(-1/dimhv) * distance_factor
  
  if (verbose==TRUE)
  {
    cat(sprintf('Estimating probability... '))
  }
  
  
  chunksize <- min(chunksize, nrow(points))  
  
  # figure out which points are 'close enough' to other points
  # (within a n-ball of the critical distance)
  num.samples.completed <- 0
  num.chunks <- ceiling(nrow(points)/chunksize)
  probabilities <- vector(mode="list",length=num.chunks) # the randomuniform points

  for (i in 1:num.chunks)
  {
    cat(sprintf("%.3f ", i/num.chunks))  
    
    points_in_hv_all_list = evalfspherical(hv_points_ss, cutoff_dist, points[(num.samples.completed+1):min(num.samples.completed+chunksize, nrow(points)),,drop=FALSE], verbose=verbose,getid.nearestneighbor=TRUE)
    points_in_hv_all_list_probs <- points_in_hv_all_list[[2]]
    
    df_probs <- as.data.frame(do.call("rbind",lapply(1:length(points_in_hv_all_list_probs), function(i) { cbind(index.point=i, index.prob=points_in_hv_all_list_probs[[i]])})))
    df_probs$index.prob[df_probs$index.prob<=0] <- NA

    df_probs$prob <- hv@ProbabilityDensityAtRandomUniformPoints[df_probs$index.prob]
    df_probs$prob[is.na(df_probs$prob)] <- 0
    
    mean_probs <- tapply(df_probs$prob, df_probs$index.point, mean)
    
    probabilities[[i]] <- mean_probs
    
    num.samples.completed = num.samples.completed + chunksize
    
  }
  probabilities <- do.call("c",probabilities)
  names(probabilities) <- NULL
  
  return(probabilities)
}