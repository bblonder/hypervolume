hypervolume_estimate_probability <- function(hv, points, reduction_factor=1, verbose=T, distance_factor=1.0)
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
  
  # figure out which points are 'close enough' to other points
  # (within a n-ball of the critical distance)
  points_in_hv_all_list = evalfspherical(hv_points_ss, cutoff_dist, points, verbose=verbose,getid.nearestneighbor=TRUE)
  points_in_hv_all_list_probs <- points_in_hv_all_list[[2]]
  
  if (verbose==TRUE)
  {
    cat(sprintf('Estimating probability at %d points...', length(points_in_hv_all_list_probs)))
  }
  
  df_probs <- as.data.frame(do.call("rbind",lapply(1:length(points_in_hv_all_list_probs), function(i) { cbind(index.point=i, index.prob=points_in_hv_all_list_probs[[i]])})))
  df_probs <- df_probs[df_probs$index.prob > -1,] # remove all 'out' cases
  df_probs$prob <- hv@ProbabilityDensityAtRandomUniformPoints[df_probs$index.prob+1]
  
  probabilities <- tapply(df_probs$prob, df_probs$index.point, mean)
  
  probabilities_out <- rep(0, length(points_in_hv_all_list_probs))
  print(str(probabilities_out))
  probabilities_out[as.numeric(names(probabilities))] <- probabilities
  print(str(probabilities_out))
  return(probabilities_out)
  
  
  # look up probabilities for all of the random points that are 'in' the ball
  probabilities <- sapply(points_in_hv_all_list_probs, function(indices) { 
    cat('.')
    indices_startatone <- indices+1 # C to R index conversion
    mean_probability <- mean(hv@ProbabilityDensityAtRandomUniformPoints[indices_startatone]) # any remaining zero indices are ignored
    return(mean_probability)
    })
  if (verbose==TRUE)
  {
    cat(sprintf('done.\n'))
  }
  
  # flag those points and return them
  probabilities[points_in_hv_all_list[[1]] == 0] <- 0
  
  return(probabilities)
}