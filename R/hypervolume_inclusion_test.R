hypervolume_inclusion_test <- function(hv, points, reduction.factor=1, fast.or.accurate='fast', fast.method.distance.factor=1.0, accurate.method.threshold=quantile(hv@ValueAtRandomPoints,0.5), verbose=TRUE, ...)
{  
  np = nrow(hv@RandomPoints)
  dimhv = ncol(hv@RandomPoints)
  dimp = ncol(points)
  
  if (dimp != dimhv)
  {
    stop('Dimensionality of hypervolume and points is not the same.')
  }
  
  if (reduction.factor <= 0 | reduction.factor > 1)
  {
    stop('Reduction factor is not in (0,1].')
  }

  # now pick a uniformly random subset of these points
  # assuming that the set of points is already uniformly random
  numpointstokeep_hv = floor(np * reduction.factor)        
  if (reduction.factor < 1)
  {
    hv_points_ss = hv@RandomPoints[sample(1:np,size=numpointstokeep_hv),]
  }
  else
  {
    hv_points_ss = hv@RandomPoints
  }
  
  if (fast.or.accurate=="fast")
  {
    warning('Results may have a high error rate. Consider setting fast.or.accurate=\'accurate\'.')
    if (verbose==TRUE)
    {
      cat(sprintf('Retaining %d/%d hypervolume random points for comparison with %d test points.\n', nrow(hv_points_ss), nrow(hv@RandomPoints), nrow(points)))
    }
    
    # determine the reduced hypervolume's point density
    point_density = nrow(hv_points_ss) / hv@Volume
    
    # calculate characteristic distances
    cutoff_dist = point_density^(-1/dimhv) * fast.method.distance.factor
    
    if (nrow(hv_points_ss) > 0)
    {
      # figure out which points are 'close enough' to other points
      # (within a n-ball of the critical distance)
      points_in_hv_all = evalfspherical(hv_points_ss, cutoff_dist, points, verbose=verbose)
      
      # keep any point that has at least one random point within it (somewhat liberal)
      points_in = points_in_hv_all > 0
      
      return(points_in)
    }
    else
    {
      warning('No points in hypervolume - increase reduction.factor or use non-empty input hypervolume! Returning NULL.')
      return(NULL)
    }
  }
  else if (fast.or.accurate=='accurate')
  {
    probabilities <- hypervolume_estimate_probability(hv, points, reduction.factor=reduction.factor, verbose=verbose, ...)
    
    inclusion_result <- (probabilities >= accurate.method.threshold)
    
    return(inclusion_result)
  }
  else
  {
    stop('Unsupported method for fast.or.accurate.')
  }
}