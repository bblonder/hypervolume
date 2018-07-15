hypervolume_estimate_probability <- function(hv, points, reduction.factor=1, weight.exponent=-1, set.edges.zero=TRUE, edges.zero.distance.factor=1, verbose=TRUE)
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
  index_ss = sample(1:np,size=numpointstokeep_hv,replace=FALSE)
  
  hv_points_ss = hv@RandomPoints[index_ss,,drop=FALSE]
  prob_ss = hv@ValueAtRandomPoints[index_ss]
  
  # determine the reduced hypervolume's point density
  point_density = nrow(hv_points_ss) / hv@Volume
  
  # calculate characteristic distances
  cutoff_dist = point_density^(-1/dimhv)
  
  
  if (verbose==TRUE)
  {
    cat(sprintf('Retaining %d/%d hypervolume random points for comparison with %d test points.\n', nrow(hv_points_ss), nrow(hv@RandomPoints), nrow(points)))
  }
  
  pb <- progress_bar$new(total=nrow(points))
  if (verbose==TRUE)
  {
    pb$tick(0)
  }
  
  points <- as.matrix(points)
  dimnames(points) <- list(NULL, NULL)
  
  probabilities <- rep(NA, nrow(points))
  
  
  
  probabilities <- sapply(1:nrow(points), function(i) {
    if (verbose==TRUE & (i%%10==0))
    {
      if (!pb$finished==TRUE)
      {
        pb$update(i/nrow(points))
      }
    }
    
    distances <- pdist(points[i,,drop=FALSE], hv_points_ss)@dist
    weights <- distances^weight.exponent
    
    result <- sum(prob_ss * weights)/sum(weights)
    
    if (set.edges.zero==TRUE)
    {
      if (min(distances) > cutoff_dist * edges.zero.distance.factor)
      {
        result <- 0
      }
    }
    
    return(result)
  })
  
  return(probabilities)
}
