hypervolume_estimate_probability <- function(hv, points, reduction.factor=1, weight.exponent=-1, set.edges.zero=TRUE, edges.zero.distance.factor=1, verbose=TRUE)
{  
  np = nrow(hv@RandomUniformPointsThresholded)
  dimhv = ncol(hv@RandomUniformPointsThresholded)
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
  
  hv_points_ss = hv@RandomUniformPointsThresholded[index_ss,,drop=FALSE]
  prob_ss = hv@ProbabilityDensityAtRandomUniformPoints[index_ss]
  
  # determine the reduced hypervolume's point density
  point_density = nrow(hv_points_ss) / hv@Volume
  
  # calculate characteristic distances
  cutoff_dist = point_density^(-1/dimhv)
  
  
  if (verbose==TRUE)
  {
    cat(sprintf('Retaining %d/%d hypervolume random points for comparison with %d test points.\n', nrow(hv_points_ss), nrow(hv@RandomUniformPointsThresholded), nrow(points)))
  }
  
  pb <- progress_bar$new(total=nrow(points))
  if (verbose==TRUE)
  {
    pb$tick()
  }
  
  points <- as.matrix(points)
  dimnames(points) <- list(NULL, NULL)
  
  probabilities <- rep(NA, nrow(points))
  
  
  
  probabilities <- sapply(1:nrow(points), function(i) {
    if (verbose==TRUE & (i%%10==0))
    {
      pb$update(i/nrow(points))
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


hypervolume_estimate_probability_kdtree <- function(hv, points, reduction.factor=1, weight.exponent=-1, set.edges.zero=TRUE, edges.zero.distance.factor=1, radius=5*mean(apply(hv@RandomUniformPointsThresholded,2,sd)),verbose=TRUE)
{  
  np = nrow(hv@RandomUniformPointsThresholded)
  dimhv = ncol(hv@RandomUniformPointsThresholded)
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
  
  hv_points_ss = hv@RandomUniformPointsThresholded[index_ss,,drop=FALSE]
  prob_ss = hv@ProbabilityDensityAtRandomUniformPoints[index_ss]
  
  # determine the reduced hypervolume's point density
  point_density = nrow(hv_points_ss) / hv@Volume
  
  # calculate characteristic distances
  cutoff_dist = point_density^(-1/dimhv)
  
  
  if (verbose==TRUE)
  {
    cat(sprintf('Retaining %d/%d points for %d pairwise comparisons.\n', nrow(hv_points_ss), nrow(hv@RandomUniformPointsThresholded), nrow(points)*nrow(hv_points_ss)))
  }
  
  pb <- progress_bar$new(total=nrow(points))
  if (verbose==TRUE)
  {
    pb$tick()
  }
  
  probabilities <- rep(NA, nrow(points))
  
  # find points in the hypervolume that are within a ciritical distance of the test points
  # choose a radius that is a rough heuristic of 'far enough' away
  points_in_hv_all_list = evalfspherical(data=hv_points_ss, radius=radius, points=points, getid.nearestneighbor=TRUE,verbose=verbose)
  
  for (i in 1:nrow(points))
  {
    if (verbose==TRUE)
    {
      pb$tick()
    }
    
    indexvals = points_in_hv_all_list[[2]][[i]]
    
    indexvals_ss <- indexvals[indexvals > 0]
    
    if (length(indexvals_ss) > 0)
    {
      # for all random points near this test point, find the distance
      distances = pdist(points[i,,drop=FALSE], hv_points_ss[indexvals_ss,])@dist
      
      # inverse distance weight
      weights = distances^weight.exponent
      
      # do weighted sum of probabilities based on distance weights
      prob_ss = hv@ProbabilityDensityAtRandomUniformPoints[indexvals_ss]
      
      probabilities[i] <- sum(prob_ss * weights)/sum(weights)
      
      # force far-away points to zero probability
      if (set.edges.zero==TRUE)
      {
        if (min(distances) > cutoff_dist * edges.zero.distance.factor)
        {
          probabilities[i] <- 0
        }
      }
    }
    else
    {
      probabilities[i] <- 0
    }
  }
  
  return(probabilities)
}



hypervolume_estimate_probability_old <- function(hv, points, reduction.factor=1, distance.factor=1.0, verbose=TRUE)
{  
  np = nrow(hv@RandomUniformPointsThresholded)
  dimhv = ncol(hv@RandomUniformPointsThresholded)
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
  
  hv_points_ss = hv@RandomUniformPointsThresholded[index_ss,,drop=FALSE]
  prob_ss = hv@ProbabilityDensityAtRandomUniformPoints[index_ss]
  
  
  if (verbose==TRUE)
  {
    cat(sprintf('Retaining %d/%d points for %d probability estimates.\n', nrow(hv_points_ss), nrow(hv@RandomUniformPointsThresholded), nrow(points)))
  }
  
  # determine the reduced hypervolume's point density
  point_density = nrow(hv_points_ss) / hv@Volume
  
  # calculate characteristic distances
  cutoff_dist = point_density^(-1/dimhv) * distance.factor

  
  
  if (verbose==TRUE)
  {
    cat(sprintf('Estimating probability... \n'))
  }
  
  points_in_hv_all_list = evalfspherical(data=hv_points_ss, radius=cutoff_dist, points=points, getid.nearestneighbor=TRUE,verbose=verbose)
  points_in_hv_all_list_indexes <- points_in_hv_all_list[[2]]
  
  probabilities = sapply(points_in_hv_all_list_indexes, function(x) { 
    x_this <- x[x>0] # keep only non-NA points
    median_probability = median(prob_ss[x_this])
    return(median_probability)
    })
  
  probabilities[is.na(probabilities)] <- 0
  
  if (verbose==TRUE)
  {
    cat(sprintf('done.'))
  }
  
  return(probabilities)
}