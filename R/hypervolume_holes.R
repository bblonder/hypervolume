hypervolume_holes <- function(hv.obs, hv.exp, set.num.points.max=NULL, set.check.memory=TRUE)
{		
  # initialize result
  finalresult <- NULL
  
  if (is.null(set.num.points.max))
  {
  set.num.points.max = ceiling(10^(3+sqrt(hv.obs@Dimensionality)))
    cat(sprintf('Choosing set.num.points.max=%.0f (choose a larger value for more accuracy.)\n',set.num.points.max))    
  }
  
  if (hv.obs@Dimensionality != hv.exp@Dimensionality)
  {
    stop('Observed and expected hypervolumes must have same dimensionality.')
  }
  
  # FIND THE DIFFERENCE between the convex hull shape and the real hypervolume
  cat("Beginning set operations (resampling to minimum density)...")
  hvs_overlap <- hypervolume_set(hv.obs, hv.exp, check.memory=set.check.memory, num.points.max=set.num.points.max)
  if (set.check.memory==TRUE)
  {
    stop('Set set.check.memory=F to continue.\n')
  }
  cat("Finished set operations.\n")
  
  # find the distance between all points in the difference
  randompoints <- hvs_overlap@HVList$Unique_2@RandomPoints
  if (is.null(randompoints) || nrow(randompoints) == 0)
  {		
    cat('No holes found.\n');
  }
  else
  {  	
    cat(sprintf("Retaining %.0f random points in set difference.\n", nrow(randompoints)))

    criticaldistance <- hvs_overlap@HVList$Unique_2@PointDensity ^(-1/hvs_overlap@HVList$Unique_2@Dimensionality)    
    
    # find points with minimum neighbor distance less than threshold
    distances <- as.matrix(dist(randompoints, method="euclidean"))
    diag(distances) <- NA
    isin <- (apply(distances, 1, min, na.rm=TRUE) < criticaldistance)
    
    cat(sprintf("Removing %d stray points...\n", length(which(isin==0))))
    
    randompoints_trimmed <- randompoints[isin,]
    
    thishv <- hvs_overlap@HVList$Unique_2 # copy base information
    thishv@RandomPoints <- randompoints_trimmed
    thishv@Volume <- hvs_overlap@HVList$Unique_2@Volume * nrow(randompoints_trimmed) / nrow(randompoints)
    thishv@Name <- sprintf("Hole in %s relative to %s", hv.obs@Name, hv.exp@Name)
    
    finalresult <- thishv

    # return the final hypervolumelist
    cat('Returning all holes.\n');
    
    if (is.null(finalresult))
    {
      cat('No holes found. Function will return NULL.\n')
    }
        
  }
  
  return(finalresult)
}