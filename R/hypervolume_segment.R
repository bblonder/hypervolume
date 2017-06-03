hypervolume_segment <- function(hv,distance.factor=1, num.points.max=NULL, verbose=TRUE, check.memory=TRUE)
{
  if(class(hv) != "Hypervolume")
  {
    stop("Input must be a Hypervolume class object")
  }
  
  # thin the input
  if (!is.null(num.points.max))
  {
    hv <- hypervolume_thin(hv, num.points=num.points.max)
  }
  
  hvrp <- hv@RandomPoints
  
  if (check.memory==TRUE)
  {
    stop(sprintf('Analysis will require storage of %d numbers. Re-run with check.memory=FALSE to continue.', nrow(hvrp)^2))
  }
  
  characteristicdistance <- (1/hv@PointDensity)^(1/hv@Dimensionality)
  
  if (verbose==TRUE)
  {
    cat('Performing hierarchical clustering... ')
  }
  hc <- fastcluster::hclust.vector(hvrp, method='single')
  if (verbose==TRUE)
  {
    cat('done.\n')
  }  
  
  if (verbose==TRUE)
  {
    cat('Segmenting clusters...')
  }
  membership <- cutree(hc, h = characteristicdistance * distance.factor)
  if (verbose==TRUE)
  {
    cat('Done.\n ')
  }  
  
  ngroups <- max(membership)
  
  hvs <- vector("list",ngroups)
  
  

  if (verbose==TRUE)
  {
    cat('Generating hypervolumes...')
  }  
  for (i in 1:ngroups)
  {
    hv_temp <- hv
    hv_temp@Data <- matrix(NaN,nrow=1,ncol=hv@Dimensionality)
    hv_temp@RandomPoints <- hvrp[membership==i,,drop=FALSE]
    hv_temp@ValueAtRandomPoints <- hv@ValueAtRandomPoints[membership==i]
    hv_temp@Volume <- hv@Volume * nrow(hv_temp@RandomPoints) / nrow(hv@RandomPoints)
    hv_temp@Method <- "Segmentation hypervolume"
    hv_temp@Name <- sprintf("%s (cluster %d/%d)", hv@Name, i, ngroups)
    hv_temp@Parameters = hv@Parameters
    
    hvs[[i]] <- hv_temp
  }
  
  hvs_segmented <- do.call("hypervolume_join",hvs)
  if (verbose==TRUE)
  {
    cat('done.\n')
  }
  
  return(hvs_segmented)
}