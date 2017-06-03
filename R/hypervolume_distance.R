hypervolume_distance <- function(hv1, hv2, type="centroid", num.points.max=1000, check.memory=TRUE)
{
  hv1p <- hv1@RandomPoints
  hv2p <- hv2@RandomPoints
  
  if (type=="centroid")
  {
    hv1p_center <- colMeans(hv1p, na.rm=TRUE)
    hv2p_center <- colMeans(hv2p, na.rm=TRUE)
    
    centroid_distance <- sqrt(sum((hv1p_center - hv2p_center)^2))
    
    return(centroid_distance)
  }
  else if (type=="minimum")
  {    
    hv1p_ss <- hv1p[ sample(1:nrow(hv1p), min(num.points.max, nrow(hv1p)))  ,]
    hv2p_ss <- hv2p[ sample(1:nrow(hv2p), min(num.points.max, nrow(hv2p)))  ,]
    

    if (check.memory==TRUE)
    {
      cat(sprintf('Calculation will require %d pairwise distance calculations.\n',nrow(hv1p_ss)*nrow(hv2p_ss)))
      
      message('Re-run with check.memory=FALSE to continue.')
      stop()
    }

    crossdistances <- fastPdist2(hv1p_ss, hv2p_ss)
    
    minimum_distance <- min(as.numeric(as.matrix(crossdistances)),na.rm=TRUE)
    
    return(minimum_distance)
  }
  else
  {
    stop('Argument \'type\' takes unrecognized value.')
  }
}