hypervolume_thin <- function(hv, factor=NULL, num.points=NULL)
{
  rp <- hv@RandomPoints
  nrp <- nrow(rp)  
  
  if (!is.null(factor))
  {
    if (factor <= 0 | factor >=1)
    {
      stop("Thinning factor must be in (0,1).")
    }
  }
  else if (!is.null(num.points))
  {
    if (num.points < 1)
    {
      stop("Number of points must be greater than zero.")
    }

    if (!is.null(factor))
    {
      stop("Cannot specify both factor and num.points.")
    }
    else
    {
      # make sure we don't take more points than exist in the dataset
      num.points <- min(nrp, num.points)
      
      # recast in terms of a factor
      factor <- num.points / nrp
    }
  }
  else
  {
    stop("Must specify either factor or num.points.")
  }
  
  hv_out <- hv
  
  hv_out@RandomPoints <- rp[sample(1:nrow(rp),nrow(rp)*factor),]
  
  hv_out@PointDensity <- hv@PointDensity * factor
  
  return(hv_out)
}