hypervolume_prune <- function(hvlist, num.points.min=NULL, volume.min=NULL, return.ids=FALSE)
{
  if (!is.null(num.points.min))
  {
    if (num.points.min < 0)
    {
      stop("num.points.min must be at least zero.")
    }
    if (!is.null(volume.min))
    {
      stop("Cannot specify both volume.min and num.points.min.")
    }
  }
  else if (!is.null(volume.min))
  {
    if (volume.min < 0)
    {
      stop("volume.min must be at least zero.")
    }
    
    if (!is.null(num.points.min))
    {
      stop("Cannot specify both volume.min and num.points.min.")
    }
  }
  else
  {
    stop("Must specify either volume.min or num.points.min")
  }  
  
  if(class(hvlist) != "HypervolumeList")
  {
    stop("Input hvlist must be of class HypervolumeList.")
  }
  
  # do segmentation
  dodrop <- rep(FALSE, length(hvlist@HVList))
  for (i in 1:length(hvlist@HVList))
  {
    np <- nrow(hvlist@HVList[[i]]@RandomPoints)
    vol <- hvlist@HVList[[i]]@Volume
    
    if (!is.null(num.points.min))
    {
      if (np < num.points.min)
      {
        dodrop[i] <- TRUE
      }
    }
    if (!is.null(volume.min))
    {
      if (vol < volume.min)
      {
        dodrop[i] <- TRUE
      }
    }
  }
  
  hvlist@HVList <- hvlist@HVList[!dodrop]
  
  if (return.ids)
  {
    return(list(HVList=hvlist, IDs=which(!dodrop)))
  }
  else
  {
    return(hvlist)
  }
}