hypervolume_join <- function(...)
{
  args <- list(...)
  if (length(args)==1 & is.list(args[[1]]))
  {
    result <- new("HypervolumeList",HVList=args[[1]])
    return(result)
  }
  else
  {
    hvl <- list()
    
    for (a in list(...))
    {
      if (inherits(a,"HypervolumeList"))
      {
        hvl <- c(hvl, a@HVList)
      }
      else if (inherits(a,"Hypervolume"))
      {
        hvl <- c(hvl, a)
      }
    }
    
    return(new("HypervolumeList",HVList=hvl))
  }
}