hypervolume_join <- function(...)
{
  hvl <- list()
  
  for (a in list(...))
  {
    if (class(a) == "HypervolumeList")
    {
      hvl <- c(hvl, a@HVList)
    }
    else if (class(a) == "Hypervolume")
    {
      hvl <- c(hvl, a)
    }
  }
  
  return(new("HypervolumeList",HVList=hvl))
}