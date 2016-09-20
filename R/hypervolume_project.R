hypervolume_project <- function(hv, rasters, do.probability=TRUE, normalize=TRUE)
{
  raster.values = data.frame(raster::getValues(rasters))
  
  if (do.probability==TRUE)
  {
    projected.values <- hypervolume_estimate_probability(hv, raster.values)
    if (normalize==TRUE)
    {
      projected.values <- projected.values / sum(projected.values)
    }    
  }
  else
  {
    projected.values = hypervolume_inclusion_test(hv, raster.values)
  }

  
  # copy raster properties
  raster.out = raster::raster(rasters[[1]]);
  raster::values(raster.out) <- projected.values
  
  return(raster.out)
}