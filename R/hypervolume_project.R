hypervolume_project <- function(hv, rasters, type='probability', verbose=TRUE, ...)
{
  raster.values = data.frame(raster::getValues(rasters))
  
  if (type=='probability')
  {
    projected.values <- hypervolume_estimate_probability(hv, raster.values,verbose=verbose, ...)
  }
  else if (type=='inclusion')
  {
    projected.values = hypervolume_inclusion_test(hv, raster.values, verbose=verbose, ...)
  }
  else
  {
    stop('Unsupported \'type\' argument')
  }

  
  # copy raster properties
  raster.out = rasters[[1]];
  raster::values(raster.out) <- projected.values
  
  return(raster.out)
}