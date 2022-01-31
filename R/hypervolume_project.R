hypervolume_project <- function(hv, rasters, type = "probability", 
                                verbose = TRUE, ...) {
  raster.values = data.frame(raster::getValues(rasters))
  if (type == "probability") {
    projected.values <- hypervolume_estimate_probability(hv = hv,
                                                         points = raster.values, 
                                                         verbose = verbose,
                                                         ...)
  } else if (type == "inclusion") {
    projected.values = hypervolume_inclusion_test(hv, raster.values, 
                                                  verbose = verbose, ...)
  } else {
    stop("Unsupported 'type' argument")
  }
  # predicted values to raster. Ignore NA cells (if any)
  raster.out <- rasters[[1]]
  if (length(attr(projected.values, "NAIdx")) > 0) {
    naIdx <- attr(projected.values, "NAIdx")
    raster.out[-naIdx] <- projected.values
  } else {
    raster.out <- projected.values
  }
  return(raster.out)
}
