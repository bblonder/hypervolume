kdtree_build <- function(data, verbose=TRUE)
{
  if (any(inherits(data,"data.frame")))
  {
    data <- as.matrix(data)
  }
  if (any(inherits(data,"matrix")))
  {
    if (verbose==TRUE)
    {
      cat("\nBuilding tree... \n")
    }
    
    # No transpose: convert.cpp reads R's col-major buffer directly and
    # writes a row-major flat vector<double> into the KDTree in one pass.
    # The old code did t(data) here so that convertMatrixToVector could
    # index it as row-major, which meant an R-side transpose allocation
    # plus a per-element re-copy inside C++.
    kdt <- kdtree_build_intl(data, nrow(data), ncol(data))
    
    if (verbose==TRUE)
    {
      cat("done.\n")
    }
    
    return(kdt)
  }
  else
  {
    stop("Input data not a matrix or data frame")
  }
}