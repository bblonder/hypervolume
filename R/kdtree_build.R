kdtree_build <- function(data, verbose=TRUE)
{
  if (class(data) == "data.frame")
  {
    data <- as.matrix(data)
  }
  
  if (class(data) == "matrix")
  {
    if (verbose==TRUE)
    {
      cat("\nBuilding tree... \n")
    }
    
    kdt <- kdtree_build_intl(t(data),nrow(data),ncol(data),verbose)
    
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