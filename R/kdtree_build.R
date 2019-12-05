kdtree_build <- function(data, verbose=TRUE)
{
  if (class(data) %in% "data.frame")
  {
    data <- as.matrix(data)
  }
  
  if (class(data) %in% "matrix")
  {
    if (verbose==TRUE)
    {
      cat("\nBuilding tree... \n")
    }
    
    kdt <- kdtree_build_intl(t(data),nrow(data),ncol(data))
    
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