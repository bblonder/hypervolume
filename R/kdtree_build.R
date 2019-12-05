kdtree_build <- function(data, verbose=TRUE)
{
<<<<<<< HEAD
  if (any(class(data) == "data.frame"))
=======
  if (class(data) %in% "data.frame")
>>>>>>> 38023b013d27b65e6c68e23a4d71de701338a22e
  {
    data <- as.matrix(data)
  }
  
<<<<<<< HEAD
  if (any(class(data) == "matrix"))
=======
  if (class(data) %in% "matrix")
>>>>>>> 38023b013d27b65e6c68e23a4d71de701338a22e
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