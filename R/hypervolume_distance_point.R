hypervolume_distance_point <- function(hv1, x, type="minimum", num.points.max=1000, check.memory=TRUE){
  
  # Ensure that hv1 a Hypervolume object
  if (!inherits(hv1, "Hypervolume")) {
    
    stop("hv1 must be a Hypervolume object.")
  }
  
  hv1p <- hv1@RandomPoints
  
  # Coerce `x` to a matrix
  if (is.vector(x)) {
    
    pt <- matrix(x, nrow = 1) # Coerce point to a one ligne matrix
    
  } else {
    
    tryCatch({
      
      pt <- as.matrix(x)
      
    }, error = function(e){stop("pt must be an object coercible to a matrix object.")})
    
  }
  
  if (type == "minimum") {
    
    hv1p_ss <- hv1p[sample(1:nrow(hv1p), min(num.points.max, nrow(hv1p))), , drop=FALSE]
    
    if (check.memory == TRUE) {
      
      cat(sprintf('Calculation will require %d pairwise distance calculations.\n', nrow(hv1p_ss) * nrow(pt)))
      
      message('Re-run with check.memory=FALSE to continue.')
      stop()
      
    }
    
    crossdistances <- fastPdist2(hv1p_ss, pt)
    
    minimum_distance <- min(as.numeric(as.matrix(crossdistances)), na.rm=TRUE)
    
    return(minimum_distance)
    
  } else if (type == "maximum"){
    
    hv1p_ss <- hv1p[sample(1:nrow(hv1p), min(num.points.max, nrow(hv1p))), , drop=FALSE]
    
    if (check.memory == TRUE) {
      
      cat(sprintf('Calculation will require %d pairwise distance calculations.\n', nrow(hv1p_ss) * nrow(pt)))
      
      message('Re-run with check.memory=FALSE to continue.')
      stop()
      
    }
    
    crossdistances <- fastPdist2(hv1p_ss, pt)
    
    maximum_distance <- max(as.numeric(as.matrix(crossdistances)), na.rm=TRUE)
    
    return(maximum_distance)
    
  } else {
    
    stop('Argument \'type\' takes unrecognized value.')
    
  }
}

