weight_data <- function(data, weights, jitter.sd=matrix(0, nrow=nrow(data),ncol=ncol(data)))
{
  data <- as.data.frame(data)
  jitter.sd <- as.matrix(jitter.sd)
  
  colClasses <- lapply(data, class)
  
  if (any(colClasses!="numeric"))
  {
    stop("All columns must be numeric type.")
  }
  
  if (length(weights)==1)
  {
    weights <- rep(weights, nrow(data))
  }
  
  if (nrow(data) != length(weights))
  {
    stop('Weights must be of same length as data.')
  }

  
  if (!all(weights==as.integer(weights)))
  {
    warning("Conversion of weights to integers resulted in rounding. Are all weights integers of at least one?")
  }
  weights <- as.integer(weights)
  
  if (any(weights<=0))
  {
    stop('All weights must be positive.')
  }
  
  if (length(jitter.sd)==1)
  {
    jitter.sd=matrix(jitter.sd, nrow=nrow(data),ncol=ncol(data))
  }
  else if (length(jitter.sd)==ncol(data))
  {
    jitter.sd=matrix(jitter.sd,nrow=nrow(data),ncol=ncol(data),byrow=TRUE)
  }
  else if (all(dim(data)==dim(jitter.sd)))
  {
    # do nothing; everything OK
  }
  else
  {
    stop('Parameter jitter.sd cannot be conformed to match data.')
  }
  
  data.expanded <- data[rep(seq(nrow(data)), weights), ]
  
  jitter.sd.expanded <- jitter.sd[rep(seq(nrow(jitter.sd)), weights), ]
  
  jitter.sd.vals <- matrix(rnorm(length(jitter.sd.expanded),mean=0,sd=jitter.sd.expanded),
                           nrow=nrow(jitter.sd.expanded),
                           ncol=ncol(jitter.sd.expanded))
  
  data.expanded.jittered <- data.expanded + jitter.sd.vals
  
  return(data.expanded.jittered)
  
}