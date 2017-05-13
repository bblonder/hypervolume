
rwmetro <- function(target,N,x,VCOV, verbose=TRUE) # target = function to sample, N = number of samples, x starting point, VCOV the jump neighborhood for a gaussian
{
  if (verbose==TRUE)
  {
    pb <- progress_bar$new(total = 100)
    pb$tick()
  }
  
  # initialize storage
  samples = matrix(NA, nrow=N,ncol=length(x)+1)
  
  # do Metropolis sampling
  for (i in 1:N)
  {   
    if (i %% (N/100) == 0 & verbose==TRUE)
    {
      pb$tick()
    }
    
    prop <- rmvnorm(n = 1, x, VCOV)
    tp <- target(prop)
    tx <- target(x)
    if (!(tp==0 & tx==0)) # effectively assumes 0/0 = 0 and discards proposal
    {
      if (runif(1) < min(1, tp/tx))
      {
        x <- prop
      }
    }
    
    samples[i,] <- c(x,tx)
  }
  
  return(samples)
}

sample_model_metropolis <- function(model, data, start=NULL, sd = NULL, N.samples=1e4, burnin=0, thin.by=1, verbose=TRUE)
{
  # set dimensionality
  d = ncol(data)
  
  # choose sd vector
  if (is.null(sd))
  {
    sd = as.numeric(0.1*diff(apply(data,2,quantile,c(0.05,0.95))))
  }
  
  normalization = prod(sqrt(2*pi*sd^2))
  
  # choose a starting point
  if (is.null(start))
  {
    start = apply(data, 2, median, na.rm=T)
  }
  
  target <- function(x) { 
    predict(model, newdata=matrix(x,ncol=d)) 
  }
  
  # perform Metropolis sampling
  if (verbose == TRUE)
  {
    cat('Metropolis sampling of model...\n')
  }
  samples = rwmetro(target, N=N.samples, x=start, VCOV=sd*diag(d), verbose=verbose)
  if (verbose == TRUE)
  {
    cat(' done.\n')
  }
  
  if (verbose == TRUE)
  {
    cat('Thinning samples...')
  }
  # remove burn-in samples
  samples = samples[(1+burnin):nrow(samples),]
  
  # thin results to reduce autocorrelation
  samples = samples[seq(1,nrow(samples),by=thin.by),]
  if (verbose == TRUE)
  {
    cat(' done...\n')
  }
  
  dimnames(samples) <- list(NULL, c(names(start),"value"))
  
  # return result
  return(samples)
}

sample_model_rejection <- function(model, range, N.samples, chunksize=1e3, verbose=TRUE, min.value=0, ...) # range should be the output of e.g. apply(data, 2, range) (2 x n matrix)
{
  # determine dimensionality
  d = ncol(range)
  
  if (is.null(dimnames(range)[[2]]))
  {
    dimnames(range)[[2]] <- paste("X",1:d,sep="")
  }
  
  if (verbose==TRUE)
  {
    pb <- progress_bar$new(total = N.samples)
  }
  
  samples = list()
  total_accepted <- 0
  while(total_accepted < N.samples)
  {
    pb$update(total_accepted/N.samples)
    
    # generate random points
    random_points <- data.frame(matrix(data=runif(chunksize*d),ncol=d,nrow=chunksize,dimnames=list(NULL, dimnames(range)[[2]])))
    
    # rescale points
    for (i in 1:d)
    {
      random_points[,i] <- random_points[,i]*(range[2,i] - range[1,i]) + range[1,i]
    }
    
    # predict values at these points
    samples_this = predict(model, newdata=random_points, ...)
    
    # store the new values
    samples = c(samples, list(cbind(random_points, samples_this)))
    
    # store the number of accepted points
    total_accepted = total_accepted + length(which(samples_this>min.value))
  }
  if (verbose==TRUE)
  {
    pb$update(1)
    pb$tick()
  }
  
  samples <- as.matrix(rbindlist(samples))

  dimnames(samples) <- list(NULL, c(dimnames(range)[[2]],"value"))
  
  # cutoff the samples to return the correct number of desired points
  index_cutoff = which(cumsum(samples[,ncol(samples)]>min.value) >= N.samples)[[1]]
  samples <- samples[1:index_cutoff,]
  
  return(samples)
}
