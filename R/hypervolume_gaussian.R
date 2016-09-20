mvnorm <- function(x, mu, Sigma, k=length(x),fast=FALSE,diagonalSigma=FALSE) # fast does quadratic approximation
{
  xminusmu <- as.numeric(x - mu)
  if (diagonalSigma==TRUE)
  {
    Sigmainv = 1/Sigma
    Sigmadet = Sigma[diag(1,nrow(Sigma),ncol(Sigma))]	
  }
  else
  {
    Sigmainv <- solve(Sigma)
    Sigmadet <- det(Sigma)
  }
  
  argument <- (-1/2) * (t(xminusmu) %*% Sigmainv %*% (xminusmu))  
  
  
  if (fast==TRUE)
  {
    exparg <- 1 + argument
    if (exparg < 0)
    {
      exparg <- 0
    }
  }
  else
  {
    exparg <- exp(argument)
  }
  
  prefactor=1/sqrt(((2*pi)^k)*Sigmadet)
  
  return(as.numeric(exparg*prefactor))
}

# choose a probability value equal to the prob density of a single point at 1 s.d. distance in dimensionless coordinates
estimate_threshold_gaussian <- function(sd.count=1, dim)
{
  mvnorm(rep(sd.count,dim),rep(0,dim),Sigma=diag(1,dim,dim))
}


hypervolume_gaussian <- function(data, name=NULL, verbose=TRUE, output.density=10^(ncol(data)), expectation.num.shifts=1, expectation.bin.widths=2*estimate_bandwidth(data), kde.bandwidth=estimate_bandwidth(data)/2, kde.chunksize=1e4, output.threshold=estimate_threshold_gaussian(sd.count=1,dim=ncol(data)))
{
  data <- as.matrix(data)
  
  if (length(kde.bandwidth) == 1)
  {
    kde.bandwidth <- rep(kde.bandwidth, ncol(data))
  }  
  
  if (length(kde.bandwidth) != ncol(data))
  {
    stop('Input bandwidth vector not same length as dimensionality of dataset.')
  }
  
  if (any(kde.bandwidth==0))
  {
    stop('Bandwidth must be non-zero.')
  }
  
  names(kde.bandwidth) <- paste("kde.bandwidth",dimnames(data)[[2]],sep=".")
  
  if (verbose == TRUE)
  {
    cat('Generating adaptive grid...\n')
  }
  print(data)
  print(expectation.num.shifts)
  print(expectation.bin.widths)
  
  expectation <- expectation_adaptive_box(data, density= output.density, num.shifts=expectation.num.shifts, bin.widths=expectation.bin.widths)
  if (verbose == TRUE)
  {
    cat('Done...\n')
  }
  
  if (ncol(data)==1) # kde package defined bandwidth as s.d. in 1 dim and var in >= 2 dim...
  {
    Hmatrix <- diag(1,nrow=ncol(data),ncol=ncol(data))*(kde.bandwidth)
  }
  else
  {
    Hmatrix <- diag(1,nrow=ncol(data),ncol=ncol(data))*(kde.bandwidth)^2
  }
  
  np <- nrow(expectation@RandomUniformPointsThresholded)
  if (verbose == TRUE)
  {
    cat(sprintf('Sampling %d random points from kernel density estimate, %d per chunk...\n',np, kde.chunksize))
  }
  num.samples.completed <- 0
  num.chunks <- ceiling(np/kde.chunksize)
  kde.probs <- vector(mode="list",length=num.chunks)
  for (i in 1:num.chunks)
  {
    if (verbose == TRUE)
    {
      cat(sprintf('Running chunk %d / %d\n', i, num.chunks))
    }
    num.samples.to.take <- min(kde.chunksize, np - num.samples.completed)
    kde.probs.this <- ks::kde(x=data, 
                              H=Hmatrix,
                              eval.points= expectation@RandomUniformPointsThresholded[num.samples.completed:(num.samples.completed+num.samples.to.take-1),],
                              verbose=TRUE)$estimate
    kde.probs[[i]] <- kde.probs.this
    num.samples.completed <- num.samples.completed + num.samples.to.take
  }
  kde.probs <- do.call("c",kde.probs)
  if (verbose == TRUE)
  {
    cat('...done.\n')
  }
  points_thresholded <- expectation@RandomUniformPointsThresholded[kde.probs > output.threshold, , drop=F]
  probs_thresholded <- kde.probs[kde.probs > output.threshold]
  
  numpoints_kde <- nrow(points_thresholded)
  vol_kde <- numpoints_kde / expectation@PointDensity
  
  finalparams <- c(kde.bandwidth,output.threshold=output.threshold, expectation.num.shifts=expectation.num.shifts, expectation.bin.widths=expectation.bin.widths)

  hv_kde <- new("Hypervolume",
                Data=as.matrix(data),
                RandomUniformPointsThresholded= points_thresholded,
                PointDensity= expectation@PointDensity,
                Volume= vol_kde,
                Dimensionality=ncol(data),
                ProbabilityDensityAtRandomUniformPoints= normalize_probability(probs_thresholded, expectation@PointDensity),
                Name=ifelse(is.null(name), "untitled", toString(name)),
                Method = "gaussian",
                Parameters=finalparams)
  
  if (nrow(hv_kde@RandomUniformPointsThresholded) < 10^ncol(data))
  {
    warning(sprintf("Hypervolume is represented by a low number of random points (%d) - suggested minimum %d.\nConsider increasing point density to improve accuracy.",nrow(hv_kde@RandomUniformPointsThresholded),10^ncol(data)))
  }
  
  return(hv_kde)	
}