# choose a probability value equal to the prob density of a single point at 1 s.d. distance
estimate_threshold_gaussian <- function(sd.count=1, bandwidth)
{
  k = length(bandwidth) # dimensionality
  # standardize distance to a given number of s.d.s away
  probs <- rep(NA, length(bandwidth))
  for (i in 1:length(bandwidth))
  {
    probs[i] <- 1/sqrt(2*pi*(2*bandwidth[i])^2) * exp(-1/2*sd.count^2)
  }
  # total probability is product of probability along each axis (no covariance)
  prob_final <- prod(probs)
  return(prob_final)
}

hypervolume_gaussian <- function(data, name=NULL, verbose=TRUE, output.density=10^(ncol(data)), expectation.num.shifts=2, expectation.bin.widths=2*estimate_bandwidth(data), kde.bandwidth=estimate_bandwidth(data)/2, kde.chunksize=1e4, threshold=estimate_threshold_gaussian(sd.count=expectation.num.shifts, bandwidth=kde.bandwidth))
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
  
  # threshold the edges within the adaptive grid
  points_final <- expectation@RandomUniformPointsThresholded[kde.probs > threshold,,drop=FALSE]
  probs_final <- kde.probs[kde.probs > threshold]
  
  numpoints_kde <- nrow(points_final)
  vol_kde <- numpoints_kde / expectation@PointDensity
  
  finalparams <- c(kde.bandwidth, 
                   expectation.num.shifts=expectation.num.shifts, 
                   expectation.bin.widths=expectation.bin.widths,
                   threshold = threshold)

  hv_kde_thresholded <- new("Hypervolume",
                Data=as.matrix(data),
                RandomUniformPointsThresholded= points_final,
                PointDensity= expectation@PointDensity,
                Volume= vol_kde,
                Dimensionality=ncol(data),
                ProbabilityDensityAtRandomUniformPoints= normalize_probability(probs_final, expectation@PointDensity),
                Name=ifelse(is.null(name), "untitled", toString(name)),
                Method = "gaussian",
                Parameters=finalparams)
  
  # apply quantile threshold
  #hv_kde_thresholded <- hypervolume_quantile_threshold(hv_kde, quantile.requested=output.quantile.threshold, plot=FALSE)
  
  if (nrow(hv_kde_thresholded@RandomUniformPointsThresholded) < 10^ncol(data))
  {
    warning(sprintf("Hypervolume is represented by a low number of random points (%d) - suggested minimum %d.\nConsider increasing point density to improve accuracy.",nrow(hv_kde_thresholded@RandomUniformPointsThresholded),10^ncol(data)))
  }
  
  return(hv_kde_thresholded)	
}