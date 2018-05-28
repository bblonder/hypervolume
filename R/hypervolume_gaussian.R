# s is a test point
# means is the centers of each data point
# kernel_sd is the kernel bandwidth (in natural units, not squared, diagonal only)
calculate_density <- function(s, means, kernel_sd, weight=NULL, chunksize=10, verbose=TRUE) {
  
  d = ncol(means)
  n_points = nrow(means)
  density = 0
  
  chunks = ceiling((1:n_points)/chunksize)
  
  num_chunks = max(chunks)
  
  if (verbose==TRUE)
  {
    pb <- progress_bar$new(total = num_chunks)
    pb$tick(0)
  }
  
  for (chunk_this in 1:num_chunks)
  {
    index_vals = (1:n_points)[which(chunks==chunk_this)]
    
    for (i in index_vals) 
    { 
      density = density + weight[i] * dmvnorm(s, mean = means[i, ], sigma = kernel_sd^2 * diag(d))
    }
    
    if (verbose==TRUE) 
    {
      pb$tick()
    }
  }
  
  density = density / n_points
  return(density)
}


hypervolume_gaussian <- function(data, name=NULL, weight=NULL, samples.per.point=ceiling((10^(3+sqrt(ncol(data))))/nrow(data)), kde.bandwidth=estimate_bandwidth(data), sd.count=3, quantile.requested=0.95, quantile.requested.type="probability", chunk.size=1e3, verbose=TRUE, ...)
{
  data = as.matrix(data)
  d = ncol(data)
  np = nrow(data)
  
  if (is.null(weight))
  {
    weight <- rep(1/nrow(data),nrow(data))
  }
  
  if (is.null(dimnames(data)[[2]]))
  {
    dimnames(data)[[2]] <- paste("X",1:ncol(data))
  }
  
  if (sd.count<3)
  {
    warning(sprintf("Values of sd.count (%d) is low.\nRecommended minimum value is 3, with higher values giving better performance.\nBoundaries and volumes may be inaccurate.",sd.count))
  }
  
  # do error check on kde.bandwidth vector
  if (length(kde.bandwidth)==1)
  {
    kde.bandwidth = rep(kde.bandwidth, ncol(data))
  }
  if(ncol(data)!=length(kde.bandwidth))
  {
    stop('data and kde.bandwidth must have same dimensionality')
  }
  
  if (any(kde.bandwidth==0))
  {
    stop('Bandwidth must be non-zero.')
  }
  
  names(kde.bandwidth) <- dimnames(data)[[2]]
  
  if (length(weight)!=nrow(data))
  {
    stop("The length of the weights must be equal to the number of observations.")
  }
  if (abs(sum(weight)-1)>10^-10)
  {
    warning("The sum of the weights must be equal to 1. Normalizing the weights.")
    weight <- weight / sum(weight)
  }
  
  predict_function_gaussian <- function(x)
  {
    return(calculate_density(s=x,means=data,kernel_sd=kde.bandwidth,weight=weight))
  }
  
  # do elliptical sampling
  samples_all = sample_model_ellipsoid(
    predict_function = predict_function_gaussian,
    data = data,
    scales = kde.bandwidth * (sd.count), # distance out to sample
    min.value = 0,
    samples.per.point = samples.per.point, 
    chunk.size=chunk.size, 
    verbose=verbose)
  
  random_points = samples_all$samples[,1:d,drop=FALSE]
  values_accepted = samples_all$samples[,d+1]
  volume <- samples_all$volume
  
  point_density = nrow(random_points) / volume

  hv_gaussian <- new("Hypervolume",
                Data=data,
                Method = 'Gaussian kernel density estimate',
                RandomPoints= random_points,
                PointDensity= point_density,
                Volume= volume,
                Dimensionality=d,
                ValueAtRandomPoints=values_accepted,
                Name=ifelse(is.null(name), "untitled", toString(name)),
                Parameters = list(kde.bandwidth=kde.bandwidth, samples.per.point=samples.per.point, sd.count=sd.count, quantile.requested=quantile.requested, quantile.requested.type=quantile.requested.type))
  
  # apply desired threshold type
  hv_gaussian_thresholded = hypervolume_threshold(hv_gaussian, 
                                                  quantile.requested=quantile.requested,
                                                  quantile.requested.type=quantile.requested.type, 
                                                  uniform.density = FALSE, # keep probability values
                                                  verbose=verbose,
                                                  num.thresholds = 1000, # give higher resolution to obtain the desired quantile
                                                  plot=FALSE, ...)[["HypervolumesThresholded"]]
  hv_gaussian_thresholded@Name = hv_gaussian@Name
  
  return(hv_gaussian_thresholded)	
}