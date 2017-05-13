
calculate_density <- function(s, means, kernel_sd, chunksize=10, verbose=TRUE) {

  d = ncol(means)
  n_points = nrow(means)
  density = 0
  
  chunks = ceiling((1:n_points)/chunksize)
  
  num_chunks = max(chunks)
  
  if (verbose==TRUE)
  {
	  pb <- progress_bar$new(total = num_chunks)
	  pb$tick()
  }
  
  for (chunk_this in 1:num_chunks)
  {
  	index_vals = (1:n_points)[which(chunks==chunk_this)]
  	
  	for (i in index_vals) 
  	{ 
    	density = dmvnorm(s, mean = means[i, ], sigma = kernel_sd^2 * diag(d)) + density
  	}
  	
  	if (verbose==TRUE) 
  	{
  		pb$tick()
  	}
  }

  density = density / n_points
  return(density)
}

importance_weights = function(q, min_density) {
  return((q > min_density) / q)
}

hypervolume_gaussian <- function(data, kde.bandwidth=estimate_bandwidth(data)/2, samples.per.point=10^ncol(data), threshold.sd.count=3, name=NULL, verbose=TRUE)
{
	data = as.matrix(data)
	d = ncol(data)
	np = nrow(data)
	
	# do error check on kde.bandwidth vector
	if (length(kde.bandwidth)==1)
	{
		kde.bandwidth = rep(kde.bandwidth, ncol(data))
	}
	if(ncol(data)!=length(kde.bandwidth))
	{
		Stop('Data and kde.bandwidth must have same dimensionality')
	}
	
	if (any(kde.bandwidth==0))
	  {
	    stop('Bandwidth must be non-zero.')
	  }
	  
  names(kde.bandwidth) <- paste("kde.bandwidth",dimnames(data)[[2]],sep=".")
  
	
	# generate uniform random points around each data point out to a certain number of standard deviations (add one to make sure we 'see' enough of the shape)
	s = lapply(1:nrow(data),
           function(x){
             result = rmvnorm(samples.per.point, data[x, ], kde.bandwidth^2 * diag(d))
             return(data.frame(result))
             
             })
	s = as.matrix(rbindlist(s))

	# calculate probability density at each point
	if (verbose==TRUE)
	{
		cat('Calculating density...\n')
	}
	density = calculate_density(s, data, kde.bandwidth, verbose)
	if (verbose==TRUE)
	{
		cat('done.\n')
	}
	
	# calculate the volume for points above a minimum threshold value corresponding to the probability value a certain distance from an isolated point
	min_density_threshold = estimate_threshold_gaussian(sd.count=threshold.sd.count, kde.bandwidth=kde.bandwidth)
	min_density = min_density_threshold
	# do importance sampling
	if (verbose==TRUE)
	{
	  cat('Doing importance sampling...\n')
	}
	w = importance_weights(density, min_density)
	# calculate volume
	volume = mean(w) # Estimated volume
	if (verbose==TRUE)
	{
		cat('done.\n')
	}

	# resample to a uniform density (possibly repeating points)
	if (verbose==TRUE)
	{
		cat('Resampling to uniform density...\n')
	}
	inside_sample_rows = sample.int(nrow(s), 
                          samples.per.point * np,
                          replace = TRUE, 
                          prob = w)
	if (verbose==TRUE)
	{
		cat('done.\n')
	}
   
	# apply threshold
  random_points_thresholded = s[inside_sample_rows, ,drop=FALSE]
  # some introduction of very small noise is needed to avoid problems with duplicated points
  random_points_thresholded = random_points_thresholded + matrix(rnorm(n=prod(dim(random_points_thresholded)),mean=0,sd=mean(kde.bandwidth)*1e-6),nrow=nrow(random_points_thresholded),ncol=ncol(random_points_thresholded))
  probability_thresholded = density[inside_sample_rows]
  
  point_density = nrow(random_points_thresholded) / volume 	
    
  finalparams <- c(kde.bandwidth, 
                   samples.per.point = samples.per.point, 
                   threshold.sd.count = threshold.sd.count)

  hv_kde_thresholded <- new("Hypervolume",
                Data=as.matrix(data),
                RandomUniformPointsThresholded= random_points_thresholded,
                PointDensity= point_density,
                Volume= volume,
                Dimensionality=ncol(data),
                ProbabilityDensityAtRandomUniformPoints= normalize_probability(probability_thresholded, point_density),
                Name=ifelse(is.null(name), "untitled", toString(name)),
                Method = "Gaussian kernel density estimate",
                Parameters=finalparams)   
  
  if (nrow(hv_kde_thresholded@RandomUniformPointsThresholded) < 10^ncol(data))
  {
    warning(sprintf("Hypervolume is represented by a low number of random points (%d) - suggested minimum %d.\nConsider increasing point density to improve accuracy.",nrow(hv_kde_thresholded@RandomUniformPointsThresholded),10^ncol(data)))
  }
                
  return(hv_kde_thresholded) 
}



