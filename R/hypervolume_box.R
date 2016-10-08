hypervolume_box <- function(data, name=NULL, verbose=TRUE, output.density=10^(ncol(data)), repsperpoint=NULL, bandwidth=estimate_bandwidth(data), kdtree.chunksize=1e4)
{
  data <- as.matrix(data)
  
  dim = ncol(data)
  np = nrow(data)
  
  if (length(bandwidth) == 1)
  {
    bandwidth <- rep(bandwidth, ncol(data))
  }  
    
  if (length(bandwidth) != ncol(data))
  {
    stop('Input bandwidth vector not same length as dimensionality of dataset.')
  }
  
  if (any(bandwidth==0))
  {
    stop('Bandwidth must be non-zero.')
  }
  
  names(bandwidth) <- paste("bandwidth",dimnames(data)[[2]],sep=".")
  
  # figure out the hypervolume of one kernel and the random point density within it
  hyperbox_volume = prod(2*bandwidth) # hyperbox is in both + and - dimensions
  
  if (is.null(repsperpoint))
  {
    repsperpoint = ceiling(output.density * hyperbox_volume)
  }
  else
  {
    output.density = repsperpoint / hyperbox_volume
    warning('Argument repsperpoint is deprecated - please use output.density instead.\n')
  }
  repsperpoint = floor(repsperpoint)
  point_density = floor(output.density) # standardize argument names
  cat(sprintf('Setting repsperpoint=%.0f, output.density=%.0f.\n', repsperpoint, output.density))
  
  # generate random point cloud around each data point
  num.samples.completed <- 0
  num.chunks <- ceiling(np/kdtree.chunksize)
  data_points <- vector(mode="list",length=num.chunks) # the randomuniform points
  points.all <- vector(mode="list",length=num.chunks) # the probability at each point
  for (i in 1:num.chunks)
  {
    if (verbose == TRUE)
    {
      cat(sprintf('Running chunk %d / %d\n', i, num.chunks))
    }
    num.samples.to.take <- min(kdtree.chunksize, np - num.samples.completed)
    random_points = 2*(matrix(runif(repsperpoint*num.samples.to.take*dim,min=0,max=1),nrow=repsperpoint*num.samples.to.take, ncol=dim) - 0.5) * repmat(t(as.matrix(bandwidth)), repsperpoint*num.samples.to.take, 1)
    offset = repmat(as.matrix(data[(num.samples.completed+1):(num.samples.completed+num.samples.to.take),]), repsperpoint, 1)  
    data_points[[i]] = random_points + offset 
    
    # determine the probability density at each random point
    if (verbose == TRUE)
    {
      cat('Evaluating probability density...\n')
    }
    points.all[[i]] = evalfrectangular(data, bandwidth, data_points[[i]],verbose=verbose)
    if (verbose == TRUE)
    {
      cat('Finished evaluating probability density.\n')
    }
    
    num.samples.completed <- num.samples.completed + num.samples.to.take
  }
  # put everything back together
  point_counts <- do.call("c",points.all)
  data_points <- do.call("rbind",data_points)
  if (verbose == TRUE)
  {
    cat('...done.\n')
  }
  
  
  # infer the total volume based on the random point counts and the quantile (now set to zero)
  # note that it may not be possible to achieve the exact quantile specified
  if (verbose == TRUE)
  {
    cat('Doing volume calculation...\n')
  }
  vc = volume_calculation(point_counts, point_density, quantile=0, verbose)
  if (verbose == TRUE)
  {
    cat('done.\n')
  }  
  disjunctfactor <- vc$final_volume / (nrow(as.matrix(data)) * prod(2*bandwidth))
  
  
  # see if the points are disjunct (if the total volume equals the summed volume of the hyperbox around each point)
  if (disjunctfactor > 0.9)
  {
    warning(sprintf('Disjunct factor: %.2f\nRatio of inferred volume to summed volume of hyperbox kernels around each data point is near unity. Points are likely disjunct. Consider increasing the bandwidth.',disjunctfactor))
  }
  # threshold the random points to include only those above the chosen quantile threshold
  point_counts_final = cbind(data_points, point_counts)[point_counts >= vc$index,]
  

  # downweight by the number of times the point intersect and standardize the point density
  invweight = 1 / point_counts_final[,ncol(point_counts_final)]
  ow <- getOption('warn')
  options(warn=-1)
  weightedsample = sample(x=1:nrow(point_counts_final),size=floor(vc$final_volume * point_density),replace=TRUE,prob=invweight)
  options(warn=ow)
  # keep only unique points
  weightedsample = unique(weightedsample)

  # prepare object for output
  points_uniform_final = as.matrix(point_counts_final[weightedsample,1:(ncol(point_counts_final)-1)])
  dimnames(points_uniform_final) <- list(NULL, dimnames(data)[[2]])
  density_uniform_final = point_counts_final[weightedsample,ncol(point_counts_final)]
  point_density_final = nrow(points_uniform_final) / vc$final_volume
  
  hv_box = new("Hypervolume", Name=ifelse(is.null(name), "untitled", toString(name)))
  hv_box@Method = "box"
  hv_box@Data = as.matrix(data)
  hv_box@Dimensionality = dim
  hv_box@Volume = vc$final_volume
  hv_box@PointDensity = point_density_final
  hv_box@Parameters = c(bandwidth, RepsPerPoint=repsperpoint)
  hv_box@RandomUniformPointsThresholded = as.matrix(points_uniform_final);  
  hv_box@ProbabilityDensityAtRandomUniformPoints = normalize_probability(density_uniform_final, point_density_final)

  if (nrow(hv_box@RandomUniformPointsThresholded) < 10^ncol(data))
  {
    warning(sprintf("Hypervolume is represented by a low number of random points (%d) - suggested minimum %d.\nConsider increasing point density to improve accuracy.",nrow(hv_box@RandomUniformPointsThresholded),10^ncol(data)))
  }
  
  return(hv_box)  

}
