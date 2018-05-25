hypervolume_box <- function(data, name=NULL, verbose=TRUE, samples.per.point=ceiling((10^(3+sqrt(ncol(data))))/nrow(data)), kde.bandwidth=2*estimate_bandwidth(data), tree.chunksize=1e4)
{
  data <- as.matrix(data)
  if (is.null(dimnames(data)[[2]]))
  {
    dimnames(data)[[2]] <- paste("X",1:ncol(data))
  }
  
  dim = ncol(data)
  np = nrow(data)
  
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
  
  names(kde.bandwidth) <- dimnames(data)[[2]]
  
  # double the bandwidth as for other functions it is interpreted as a box half-width
  # but in this function is interpreted as a box full-width.
  #kde.bandwidth = 2*kde.bandwidth
  
  # figure out the hypervolume of one kernel and the random point density within it
  # note this used to have a factor of 2 in previous versions
  hyperbox_volume = prod(2*kde.bandwidth) # hyperbox is in both + and - dimensions
  
  point_density = ceiling(samples.per.point / hyperbox_volume)

  # generate random point cloud around each data point
  num.samples.completed <- 0
  num.chunks <- ceiling(np/tree.chunksize)
  data_points <- vector(mode="list",length=num.chunks) # the randomuniform points
  points.all <- vector(mode="list",length=num.chunks) # the probability at each point
  
  if (verbose==TRUE)
  {
    pb <- progress_bar$new(total = num.chunks)
    pb$tick(0)
  }
  
  for (i in 1:num.chunks)
  {
    if (verbose == TRUE)
    {
      pb$tick()
    }
    num.samples.to.take <- min(tree.chunksize, np - num.samples.completed)
    random_points = 2*(matrix(runif(samples.per.point*num.samples.to.take*dim,min=0,max=1),nrow=samples.per.point*num.samples.to.take, ncol=dim) - 0.5) * repmat(t(as.matrix(kde.bandwidth)), samples.per.point*num.samples.to.take, 1)
    offset = repmat(as.matrix(data[(num.samples.completed+1):(num.samples.completed+num.samples.to.take),,drop=FALSE]), samples.per.point, 1)  
    data_points[[i]] = random_points + offset 
    
    # determine the probability density at each random point
    points.all[[i]] = evalfrectangular(data=data, bandwidth=kde.bandwidth, points=data_points[[i]])
    
    num.samples.completed <- num.samples.completed + num.samples.to.take
  }
  # put everything back together
  if (verbose == TRUE)
  {
    cat('Binding random points...')
  }  
  point_counts <- do.call("c",points.all)
  data_points <- do.call("rbind",data_points)
  if (verbose == TRUE)
  {
    cat(' done.\n')
  }
  
  
  # infer the total volume based on the random point counts and the quantile (now set to zero)
  # note that it may not be possible to achieve the exact quantile specified
  vc = volume_calculation(point_counts, point_density, quantile=0, verbose)
  disjunctfactor <- vc$final_volume / (nrow(as.matrix(data)) * prod(2*kde.bandwidth))
  
  
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
  weightedsample = sample(x=1:nrow(point_counts_final),size=ceiling(vc$final_volume * point_density),replace=TRUE,prob=invweight)
  options(warn=ow)
  # keep only unique points
  weightedsample = unique(weightedsample)

  # prepare object for output
  points_uniform_final = as.matrix(point_counts_final[weightedsample,1:(ncol(point_counts_final)-1)])
  dimnames(points_uniform_final) <- list(NULL, dimnames(data)[[2]])
  density_uniform_final = point_counts_final[weightedsample,ncol(point_counts_final)]
  point_density_final = nrow(points_uniform_final) / vc$final_volume
  
  hv_box = new("Hypervolume", Name=ifelse(is.null(name), "untitled", toString(name)))
  hv_box@Method = "Box kernel density estimate"
  hv_box@Data = as.matrix(data)
  hv_box@Dimensionality = dim
  hv_box@Volume = vc$final_volume
  hv_box@PointDensity = point_density_final
  hv_box@Parameters = list(kde.bandwidth=kde.bandwidth, samples.per.point=samples.per.point)
  hv_box@RandomPoints = as.matrix(points_uniform_final);  
  hv_box@ValueAtRandomPoints = density_uniform_final
  
  return(hv_box)  

}
