## should reference padded_range in the 
hypervolume_general_model <- function(model, name=NULL, verbose=TRUE, data=NULL, range.box=NULL, num.samples=10^(1+ncol(range.box)), chunksize=1e4, min.value=0, ...)
{
  if (!is.null(data) & is.null(range.box))
  {
    range.box <- padded_range(data)
  }
  else if (is.null(data) & is.null(range.box))
  {
    stop('Must specify data or range.box')
  }
    
  
  d = ncol(range.box)
  
  if(is.null(dimnames(range.box)[[2]]))
  {
    dimnames(range.box)[[2]] <- paste("X",1:d,sep="")
  }
  
  # delineate the hyperbox over which the function will be evaluated
  range.box_volume <- prod(apply(range.box, 2, diff))
  
  samples <- sample_model_rejection(model, range=range.box, verbose=verbose, N.samples = num.samples, min.value=min.value, chunksize=chunksize, ...)
  samples_accepted <- samples[samples[,ncol(samples)]>min.value,1:(ncol(samples)-1),drop=FALSE] # last column contains the values
  dimnames(samples_accepted) <- list(NULL, dimnames(range.box)[[2]])
  values_accepted <- samples[samples[,ncol(samples)]>min.value,ncol(samples)]
  
  stopifnot(nrow(samples_accepted) == num.samples)
  
  # Monte Carlo integration step
  point_density <- nrow(samples) / range.box_volume # all the uniformly random points we tried (both rejected & accepted)
  volume <- nrow(samples_accepted) / point_density # only the accepted points
  
  hv <- new("Hypervolume",
                Data=matrix(NA,ncol=d,nrow=1),
                Method = 'Rejection sampling',
                RandomUniformPointsThresholded= samples_accepted,
                PointDensity= point_density,
                Volume= volume,
                Dimensionality=d,
                ProbabilityDensityAtRandomUniformPoints=normalize_probability(values_accepted,point_density),
                Name=ifelse(is.null(name), "untitled", toString(name))) # no parameters
  hv@Parameters = c(num.samples=num.samples, min.value=min.value)
  
  if (!is.null(data))
  {
    hv@Data = data
  }
  
  if (nrow(hv@RandomUniformPointsThresholded) < 10^d)
  {
    warning(sprintf("Hypervolume is represented by a low number of random points (%d) - suggested minimum %d.\nConsider increasing point density to improve accuracy.",nrow(hv@RandomUniformPointsThresholded),10^ncol(data)))
  }
  
  return(hv)	
}
