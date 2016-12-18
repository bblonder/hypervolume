pad_range <- function(vals, frac=0)
{
  v_min <- min(vals,na.rm=T)
  v_max <- max(vals,na.rm=T)
  
  range <- v_max - v_min
  
  return(c(v_min - range*frac, v_max + range*frac))
}

hypervolume_other <- function(model, range, data=NULL, data_pad_fraction=0, name=NULL, verbose=TRUE, output.density=10^ncol(data), threshold_initial=0, normalize=TRUE, chunksize=1e4, check_memory=TRUE, ...)
{
  if (!is.null(data))
  {
    data <- as.matrix(data)
    range <- apply(as.matrix(data), 2, pad_range, data_pad_fraction)
  }
  else
  {
    data <- matrix(NA,nrow=1,ncol=ncol(range))
  }
  
  box_volume <- prod(diff(range))
  np <- ceiling(box_volume * output.density)
  ndim <- ncol(range)
  
  if (check_memory==TRUE)
  {
    stop(sprintf('Analysis will require function evaluation at %d points and storage of %d values. To continue set check_memory=FALSE.', np, np*ndim))
  }
  
  rp <- matrix(runif(np*ndim),nrow=np,ncol=ndim, dimnames=list(NULL,dimnames(range)[[2]]))
  for (i in 1:ncol(rp))
  {
    rp[,i] <- rp[,i] * (range[2,i] - range[1,i]) + range[1,i]
  }
  
  if (chunksize <= 100 | chunksize > 1e7)
  {
    stop("Parameter chunksize may be in range that will produce poor performance.")
  }
  if (verbose == TRUE)
  {
    cat('Evaluating model')
  }
  
  if (verbose == TRUE)
  {
    cat(sprintf('Sampling %d random points from padded range hyperbox, %d per chunk.\n',np, chunksize))
  }
  num.samples.completed <- 0
  num.chunks <- ceiling(np/chunksize)

  model.probs <- rep(NA, np)
  for (i in 1:num.chunks)
  {
    if (verbose == TRUE)
    {
      cat(sprintf('Running chunk %d / %d\n', i, num.chunks))
    }
    num.samples.to.take <- min(chunksize, np - num.samples.completed)
    
    inputdata.this <- rp[num.samples.completed:(num.samples.completed+num.samples.to.take-1),]
    
    model.probs[num.samples.completed:(num.samples.completed+num.samples.to.take-1)] <- predict(model,type="response",as.data.frame(inputdata.this), ...)

    num.samples.completed <- num.samples.completed + num.samples.to.take
  }
  if (verbose == TRUE)
  {
    cat('done.\n')	
  }

  ids.tokeep <- which(model.probs > threshold_initial)
  rp_thresholded <- rp[ids.tokeep,]
  probs_thresholded <- model.probs[ids.tokeep]
  
  if (normalize==TRUE)
  {
    probs_thresholded <- normalize_probability(probs_thresholded, output.density)
  }
  
  numpoints_model <- nrow(probs_thresholded)
  vol_model <- numpoints_model / output.density
  
  hv_other <- new("Hypervolume",
                Data=data,
                Method = "other",
                RandomUniformPointsThresholded=rp_thresholded,
                PointDensity= output.density,
                Volume= vol_model,
                Dimensionality=ndim,
                ProbabilityDensityAtRandomUniformPoints=probs_thresholded,
                Name=ifelse(is.null(name), "untitled", toString(name)),
                Parameters=c(range.lo=range[1,], range.hi=range[2,], threshold_initial=threshold_initial))
  
  if (nrow(hv_other@RandomUniformPointsThresholded) < 10^ncol(range))
  {
    warning(sprintf("Hypervolume is represented by a low number of random points (%d) - suggested minimum %d.\nConsider increasing point density to improve accuracy.",nrow(hv_other@RandomUniformPointsThresholded),10^ncol(range)))
  }
  
  return(hv_other)	
}
