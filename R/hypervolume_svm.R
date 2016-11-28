hypervolume_svm <- function(data, name=NULL, verbose=TRUE, output.density=10^ncol(data), expectation.num.shifts=1, expectation.bin.widths=2*estimate_bandwidth(data), svm.nu=0.01, svm.gamma=0.5, svm.chunksize=1e4)
{
  data <- as.matrix(data)
  
  expectation <- expectation_adaptive_box(data, density= output.density, num.shifts=expectation.num.shifts, bin.widths=expectation.bin.widths,verbose=verbose)
  
  if (svm.nu <= 0)
  {
    stop("Parameter svm.nu must be positive")
  }
  
  if (svm.gamma <= 0)
  {
    stop("Parameter svm.gamma must be positive")
  }
  
  if (svm.chunksize <= 100 | svm.chunksize > 1e7)
  {
    stop("Parameter svm.chunksize may be in range that will produce poor performance.")
  }
  if (verbose == TRUE)
  {
    cat('Building support vector machine model')
  }
  svm.model<-e1071::svm(data,
                        y=NULL,
                        type='one-classification',
                        nu= svm.nu,
                        gamma= svm.gamma,
                        scale=TRUE,
                        kernel="radial")	
  if (verbose == TRUE)
  {
    cat('...done\n')
  }
  
  np <- nrow(expectation@RandomUniformPointsThresholded)
  
  if (verbose == TRUE)
  {
    cat(sprintf('Sampling %d random points from support vector machine, %d per chunk.\n',np, svm.chunksize))
  }
  num.samples.completed <- 0
  num.chunks <- ceiling(np/svm.chunksize)
  svm.probs <- vector(mode="list",length=num.chunks)
  for (i in 1:num.chunks)
  {
    if (verbose == TRUE)
    {
      cat(sprintf('Running chunk %d / %d\n', i, num.chunks))
    }
    num.samples.to.take <- min(svm.chunksize, np - num.samples.completed)
    
    inputdata.this <- expectation@RandomUniformPointsThresholded[num.samples.completed:(num.samples.completed+num.samples.to.take-1),]
    
    svm.pred.this <- predict(svm.model,inputdata.this)
    svm.pred.this.pos <- inputdata.this[svm.pred.this==TRUE,]
    
    svm.probs[[i]] <- svm.pred.this.pos
    num.samples.completed <- num.samples.completed + num.samples.to.take
  }
  svm.probs <- do.call("rbind",svm.probs)
  dimnames(svm.probs) <- list(NULL,dimnames(data)[[2]])
  if (verbose == TRUE)
  {
    cat('done.\n')	
  }
  
  numpoints_svm <- nrow(svm.probs)
  vol_svm <- numpoints_svm / expectation@PointDensity
  
  hv_svm <- new("Hypervolume",
                Data=data,
                Method = "svm",
                RandomUniformPointsThresholded= svm.probs,
                PointDensity= expectation@PointDensity,
                Volume= vol_svm,
                Dimensionality=ncol(svm.probs),
                ProbabilityDensityAtRandomUniformPoints=normalize_probability(rep(1,nrow(svm.probs)),expectation@PointDensity),
                Name=ifelse(is.null(name), "untitled", toString(name)),
                Parameters=c(svm.nu=svm.nu, svm.gamma=svm.gamma, expectation.num.shifts=expectation.num.shifts, expectation.bin.widths=expectation.bin.widths))
  
  if (nrow(hv_svm@RandomUniformPointsThresholded) < 10^ncol(data))
  {
    warning(sprintf("Hypervolume is represented by a low number of random points (%d) - suggested minimum %d.\nConsider increasing point density to improve accuracy.",nrow(hv_svm@RandomUniformPointsThresholded),10^ncol(data)))
  }
  
  return(hv_svm)	
}
