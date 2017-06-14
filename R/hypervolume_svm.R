hypervolume_svm <- function(data, name=NULL, samples.per.point=ceiling((10^(3+sqrt(ncol(data))))/nrow(data)), svm.nu=0.01, svm.gamma=0.5, scale.factor=1, chunk.size=1e3, verbose=TRUE)
{
  data <- as.matrix(data)
  d <- ncol(data)
  if (is.null(dimnames(data)[[2]]))
  {
    dimnames(data) <- list(NULL,paste("X",1:d,sep=""))
  }
  
  if (svm.nu <= 0)
  {
    stop("Parameter svm.nu must be >0.")
  }
  
  if (svm.gamma <= 0)
  {
    stop("Parameter svm.gamma must be >0.")
  }
  
  if (verbose == TRUE)
  {
    cat('Building support vector machine model...')
  }
  svm_this <- e1071::svm(data,
                         y=NULL,
                         type='one-classification',
                         nu= svm.nu,
                         gamma= svm.gamma,
                         scale=TRUE,
                         kernel="radial")	
  if (verbose == TRUE)
  {
    cat(' done\n')
  }
  
  # Assuming all of the support vectors were right on top of each other,
  # how far could we move away from that point before we left the interior
  # of the hypervolume? 
  # `solve b * exp(-g d_2) == r for d_2` into Wolfram Alpha, where
  #    * "b" is the sum of the SVM coefficients
  #    * "g" is the SVM kernel bandwidth (gamma)
  #    * "r" is the minimum kernel density of the hypervolume (rho)
  #    * "d_2" is squared distance  
  squared_scaled_dist = log(sum(svm_this$coefs) / svm_this$rho) / svm.gamma
  scales = scale.factor * apply(data, 2, sd) * sqrt(squared_scaled_dist)  
  
  predict_function_svm <- function(x)
  {
    return(predict(svm_this, x))
  }  
  
  # do elliptical sampling over the hyperbox
  samples_all = sample_model_ellipsoid(
    predict_function = predict_function_svm,
    data = data[svm_this$index,,drop=FALSE],
    scales = scales, 
    min.value = 0, # because output is binary
    samples.per.point = samples.per.point, 
    chunk.size=chunk.size, 
    verbose=verbose)
  
  random_points = samples_all$samples[,1:d]
  values_accepted <- samples_all$samples[,d+1]
  volume <- samples_all$volume
  point_density = nrow(random_points) / volume
  
  hv_svm <- new("Hypervolume",
                Data=data,
                Method = 'One-class support vector machine',
                RandomPoints= random_points,
                PointDensity= point_density,
                Volume= volume,
                Dimensionality=d,
                ValueAtRandomPoints=values_accepted,
                Name=ifelse(is.null(name), "untitled", toString(name)),
                Parameters = list(svm.nu=svm.nu, svm.gamma=svm.gamma, samples.per.point=samples.per.point))
  
  return(hv_svm)	
}