hypervolume_svm <- function(data, name=NULL, verbose=TRUE, samples.per.point=ceiling((10^(1+ncol(data)))/nrow(data)), range.padding.multiply.interval.amount=0.5, range.padding.add.amount=0, svm.nu=0.01, svm.gamma=0.5, chunksize=1e4)
{
  data <- as.matrix(data)
  if (is.null(dimnames(data)[[2]]))
  {
    dimnames(data)[[2]] <- paste("X",1:ncol(data),sep="")
  }
  
  num.samples = nrow(data) * samples.per.point

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
  svm.model<-e1071::svm(data,
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
  
  # delineate the hyperbox over which the function will be evaluated
  range.box <- padded_range(data, multiply.interval.amount = range.padding.multiply.interval.amount, add.amount = range.padding.add.amount)
  dimnames(range.box) <- list(NULL, dimnames(data)[[2]])
  
  # do rejection sampling over the hyperbox
  hv_svm = hypervolume_general_model(svm.model, name=name, verbose=verbose, range.box = range.box, num.samples=num.samples, chunksize=chunksize, min.value=0)
  hv_svm@Parameters = c(svm.nu=svm.nu, svm.gamma=svm.gamma, range.padding.multiply.interval.amount=range.padding.multiply.interval.amount, range.padding.add.amount=range.padding.add.amount, samples.per.point=samples.per.point)
  hv_svm@Method = 'One-class support vector machine'
  hv_svm@Data = data
  
  #samples <- sample_model_metropolis(svm.model, data, verbose=verbose, Nsamples = num_samples*(1+burnin.fraction), burnin = burnin.fraction*num_samples, thin.by=thin.by)
  
  return(hv_svm)	
}
