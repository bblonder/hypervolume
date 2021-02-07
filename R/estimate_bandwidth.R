estimate_bandwidth <- function(data,method="silverman",value=NULL)
{
  # remove missing cases
  data <- na.omit(data)
  npoints = nrow(data)
  ndim <- ncol(data)
  
  bandwidth_final = NULL
  
  if (method=='fixed')
  {
    if (is.null(value))
    {
      stop('When using fixed bandwidth, value must be non-NULL')
    }
    else
    {
      if (length(value)==1)
      {
        bandwidth_final = rep(value, ncol(data))
      }
      else if (length(value)==ndim)
      {
        bandwidth_final = value
      }
      else
      {
        stop('value must have either length=1 or =number of columns of input data')
      }
    }
  }
  else if (method=="silverman")
  {
    stdev = apply(data,2,sd)
    
    message("Note that the formula used for the Silverman estimator differs in version 3 compared to prior versions of this package.\nUse method=\'silverman-1d\' to replicate prior behavior.")
    
    bandwidth_final = (4/(ndim+2))^(1/(ndim+4)) * npoints^(-1/(ndim+4))*stdev
  }
  else if (method=="silverman-1d")
  {
    stdev = apply(data,2,sd)
    
    bandwidth_final = 1.06 * stdev * npoints ^ (-1/5)
  }
  else
  {
    if (ndim <= 6)
    {
     if (method=="plug-in")
     {
       bw_plugin <- ks::Hpi.diag(data, nstage=2, pilot="samse",pre="scale")
       
       # convert estimated variances to standard deviations
       bandwidth_final <- sqrt(diag(bw_plugin))
     }
     else if (method=="cross-validation")
     {
       bw_crossvalidation <- ks::Hscv.diag(data, pilot="samse", pre="scale")
       
       # convert estimated variances to standard deviations
       bandwidth_final <- sqrt(diag(bw_crossvalidation))
     }
     else
     {
       stop("Argument 'method' not recognized")
     }
    }
    else
    {
      stop("Non-Silverman bandwidth estimators not available for n>6 dimensions.")
    }
  }
  
  # set attribute on method
  attr(bandwidth_final,"method") <- method
  
  return(bandwidth_final)
}