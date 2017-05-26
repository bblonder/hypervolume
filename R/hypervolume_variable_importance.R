hypervolume_variable_importance <- function(hv, verbose=TRUE)
{
  method <- hv@Method
  data <- hv@Data
  params <- hv@Parameters
  
  if (is.null(dimnames(data)[[2]]))
  {
    dimnames(data) <- list(NULL, paste("X",1:ncol(data),sep=""))
  }
  
  hv_others <- rep(NA, ncol(data))
  
  for (i in 1:length(hv_others) )
  {
    if (verbose==TRUE)
    {
      cat(sprintf('* Repeating hypervolumes for dimension %s\n', dimnames(data)[[2]][i]))
    }
    vars <- setdiff(1:ncol(data),i)

    if (method=="Box kernel density estimate")
    {
      hv_this <- hypervolume_box(data[,vars,drop=FALSE], 
                                 samples.per.point=params[["samples.per.point"]], 
                                 name=NULL, 
                                 verbose=verbose, 
                                 kde.bandwidth=params[["kde.bandwidth"]][vars]
                                 )
    }
    else if (method=="One-class support vector machine")
    {
      hv_this <- hypervolume_svm(data[,vars,drop=FALSE], 
                                 samples.per.point=params[["samples.per.point"]], 
                                 name=NULL, 
                                 verbose=verbose, 
                                 svm.nu=params[["svm.nu"]],
                                 svm.gamma=params[["svm.gamma"]]
                                 )
    }
    else if (method=="Gaussian kernel density estimate")
    {
      hv_this <- hypervolume_gaussian(data[,vars,drop=FALSE], 
                                      samples.per.point=params[["samples.per.point"]], 
                                      name=NULL, 
                                      verbose=verbose, 
                                      kde.bandwidth=params[["kde.bandwidth"]][vars], 
                                      sd.count=params[["sd.count"]],
                                      quantile.requested=params[["quantile.requested"]],
                                      quantile.requested.type=params[["quantile.requested.type"]]
                                      )
    }
    else
    {
      stop('Importance calculation only possible for restricted set of hypervolume types and for those with @Data slot.')
    }
    
    hv_others[i] <- hv_this@Volume
  }
  
  scaledvalues <- hv@Volume / hv_others
  
  names(scaledvalues) <- colnames(data)
  
  return(scaledvalues)
}
