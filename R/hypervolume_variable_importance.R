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
      bw <- params[grep("kde.bandwidth",names(params),fixed=TRUE)]
      names(bw) <- gsub("kde.bandwidth.","",names(bw),fixed=TRUE)

      hv_this <- hypervolume_box(data[,vars,drop=FALSE], samples.per.point=params["samples.per.point"], name=NULL, verbose=verbose, kde.bandwidth=bw[vars])
    }
    else if (method=="One-class support vector machine")
    {
      hv_this <- hypervolume_svm(data[,vars,drop=FALSE], samples.per.point=params["samples.per.point"], name=NULL, verbose=verbose, svm.nu=params["svm.nu"],svm.gamma=params["svm.gamma"], range.padding.multiply.interval.amount=params["range.padding.multiply.interval.amount"], range.padding.add.amount=params["range.padding.add.amount"])
    }
    else if (method=="Gaussian kernel density estimate")
    {
      bw <- params[grep("kde.bandwidth",names(params))]
      names(bw) <- gsub("kde.bandwidth.","",names(bw),fixed=TRUE)

      hv_this <- hypervolume_gaussian(data[,vars,drop=FALSE], samples.per.point=params["samples.per.point"], name=NULL, verbose=verbose, kde.bandwidth=bw[vars], threshold.sd.count=params["threshold.sd.count"])
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
