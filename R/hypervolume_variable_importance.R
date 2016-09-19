hypervolume_variable_importance <- function(hv, verbose=TRUE)
{
  method <- hv@Method
  data <- hv@Data
  params <- hv@Parameters
  density <- hv@PointDensity
  
  hv_others <- rep(NA, ncol(data))
  
  for (i in 1:length(hv_others) )
  {
    if (verbose==TRUE)
    {
      cat('* Repeating hypervolumes for ')
    }
    vars <- setdiff(1:ncol(data),i)

    if (method=="box")
    {
      bw <- params[grep("bandwidth",names(params))]
      names(bw) <- gsub("bandwidth.","",names(bw),fixed=TRUE)

      hv_this <- hypervolume_box(data[,vars], output.density=density, name=NULL, verbose=verbose, bandwidth=bw[vars])
    }
    else if (method=="svm")
    {
      binwidths <- params[grep("expectation.bin.widths",names(params))]
      
      hv_this <- hypervolume_svm(data[,vars], output.density=density, name=NULL, verbose=verbose, svm.nu=params["svm.nu"],svm.gamma=params["svm.gamma"],expectation.num.shifts = params["expectation.num.shifts"], expectation.bin.widths = binwidths)
    }
    else if (method=="gaussian")
    {
      bw <- params[grep("kde.bandwidth",names(params))]
      names(bw) <- gsub("kde.bandwidth.","",names(bw),fixed=TRUE)
      
      binwidths <- params[grep("expectation.bin.widths",names(params))]
      
      hv_this <- hypervolume_gaussian(data[,vars], output.density=density, name=NULL, verbose=verbose, kde.bandwidth=bw[vars],output.threshold=params["output.threshold"], expectation.num.shifts = params["expectation.num.shifts"], expectation.bin.widths = binwidths)
    }
    else
    {
      stop('Importance calculation only possible when hypervolume has @Data and @Method available.')
    }
    
    hv_others[i] <- hv_this@Volume
  }
  
  scaledvalues <- hv@Volume / hv_others
  
  names(scaledvalues) <- colnames(data)
  
  return(scaledvalues)
}
