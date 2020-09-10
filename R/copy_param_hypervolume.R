copy_param_hypervolume <- function(hv, data, name = NULL) {
  if(hv@Method == 'One-class support vector machine') {
    return(hypervolume_svm(data, name,
                           svm.nu = hv@Parameters$svm.nu,
                           svm.gamma = hv@Parameters$svm.gamma,
                           samples.per.point = hv@Parameters$samples.per.point,
                           verbose = FALSE))
  } else if(hv@Method == 'Box kernel density estimate') {
    return(hypervolume_box(data, name,
                           samples.per.point = hv@Parameters$samples.per.point,
                           verbose = FALSE))
  } else if(hv@Method == 'Gaussian kernel density estimate') {
    return(hypervolume_gaussian(data, name,
                                samples.per.point = hv@Parameters$samples.per.point,
                                sd.count = hv@Parameters$sd.count,
                                quantile.requested = hv@Parameters$quantile.requested,
                                quantile.requested.type = hv@Parameters$quantile.requested.type,
                                verbose = FALSE))
  } else {
    stop("Hypervolume does not have a valid construction method")
  }
}