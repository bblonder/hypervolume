copy_param_hypervolume <- function(hv, data, name = NULL) {
  if(hv@Method == 'One-class support vector machine') {
    return(hypervolume_svm(data, name,
                           svm.nu = hv@Parameters$svm.nu,
                           svm.gamma = hv@Parameters$svm.gamma,
                           samples.per.point = hv@Parameters$samples.per.point,
                           verbose = FALSE))
  } else if(hv@Method == 'Box kernel density estimate') {
    method = attr(hv@Parameters$kde.bandwidth, "method")
    if(!(method %in% c("fixed", "silverman", "silverman-1d", "plug-in", "cross-validation"))) {
      stop("Hypervolume does not have valid bandwidth estimation method")
    }
    if(method == "fixed") {
      value = hv@Parameters$kde.bandwidth
    } else {
      value = NULL
    }
    return(hypervolume_box(data, name,
                           samples.per.point = hv@Parameters$samples.per.point,
                           kde.bandwidth = estimate_bandwidth(data, method = method, value = value),
                           verbose = FALSE))
  } else if(hv@Method == 'Gaussian kernel density estimate') {
    method = attr(hv@Parameters$kde.bandwidth, "method")
    if(!(method %in% c("fixed", "silverman", "silverman-1d", "plug-in", "cross-validation"))) {
      stop("Hypervolume does not have valid bandwidth estimation method")
    }
    if(method == "fixed") {
      value = hv@Parameters$kde.bandwidth
    } else {
      value = NULL
    }
    return(hypervolume_gaussian(data, name,
                                samples.per.point = hv@Parameters$samples.per.point,
                                kde.bandwidth = estimate_bandwidth(data, method = method, value = value),
                                sd.count = hv@Parameters$sd.count,
                                quantile.requested = hv@Parameters$quantile.requested,
                                quantile.requested.type = hv@Parameters$quantile.requested.type,
                                verbose = FALSE))
  } else {
    stop("Hypervolume does not have a valid construction method")
  }
}
