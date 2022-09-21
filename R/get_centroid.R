get_centroid <- function(hv) {
  if (inherits(hv, "Hypervolume")) {
    if (grepl("occupancy", hv@Method)) {
      warning("get_centroid is probably not the best choice for occupancy objects. Please consider get_centroid_weighted instead.")
    }
    result <- apply(hv@RandomPoints, 2, mean)
    names(result) <- dimnames(hv@RandomPoints)[[2]]
    return(result)
  } else if (inherits(hv, "HypervolumeList")) {
    method_list <- unique(unlist(lapply(hv@HVList, function(x) x@Method)))
    if (any(grepl("occupancy", method_list))) {
      warning("get_centroid is probably not the best choice for occupancy objects. Please consider get_centroid_weighted instead.")
    }
    result <- (sapply(hv@HVList, function(x) {
      apply(x@RandomPoints, 2, mean)
    }))
    dimnames(result)[[2]] <- sapply(hv@HVList, function(x) {
      x@Name
    })
    return(t(result))
  } else {
    stop("Wrong class input.")
  }
}
