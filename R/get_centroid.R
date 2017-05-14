get_centroid <- function(hv)
{
  if (class(hv)=="Hypervolume")
  {
    result = apply(hv@RandomUniformPointsThresholded,2,mean)
    names(result) <- dimnames(hv@RandomUniformPointsThresholded)[[2]]
    return(result)
  }
  else if (class(hv)=="HypervolumeList")
  {
    result = (sapply(hv@HVList, function(x) { apply(x@RandomUniformPointsThresholded,2,mean) }))
    dimnames(result)[[2]] <- sapply(hv@HVList, function(x) {x@Name})
    return(t(result))
  }
  else
  {
    stop('Wrong class input.')
  }
}