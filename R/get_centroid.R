get_centroid <- function(hv)
{
  if (class(hv)=="Hypervolume")
  {
    result = apply(hv@RandomPoints,2,mean)
    names(result) <- dimnames(hv@RandomPoints)[[2]]
    return(result)
  }
  else if (class(hv)=="HypervolumeList")
  {
    result = (sapply(hv@HVList, function(x) { apply(x@RandomPoints,2,mean) }))
    dimnames(result)[[2]] <- sapply(hv@HVList, function(x) {x@Name})
    return(t(result))
  }
  else
  {
    stop('Wrong class input.')
  }
}