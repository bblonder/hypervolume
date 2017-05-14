
get_volume.Hypervolume <- function(object)
{
  result = object@Volume
  names(result) = object@Name
  
  return(result)
}

get_volume.HypervolumeList <- function(object)
{
  result = sapply(object@HVList, get_volume.Hypervolume)
  names(result) = sapply(object@HVList, function(x) {x@Name})
  return(result)
} 


setGeneric("get_volume", function(object) {})
setMethod("get_volume","Hypervolume", function(object) {get_volume.Hypervolume(object)})
setMethod("get_volume","HypervolumeList", function(object) {get_volume.HypervolumeList(object)})
