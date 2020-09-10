to_hv_list <- function(path) {
  # Returns hypervolume list object
  hvs = new("HypervolumeList")
  
  # adds every hypervolume object in the directory to list
  hvs@HVList = foreach(file = list.files(path), .combine = c) %do% {
    readRDS(file.path(path, file))
  }
  return(hvs)
}