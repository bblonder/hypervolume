hypervolume_sorensen_overlap <- function(hvlist)
{
  value <- 2 * hvlist@HVList$Intersection@Volume / (hvlist@HVList$HV1@Volume + hvlist@HVList$HV2@Volume)
  return(value)
}