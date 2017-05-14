
hypervolume_overlap_statistics <- function(hvlist)
{
  sorensen <- 2 * hvlist@HVList$Intersection@Volume / (hvlist@HVList$HV1@Volume + hvlist@HVList$HV2@Volume)

  frac_unique_1 <- hvlist@HVList$Unique_1@Volume / (hvlist@HVList$HV1@Volume)
  frac_unique_2 <- hvlist@HVList$Unique_2@Volume / (hvlist@HVList$HV2@Volume)
  
  jaccard <- hvlist@HVList$Intersection@Volume / hvlist@HVList$Union@Volume
  
  return(c(jaccard=jaccard,sorensen=sorensen,frac_unique_1=frac_unique_1, frac_unique_2=frac_unique_2))
}

hypervolume_sorensen_overlap <- function(hvlist)
{
  warning("This function is deprecated; use hypervolume_overlap_statistics instead.");
  warning("The definition of overlap has changed since version 1.2 - a prefactor of 2 instead of 1.")
  
  return(hypervolume_overlap_statistics(hvlist)$sorensen)
}
