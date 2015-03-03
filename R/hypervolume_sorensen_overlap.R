hypervolume_sorensen_overlap <- function(hvlist)
{
  hvlist@HVList$Intersection@Volume / hvlist@HVList$Union@Volume
}