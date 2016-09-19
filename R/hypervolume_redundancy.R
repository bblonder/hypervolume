hypervolume_redundancy <- function(...)
{
  p <- hypervolume_estimate_probability(...)
  return(p^2)
}