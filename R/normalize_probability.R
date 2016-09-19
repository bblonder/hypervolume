normalize_probability <- function(probability_raw, point_density)
{
  k <- point_density / sum(probability_raw)
  
  return(probability_raw * k)
}