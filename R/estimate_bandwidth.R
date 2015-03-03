estimate_bandwidth <- function(data)
{
  n = nrow(data)
  stdev = apply(data,2,sd)
  
  normal_silverman_bw = 1.06 * stdev * n ^ (-1/5)
  
  rect_bw = 2*normal_silverman_bw
  
  return(rect_bw)
}