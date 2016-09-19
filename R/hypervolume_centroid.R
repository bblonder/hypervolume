hypervolume_centroid <- function(hv)
{
  centroid <- colMeans(hv@RandomUniformPointsThresholded, na.rm=T)
  
  return(centroid)
}