hypervolume_centroid <- function(hv)
{
  centroid <- colMeans(hv@RandomUniformPointsThresholded, na.rm=TRUE)
  
  return(centroid)
}