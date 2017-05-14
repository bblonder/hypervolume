evalfrectangular <- function(data, bandwidth, points, verbose=TRUE)
{
  result = rep(NA, nrow(points))
  
  tree <- kdtree_build(data,verbose=verbose)
  
  pointsfinalmax <- points + repmat(t(as.matrix(bandwidth,nrow=1)), nrow(points), 1)
  pointsfinalmin <- points - repmat(t(as.matrix(bandwidth,nrow=1)), nrow(points), 1)
  
  pointsmax_numeric = t(as.matrix(pointsfinalmax))
  pointsmin_numeric = t(as.matrix(pointsfinalmin))
  
  nr = nrow(pointsfinalmax)
  nc = ncol(pointsfinalmax)

  result <- kdtree_range_query_multiple(tree, pointsmin_numeric, pointsmax_numeric, nr, nc, verbose)

  rm(tree)
  
  return(result)
}