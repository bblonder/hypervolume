evalfspherical <- function(data, radius, points, getid.nearestneighbor=FALSE, verbose=TRUE)
{
  result = rep(NA, nrow(points))
  
  points_numeric = t(as.matrix(points))
  nr = nrow(points)
  nc = ncol(points)
  
  tree <- kdtree_build(data,verbose=verbose)

  result <- kdtree_ball_query_multiple(tree, points_numeric, nr, nc, radius, verbose)

  if (getid.nearestneighbor==TRUE)
  {
    result_id <- kdtree_ball_query_id_multiple(tree, points_numeric, nr, nc, radius, verbose)
  }
  
  rm(tree)

  if (getid.nearestneighbor==TRUE)
  {
    return(list(result, result_id))
  }
  else
  {
    return(result)
  }
  
  
}