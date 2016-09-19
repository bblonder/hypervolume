evalfspherical <- function(data, radius, points, verbose=T, getid.nearestneighbor=F)
{
  if (verbose==TRUE) {cat('Building tree...')}
  tree = kdtree_build(data)
  if (verbose==TRUE) {cat(' done.\n')}
  result = rep(NA, nrow(points))
  
  points_numeric = t(as.matrix(points))
  nr = nrow(points)
  nc = ncol(points)
  
  if (verbose==TRUE) {cat('Querying tree...\n')}

  result <- kdtree_ball_query_multiple(tree, points_numeric, nr, nc, radius, verbose)

  if (getid.nearestneighbor==TRUE)
  {
    result_id <- kdtree_ball_query_id_multiple(tree, points_numeric, nr, nc, radius, verbose)
  }
  
  if (verbose==TRUE) {cat(' done.\n')}
  
  rm(tree); gc(reset=T); # make sure the memory is released


  if (getid.nearestneighbor==TRUE)
  {
    return(list(result, result_id))
  }
  else
  {
    return(result)
  }
}