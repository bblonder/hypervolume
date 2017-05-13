evalfrectangular_multiple <- function(data, bandwidth_list, points, verbose=TRUE)
{
  if (verbose==TRUE) {cat('\nBuilding tree...')}
  tree = kdtree_build(data)
  if (verbose==TRUE) {cat(' done.\n')}
  result = rep(NA, nrow(points))
  
  counts <- sapply(bandwidth_list, function(bandwidth) {
    pointsfinalmax <- points + repmat(t(as.matrix(bandwidth,nrow=1)), nrow(points), 1)
    pointsfinalmin <- points - repmat(t(as.matrix(bandwidth,nrow=1)), nrow(points), 1)
    
    pointsmax_numeric = t(as.matrix(pointsfinalmax))
    pointsmin_numeric = t(as.matrix(pointsfinalmin))
    
    nr = nrow(pointsfinalmax)
    nc = ncol(pointsfinalmax)
    
    if (verbose==TRUE) {cat('Querying tree...')}
    result <- kdtree_range_query_multiple(tree, pointsmin_numeric, pointsmax_numeric, nr, nc, verbose)
    if (verbose==TRUE) {cat(' done.\n')}
    

    return(result)
  
  })
  
  rm(tree); gc(reset=TRUE); # make sure the memory is released
  
  return(counts)
}

evalfrectangular <- function(data, bandwidth, points, verbose=TRUE)
{
  if (verbose==TRUE) {cat('\nBuilding tree...')}
  tree = kdtree_build(data)
  if (verbose==TRUE) {cat(' done.\n')}
  result = rep(NA, nrow(points))
  
  pointsfinalmax <- points + repmat(t(as.matrix(bandwidth,nrow=1)), nrow(points), 1)
  pointsfinalmin <- points - repmat(t(as.matrix(bandwidth,nrow=1)), nrow(points), 1)
  
  pointsmax_numeric = t(as.matrix(pointsfinalmax))
  pointsmin_numeric = t(as.matrix(pointsfinalmin))
  
  nr = nrow(pointsfinalmax)
  nc = ncol(pointsfinalmax)
  
  if (verbose==TRUE) {cat('Querying tree...')}
  result <- kdtree_range_query_multiple(tree, pointsmin_numeric, pointsmax_numeric, nr, nc, verbose)
  if (verbose==TRUE) {cat(' done.\n')}
  
  rm(tree); gc(reset=TRUE); # make sure the memory is released
  
  return(result)
}