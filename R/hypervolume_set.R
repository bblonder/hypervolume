hypervolume_set <- function(hv1, hv2, num.points.max=NULL, verbose=TRUE, check.memory=TRUE, distance.factor=1.0)
{
  # determine dataset sizes and dimensionality
  np1 = nrow(hv1@RandomPoints)
  np2 = nrow(hv2@RandomPoints)
  
  hv1_point_density = np1 / hv1@Volume
  if(is.nan(hv1_point_density)) { hv1_point_density <- NA }
  hv2_point_density = np2 / hv2@Volume
  if(is.nan(hv2_point_density)) { hv2_point_density <- NA }
  
  dimhv1 = ncol(hv1@RandomPoints)
  dimhv2 = ncol(hv2@RandomPoints)
  
  if (dimhv1 != dimhv2)
  {
    stop('Dimensionality of hypervolumes is not the same.')
  }
  
  # define the consistent dimensionality
  dim = dimhv1
  
  if (is.null(num.points.max))
  {
    num.points.max = ceiling(10^(3+sqrt(hv1@Dimensionality)))
    if (verbose==TRUE)
    {
      cat(sprintf('Choosing num.points.max=%.0f (use a larger value for more accuracy.)\n',num.points.max))  
    }
  }
  
  # sample both hypervolumes down to the minimum point density
  mindensity = min(c(hv1_point_density, hv2_point_density, num.points.max / hv1@Volume , num.points.max / hv2@Volume), na.rm=T)
  if (verbose==TRUE)
  {
    cat(sprintf('Using minimum density of %f\n', mindensity))
  }
  
  numpointstokeep_hv1 = floor(mindensity * hv1@Volume)
  numpointstokeep_hv2 = floor(mindensity * hv2@Volume)
  
  hv_empty <- hv1
  hv_empty@Method = "Set operations"
  hv_empty@Data = matrix(NaN,nrow=1,ncol=dim)
  hv_empty@Dimensionality = dim
  hv_empty@Volume = 0
  hv_empty@PointDensity = mindensity
  hv_empty@Parameters = list()
  hv_empty@RandomPoints = matrix(NA,nrow=0,ncol=dim,dimnames=dimnames(hv1@RandomPoints))  
  hv_empty@ValueAtRandomPoints = 0
  
  
  # handle the cases when one of the hypervolumes is empty
  # hv1 empty, hv2 not
  if ((numpointstokeep_hv1 == 0 | is.null(numpointstokeep_hv1)) & !(numpointstokeep_hv2 == 0 | is.null(numpointstokeep_hv2)))
  {
    warning('hv1 has no random points and is empty.')
    
    result = new("HypervolumeList")
    hv_int = hv_empty
    hv_int@Name = sprintf("Intersection of (%s, %s)", hv1@Name, hv2@Name)
    
    hv_union = hv2
    hv_union@Name = sprintf("Union of (%s, %s)", hv1@Name, hv2@Name)
    
    hv_unique_1 = hv_empty
    hv_unique_1@Name = sprintf("Unique component of (%s) relative to (%s)", hv1@Name, hv2@Name)
    
    hv_unique_2 = hv2
    hv_unique_2@Name = sprintf("Unique component of (%s) relative to (%s)", hv2@Name, hv1@Name)
    
    result@HVList = list(
      HV1 = hv1,
      HV2 = hv2,
      Intersection = hv_int, 
      Union = hv_union, 
      Unique_1 = hv_unique_1, 
      Unique_2 = hv_unique_2
    )
    
    return(result)
  }
  # hv2 empty, hv1 not
  if (!(numpointstokeep_hv1 == 0 | is.null(numpointstokeep_hv1)) & (numpointstokeep_hv2 == 0 | is.null(numpointstokeep_hv2)))
  {
    warning('hv2 has no random points and is empty.')
    
    result = new("HypervolumeList")
    hv_int = hv_empty
    hv_int@Name = sprintf("Intersection of (%s, %s)", hv1@Name, hv2@Name)
    
    hv_union = hv1
    hv_union@Name = sprintf("Union of (%s, %s)", hv1@Name, hv2@Name)
    
    hv_unique_1 = hv1
    hv_unique_1@Name = sprintf("Unique component of (%s) relative to (%s)", hv1@Name, hv2@Name)
    
    hv_unique_2 = hv_empty
    hv_unique_2@Name = sprintf("Unique component of (%s) relative to (%s)", hv2@Name, hv1@Name)
    
    result@HVList = list(
      HV1 = hv1,
      HV2 = hv2,
      Intersection = hv_int, 
      Union = hv_union, 
      Unique_1 = hv_unique_1, 
      Unique_2 = hv_unique_2
    )
    
    return(result)
  }
  # hv1 and hv2 empty
  if ((numpointstokeep_hv1 == 0 | is.null(numpointstokeep_hv1)) & (numpointstokeep_hv2 == 0 | is.null(numpointstokeep_hv2)))
  {
    warning('hv1 and hv2 both have no random points and are empty.')
    
    result = new("HypervolumeList")
    hv_int = hv_empty
    hv_int@Name = sprintf("Intersection of (%s, %s)", hv1@Name, hv2@Name)
    
    hv_union = hv_empty
    hv_union@Name = sprintf("Union of (%s, %s)", hv1@Name, hv2@Name)
    
    hv_unique_1 = hv_empty
    hv_unique_1@Name = sprintf("Unique component of (%s) relative to (%s)", hv1@Name, hv2@Name)
    
    hv_unique_2 = hv_empty
    hv_unique_2@Name = sprintf("Unique component of (%s) relative to (%s)", hv2@Name, hv1@Name)
    
    result@HVList = list(
      HV1 = hv1,
      HV2 = hv2,
      Intersection = hv_int, 
      Union = hv_union, 
      Unique_1 = hv_unique_1, 
      Unique_2 = hv_unique_2
    )
    
    return(result)    
  }

  
  # if the algorithm passed returning early from any of the degenerate cases...
  
  if (verbose==TRUE | check.memory==TRUE)
  {
    cat(sprintf('Retaining %d points in hv1 and %d points in hv2.\n', numpointstokeep_hv1, numpointstokeep_hv2))
    if (check.memory == TRUE)
    {
      ans <- message(sprintf('This will require %.0f pairwise comparisons. Re-run function with check.memory=FALSE if acceptable; otherwise use a smaller value of reduction_factor.\n', numpointstokeep_hv1*numpointstokeep_hv2))
      return(NULL)
    }
  }

  hv1_points_ss = hv1@RandomPoints[sample(1:np1,size=numpointstokeep_hv1),,drop=FALSE]
  hv2_points_ss = hv2@RandomPoints[sample(1:np2,size=numpointstokeep_hv2),,drop=FALSE]
  
  point_density = nrow(hv1_points_ss) / hv1@Volume
  
  # calculate characteristic distances
  cutoff_dist = point_density^(-1/dim) * distance.factor
  
  
  # figure out which points are 'close enough' to other points
  # (within a n-ball of the critical distance)
  if (verbose==TRUE)
  {
    cat('Beginning ball queries... \n')
  }

  p2_in_1_all = evalfspherical(data=hv1_points_ss, radius=cutoff_dist, points=hv2_points_ss)
  p1_in_2_all = evalfspherical(data=hv2_points_ss, radius=cutoff_dist, points=hv1_points_ss)

  if (verbose==TRUE)
  {
    cat('Finished ball queries. \n')
  }
  
  # subset to retain only those 'close enough' points
  p2_in_1 = as.data.frame(hv2_points_ss)[p2_in_1_all > 0,,drop=FALSE]
  p1_in_2 = as.data.frame(hv1_points_ss)[p1_in_2_all > 0,,drop=FALSE]
  
  
  # the final volume is proportional to the fraction 
  # of points in hv1 in hv2, and vice versa
  v1 = nrow(p1_in_2) / nrow(hv1_points_ss) * hv1@Volume
  v2 = nrow(p2_in_1) / nrow(hv2_points_ss) * hv2@Volume
  
  
  # take the lower estimate as a conservative estimate
  final_volume_intersection = min(c(v1,v2))
  
  # create the intersection point cloud by merging both sets of sampled points.
  final_points_intersection = unique(rbind(p1_in_2, p2_in_1))
  final_density_intersection = nrow(final_points_intersection) / final_volume_intersection
  
  
  # now find the union point cloud
  p1_not_in_2 = hv1_points_ss[p1_in_2_all == 0,,drop=FALSE] # the points only in the first hypervolume
  p2_not_in_1 = hv2_points_ss[p2_in_1_all == 0,,drop=FALSE] # the points only in the second hypervolume
  
  num_points_to_sample_in_intersection = floor(point_density * final_volume_intersection) # choose the right number of points to keep the point density constant
  
  p_in_1_and_2 = final_points_intersection[sample(1:nrow(final_points_intersection), size=num_points_to_sample_in_intersection),,drop=FALSE] # randomly sample the intersection to grab
  
  final_volume_union = hv1@Volume + hv2@Volume - final_volume_intersection # union is sum minus intersection 
  final_points_union = unique(rbind(p1_not_in_2, p2_not_in_1, p_in_1_and_2)) 
  final_density_union = nrow(final_points_union) / final_volume_union
  
  # calculate the unique components for hv1 
  # (these will occasionally be too permissive and show some outlying points)
  final_volume_unique_hv1 = hv1@Volume - final_volume_intersection
  final_points_in_unique_1 = unique(p1_not_in_2)
  final_density_unique_1 = nrow(final_points_in_unique_1) / final_volume_unique_hv1
  
  # calculate the unique components for hv2
  final_volume_unique_hv2 = hv2@Volume - final_volume_intersection
  final_points_in_unique_2 = unique(p2_not_in_1)
  final_density_unique_2 = nrow(final_points_in_unique_2) / final_volume_unique_hv2
  
  # get column names
  cn <- dimnames(hv1@RandomPoints)[[2]]
  dn <- list(NULL, cn)
  
  # prepare final hypervolumes
  result_intersection = new("Hypervolume")
  result_intersection@Method = "Set operations"
  result_intersection@Name = sprintf("Intersection of (%s, %s)", hv1@Name, hv2@Name)
  result_intersection@Data = matrix(NaN,nrow=1,ncol=dim)
  result_intersection@Dimensionality = dim
  result_intersection@Volume = final_volume_intersection
  result_intersection@PointDensity = final_density_intersection
  result_intersection@Parameters = list()
  result_intersection@RandomPoints = matrix(as.matrix(as.data.frame(final_points_intersection)),ncol=dim)  
  dimnames(result_intersection@RandomPoints) = dn
  result_intersection@ValueAtRandomPoints = rep(1,nrow(result_intersection@RandomPoints))
  
  result_union = new("Hypervolume")
  result_union@Name = sprintf("Union of (%s, %s)", hv1@Name, hv2@Name)
  result_union@Method = "Set operations"
  result_union@Data = matrix(NaN,nrow=1,ncol=dim)
  result_union@Dimensionality = dim
  result_union@Volume = final_volume_union
  result_union@PointDensity = final_density_union
  result_union@Parameters = list()
  result_union@RandomPoints = matrix(as.matrix(as.data.frame(final_points_union)),ncol=dim)
  dimnames(result_union@RandomPoints) = dn
  result_union@ValueAtRandomPoints = rep(1,nrow(result_union@RandomPoints))
  
  result_unique_hv1 = new("Hypervolume")
  result_unique_hv1@Name = sprintf("Unique component of (%s) relative to (%s)", hv1@Name, hv2@Name)
  result_unique_hv1@Data = matrix(NaN,nrow=1,ncol=dim)
  result_unique_hv1@Method = "Set operations"
  result_unique_hv1@Dimensionality = dim
  result_unique_hv1@Volume = final_volume_unique_hv1
  result_unique_hv1@PointDensity = final_density_unique_1
  result_unique_hv1@Parameters = list()
  result_unique_hv1@RandomPoints = matrix(as.matrix(as.data.frame(final_points_in_unique_1)),ncol=dim)
  dimnames(result_unique_hv1@RandomPoints) = dn
  result_unique_hv1@ValueAtRandomPoints = rep(1,nrow(result_unique_hv1@RandomPoints))
  
  result_unique_hv2 = new("Hypervolume")
  result_unique_hv2@Name = sprintf("Unique component of (%s) relative to (%s)", hv2@Name, hv1@Name)
  result_unique_hv2@Data = matrix(NaN,nrow=1,ncol=dim)
  result_unique_hv2@Method = "Set operations"
  result_unique_hv2@Dimensionality = dim
  result_unique_hv2@Volume = final_volume_unique_hv2
  result_unique_hv2@PointDensity = final_density_unique_2
  result_unique_hv2@Parameters = list()
  result_unique_hv2@RandomPoints = matrix(as.matrix(as.data.frame(final_points_in_unique_2)),ncol=dim)
  dimnames(result_unique_hv2@RandomPoints) = dn
  result_unique_hv2@ValueAtRandomPoints = rep(1,nrow(result_unique_hv2@RandomPoints))
  
  # assemble final results into a list
  result = new("HypervolumeList")
  result@HVList = list(
    HV1 = hv1,
    HV2 = hv2,
    Intersection = result_intersection, 
    Union = result_union, 
    Unique_1 = result_unique_hv1, 
    Unique_2 = result_unique_hv2
  )
  
  return(result)
}