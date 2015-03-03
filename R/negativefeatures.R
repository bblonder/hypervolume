randsphere <- function(m,n,r)
{
  X = matrix(data=rnorm(m*n),nrow=m,ncol=n)
  s2 = rowSums(X * X)

  X = X * matrix(rep(r*(pgamma(s2/2,n/2)^(1/n))/sqrt(s2),n),ncol=n)
    
  return(X)
}

nball_volume <- function(n, r)
{
  return(pi^(n/2) * r^n / gamma(n/2+1))
}

expectation_maximal <- function(input, ...)
{
  if(class(input)=="Hypervolume")
  {
    return(input)
  }
  else
  {
    return(hypervolume(input, name=sprintf("Maximal expectation for %s", deparse(substitute(input))[1]), ...))
  }
}

expectation_box <- function(input, npoints=NULL, userandom=FALSE)
{
  if (class(input)=="Hypervolume")
  {
    if (userandom==TRUE)
    {
      data <- input@RandomUniformPointsThresholded
    }
    else 
    {
      data <- input@Data
    }
  }
  else
  {
    data <- input
    

  }
  
  minv <- apply(data,2,min)
  maxv <- apply(data,2,max)
  
  volume = prod(maxv - minv)
  
  if (is.null(npoints))
  {
    if (class(input)=="Hypervolume")
    {
      npoints = input@PointDensity * volume
    }
    else 
    {
      npoints <- 10 * 10^ncol(input)
    }
    cat(sprintf('Choosing npoints=%.0f (use a larger value for more accuracy.)\n',npoints))
  }
  
  density = npoints / volume
  
  result <- matrix(NA, nrow=npoints, ncol=length(minv), dimnames=list(NULL, dimnames(data)[[2]]))
  
  for (i in 1:length(minv)) # across each dimension
  {
    result[,i] <- runif(npoints, minv[i], maxv[i])
  }

  hv_box <- new("Hypervolume",
                Name=sprintf("Box expectation for %s", ifelse(class(input)=="Hypervolume", input@Name, deparse(substitute(data))[1])),
                Data = as.matrix(data),
                RandomUniformPointsThresholded=result, 
                Dimensionality=ncol(data), 
                Volume=volume,
                PointDensity=density, 
                DisjunctFactor=NaN,
                Bandwidth= NaN, 
                RepsPerPoint=floor(npoints / nrow(data)), 
                QuantileThresholdDesired=0, 
                QuantileThresholdObtained=0,
                ProbabilityDensityAtRandomUniformPoints = rep(1, npoints)
  )  
  
  return(hv_box)
}


expectation_ball <- function(input, npoints=NULL, userandom=FALSE)
{
  if (class(input)=="Hypervolume")
  {
    if (userandom==TRUE)
    {
      data <- input@RandomUniformPointsThresholded
    }
    else 
    {
      data <- input@Data
    }
  }
  else
  {
    data <- input
    
    
  }
  
  center <- apply(data,2,mean,na.rm=T)
  data_centered <- sweep(data, 2, center, "-")
  radii = sqrt(rowSums(data_centered^2))
  maxradius = max(radii)
  
  volume = nball_volume(n=ncol(data),r=maxradius)
  
  if (is.null(npoints))
  {
    if (class(input)=="Hypervolume")
    {
      npoints = input@PointDensity * volume
    }
    else 
    {
      npoints <- 10 * 10^ncol(data)
    }
    cat(sprintf('Choosing npoints=%.0f (use a larger value for more accuracy.)\n',npoints))
  }
  
  density = npoints / volume
  
  result <- randsphere(m=npoints, n=ncol(data), r=maxradius)
  result <- sweep(result, 2, center,"+")
  dimnames(result) <- list(NULL, dimnames(data)[[2]])
  
  hv_ball <- new("Hypervolume",
                Name=sprintf("Ball expectation for %s", ifelse(class(input)=="Hypervolume", input@Name, deparse(substitute(data))[1])),
                Data = as.matrix(data),
                RandomUniformPointsThresholded=result, 
                Dimensionality=ncol(data), 
                Volume=volume,
                PointDensity=density, 
                DisjunctFactor=NaN,
                Bandwidth= NaN, 
                RepsPerPoint=floor(npoints / nrow(data)), 
                QuantileThresholdDesired=0, 
                QuantileThresholdObtained=0,
                ProbabilityDensityAtRandomUniformPoints = rep(1, npoints)
  )  
  
  return(hv_ball)
}



expectation_convex <- function(input, npoints_inhull=NULL, npoints_onhull=NULL, check_memory=TRUE, userandom=FALSE, method="rejection", burnin=NULL, delta=NULL)
{  
  if (class(input)=="Hypervolume")
  {
    if (userandom==TRUE)
    {
      data <- input@RandomUniformPointsThresholded
    }
    else 
    {
      data <- input@Data
    }
  }
  else
  {
    data <- input

  }  

  if (is.null(npoints_onhull))
  {
    npoints_onhull <- min(nrow(data),floor(10^sqrt(ncol(data))))
    cat(sprintf('Choosing npoints_onhull=%.0f (use a larger value for more accuracy.)\n',npoints_onhull))
  }

  numconvexhullpoints <- min(nrow(data),npoints_onhull)

  data_reduced <- data[sample(nrow(data), numconvexhullpoints, prob=rowSums(scale(data, center=TRUE, scale=FALSE)^2)),]
  
  if (ncol(data) > 5)
  {
    warning(sprintf("Algorithm may crash or have very high runtime on high dimensional data (n=%d",ncol(data)))
  }
  
  # FIND THE CONVEX HULL of the reduced data	
  convexhull <- geometry::convhulln(data_reduced,options="FA")
  hull_matrix <- convexhull$hull #extract coordinates of vertices
  hull_volume <- convexhull$vol # extract volume
  
  cat(sprintf("Convex hull calculated with %.0f simplices.\nCalculation of inequality constraints will require allocation of %.0f double-precision numbers.\n", nrow(hull_matrix), nrow(hull_matrix)^2))
  
  if (check_memory)
  {
    stop('Set check_memory=F to continue.\n')
  }	
  
  if (class(input)=="Hypervolume" && is.null(npoints_inhull))
  {
    print(input@PointDensity)
    print(hull_volume)
    npoints_inhull = input@PointDensity * hull_volume
    cat(sprintf('Matching density: choosing npoints_inhull=%.0f (use a larger value for more accuracy.)\n',npoints_inhull))
  }  
  
  if (class(input)!="Hypervolume" && is.null(npoints_inhull))
  {
    npoints_inhull = floor(10*10^ncol(data))
    cat(sprintf('Choosing npoints_inhull=%.0f (use a larger value for more accuracy.)\n',npoints_inhull))
  }  
  
  
  if (method=="rejection")
  {
    # REJECTION SAMPLING FROM BOUNDING BOX TO FILL IN CONVEX POLYTOPE APPROACH  
    ntrials <- 1000 # number of candidate points to propose at a time (algorithm always works regardless of this choice)
    done = FALSE
    inpoints <- NULL
    niter <- 0
    
    # generate a bounding box, and test if each point is in the convex hull or not
    while(!done)
    {
      # try a set of test points
      testpoints <- expectation_box(data_reduced, npoints=ntrials)

      chullin <- (inhull(testpts= testpoints@RandomUniformPointsThresholded, calpts=data_reduced, hull=hull_matrix)==1)
      # figure out which are 'in'
      inpoints_temp <- testpoints@RandomUniformPointsThresholded[chullin,]
      
      niter <- niter + 1
      inpoints <- rbind(inpoints, inpoints_temp)
      
      done = nrow(inpoints) >= npoints_inhull
      
      cat(sprintf('Rejection sampling: iteration %.0f - %.0f / %.0f points accepted\n', niter, nrow(inpoints), npoints_inhull))
    }
  }
  else if (method=="metropolis")
  {
    isin <- function(candidate)
    {
      inhull(testpts=matrix(candidate, nrow=1), calpts=data_reduced, hull=hull_matrix)
    }
    
    if (is.null(burnin))
    {
      burnin <- 10^ncol(data_reduced)
    }
    if (is.null(delta))
    {
      delta <- diff(range(data_reduced))/100
    }
    cat(sprintf('Running Metropolis algorithm with burnin=%d delta=%f.\nAlgorithm NOT RECOMMENDED unless dimensionality is high (n > 6). Remember that sample is not guaranteed to be uniformly random unless burnin >> 1.\n', burnin, delta))
    
    inpoints <- rwmetro(target=isin, N=npoints_inhull, x=colMeans(data_reduced, na.rm=TRUE), burnin=burnin, delta=delta)

  }
  else
  {
    stop('Invalid convex hull sampling method specified.')
  }
  
  dimnames(inpoints) <- list(NULL,dimnames(data)[[2]])
  
  # MAKE a new hypervolume with the convex hull shape and volume
  hv_chull <- new("Hypervolume",
                  Name=sprintf("Convex expectation for %s", ifelse(class(input)=="Hypervolume", input@Name, deparse(substitute(data))[1])),
                  Data=as.matrix(data),
                  RandomUniformPointsThresholded=inpoints, 
                  Dimensionality=ncol(inpoints), 
                  Volume=hull_volume,
                  DisjunctFactor=NaN,
                  PointDensity = nrow(inpoints) / hull_volume,
                  #PointDensity=niter * testpoints@PointDensity, 
                  Bandwidth= NaN, 
                  RepsPerPoint=NaN, 
                  QuantileThresholdDesired=0, 
                  QuantileThresholdObtained=0,
                  ProbabilityDensityAtRandomUniformPoints = rep(1, nrow(inpoints))
  )
  
  return(hv_chull)
}

negative_features <- function(hv_obs, hv_exp, set_npoints_max=NULL, set_check_memory=TRUE)
{		
  # initialize result
  finalresult <- NULL
  
  if (is.null(set_npoints_max))
  {
    set_npoints_max = floor(100*10^(hv_obs@Dimensionality))
    cat(sprintf('Choosing set_npoints_max=%.0f (choose a larger value for more accuracy.)\n',set_npoints_max))    
  }
  
  # make sure we are using a hypervolume with data (i.e. not the output of set operations)
  if(all(is.nan(hv_obs@Data)))
  {
    stop('Hypervolume must be associated with datapoints')
  }
  
  if (hv_obs@Dimensionality != hv_exp@Dimensionality)
  {
    stop('Observed and expected hypervolumes must have same dimensionality.')
  }
  
  # FIND THE DIFFERENCE between the convex hull shape and the real hypervolume
  cat("Beginning set operations (resampling to minimum density)...")
  hvs_overlap <- hypervolume_set(hv_obs, hv_exp, check_memory=set_check_memory, npoints_max=set_npoints_max)
  if (set_check_memory)
  {
    stop('Set set_check_memory=F to continue.\n')
  }
  cat("Finished set operations.\n")
  
  # find the distance between all points in the difference
  randompoints <- hvs_overlap@HVList$Unique_2@RandomUniformPointsThresholded
  if (is.null(randompoints) || nrow(randompoints) == 0)
  {		
    cat('No negative features found.\n');
  }
  else
  {  	
    cat(sprintf("Retaining %.0f random points in set difference.\n", nrow(randompoints)))

    criticaldistance <- hvs_overlap@HVList$Unique_2@PointDensity ^(-1/hvs_overlap@HVList$Unique_2@Dimensionality)    
    
    # find points with minimum neighbor distance less than threshold
    distances <- as.matrix(dist(randompoints, method="euclidean"))
    diag(distances) <- NA
    isin <- (apply(distances, 1, min, na.rm=T) < criticaldistance)
    
    cat(sprintf("Removing %d stray points...\n", length(which(isin==0))))
    
    randompoints_trimmed <- randompoints[isin,]
    
    thishv <- hvs_overlap@HVList$Unique_2 # copy base information
    thishv@RandomUniformPointsThresholded <- randompoints_trimmed
    thishv@Volume <- hvs_overlap@HVList$Unique_2@Volume * nrow(randompoints_trimmed) / nrow(randompoints)
    thishv@Name <- sprintf("Negative features of %s relative to %s", hv_obs@Name, hv_exp@Name)
    
    finalresult <- thishv

    # return the final hypervolumelist
    cat('Returning all negative features.\n');
    
    if (is.null(finalresult))
    {
      cat('No negative features retained. This may indicate one or more of the following scenarios\n\t- npoints_min_cluster is too large.\n\t- no negative features\n\t- low PointDensity in input hypervolumes\n\t- reduction_factor is too low.\n\nFunction will return NULL.\n')
    }
        
  }
  
  return(finalresult)
}