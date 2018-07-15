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

expectation_box <- function(input, point.density=NULL, num.samples=NULL, use.random=FALSE)
{
  if (class(input)=="Hypervolume")
  {
    if (use.random==TRUE)
    {
      data <- input@RandomPoints
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
  
  if        (is.null(point.density) & !is.null(num.samples)) # numpoints specified, use it
  {
    npoints = num.samples
    point.density = npoints / volume
  }
  else if   (is.null(point.density) & is.null(num.samples)) # no input, use default # of points
  {
    npoints = ceiling(10^(3+sqrt(ncol(data))))
    point.density = npoints / volume
  }
  else if   (!is.null(point.density) & is.null(num.samples)) # point density specified, use it
  {
    npoints = ceiling(volume * point.density)
    point.density = point.density
  }
  else if   (!is.null(point.density) & !is.null(num.samples)) # error
  {
    stop('Cannot specify both point.density and num.samples')
  }
  
  result <- matrix(NA, nrow=npoints, ncol=length(minv), dimnames=list(NULL, dimnames(data)[[2]]))
  
  for (i in 1:length(minv)) # across each dimension
  {
    result[,i] <- runif(npoints, minv[i], maxv[i])
  }
  
  hv_box <- new("Hypervolume",
                Name=sprintf("Box expectation for %s", ifelse(class(input)=="Hypervolume", input@Name, deparse(substitute(data))[1])),
                Method="Box expectation",
                Data = as.matrix(data),
                RandomPoints=result, 
                Dimensionality=ncol(data), 
                Volume=volume,
                PointDensity=point.density, 
                Parameters= list(),
                ValueAtRandomPoints = rep(1, npoints)
  )  
  
  return(hv_box)
}


expectation_ball <- function(input, point.density=NULL, num.samples=NULL, use.random=FALSE)
{
  if (class(input)=="Hypervolume")
  {
    if (use.random==TRUE)
    {
      data <- input@RandomPoints
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

  center <- apply(data,2,mean,na.rm=TRUE)
  data_centered <- sweep(data, 2, center, "-")
  radii = sqrt(rowSums(data_centered^2))
  maxradius = max(radii)
  
  volume = nball_volume(n=ncol(data),r=maxradius)  
  
  if        (is.null(point.density) & !is.null(num.samples)) # numpoints specified, use it
  {
    npoints = num.samples
    point.density = npoints / volume
  }
  else if   (is.null(point.density) & is.null(num.samples)) # no input, use default # of points
  {
    npoints = ceiling(10^(3+sqrt(ncol(data))))
    point.density = npoints / volume
  }
  else if   (!is.null(point.density) & is.null(num.samples)) # point density specified, use it
  {
    npoints = ceiling(volume * point.density)
    point.density = point.density
  }
  else if   (!is.null(point.density) & !is.null(num.samples)) # error
  {
    stop('Cannot specify both point.density and num.samples')
  }
  
  result <- randsphere(m=npoints, n=ncol(data), r=maxradius)
  result <- sweep(result, 2, center,"+")
  dimnames(result) <- list(NULL, dimnames(data)[[2]])
  
  hv_ball <- new("Hypervolume",
                 Name=sprintf("Ball expectation for %s", ifelse(class(input)=="Hypervolume", input@Name, deparse(substitute(data))[1])),
                 Method="Ball expectation",
                 Data = as.matrix(data),
                 RandomPoints=result, 
                 Dimensionality=ncol(data), 
                 Volume=volume,
                 PointDensity=point.density, 
                 Parameters= list(), 
                 ValueAtRandomPoints = rep(1, npoints)
  )  
  
  return(hv_ball)
}

# from Keith Jewell
norm_compiled <- cmpfun(function(x) {sqrt(sum(x^2))})

inhull <- function(testpts, calpts, hull=geometry::convhulln(calpts), tol=mean(mean(abs(calpts)))*sqrt(.Machine$double.eps)) {
  #++++++++++++++++++++
  # R implementation of the Matlab code by John D'Errico 04 Mar 2006 (Updated 30 Oct 2006)
  # downloaded from http://www.mathworks.com/matlabcentral/fileexchange/10226-inhull
  # with some modifications and simplifications
  #
  # Efficient test for points inside a convex hull in n dimensions
  #
  #% arguments: (input)
  #%  testpts - nxp array to test, n data points, in p dimensions
  #%       If you have many points to test, it is most efficient to
  #%       call this function once with the entire set.
  #%
  #%  calpts - mxp array of vertices of the convex hull, as used by
  #%       convhulln.
  #%
  #%  hull - (OPTIONAL) tessellation (or triangulation) generated by convhulln
  #%       If hull is left empty or not supplied, then it will be
  #%       generated.
  #%
  #%  tol - (OPTIONAL) tolerance on the tests for inclusion in the
  #%       convex hull. You can think of tol as the distance a point
  #%       may possibly lie outside the hull, and still be perceived
  #%       as on the surface of the hull. Because of numerical slop
  #%       nothing can ever be done exactly here. I might guess a
  #%       semi-intelligent value of tol to be
  #%
  #%         tol = 1.e-13*mean(abs(calpts(:)))
  #%
  #%       In higher dimensions, the numerical issues of floating
  #%       point arithmetic will probably suggest a larger value
  #%       of tol.
  #%
  # In this R implementation default tol=mean(mean(abs(calpts)))*sqrt(.Machine$double.eps)
  #       DEFAULT: tol = 1e-6
  #
  # VALUE: Matlab returns a vector of TRUE (inside/on) or FALSE (outside)
  #       This R implementation returns an integer vector of length n
  #       1 = inside hull
  #      -1 = inside hull
  #       0 = on hull (to precision indicated by tol)
  #--------------------------------------------------------
  # ensure arguments are matrices (not data frames) and get sizes
  calpts <- as.matrix(calpts)
  testpts <- as.matrix(testpts)
  p <- dim(calpts)[2]   # columns in calpts
  cx <- dim(testpts)[1]  # rows in testpts
  nt <- dim(hull)[1]    # number of simplexes in hull
  # find normal vectors to each simplex
  nrmls <- matrix(NA, nt, p)         # predefine each nrml as NA, degenerate
  degenflag <- matrix(TRUE, nt, 1)
  for (i in  1:nt) {
    nullsp <- t(MASS::Null(t(calpts[hull[i,-1],] - matrix(calpts[hull[i,1],],p-1,p, byrow=TRUE))))
    if (dim(nullsp)[1] == 1) { nrmls[i,] <- nullsp
    degenflag[i] <- FALSE}}
  # Warn of degenerate faces, and remove corresponding normals
  if(length(degenflag[degenflag]) > 0) 
    warning(length(degenflag[degenflag])," degenerate faces in convex hull")
  nrmls <- nrmls[!degenflag,]
  nt <- dim(nrmls)[1]
  # find center point in hull, and any (1st) point in the plane of each simplex
  center = apply(calpts, 2, mean)
  a <- calpts[hull[!degenflag,1],]
  # scale normal vectors to unit length and ensure pointing inwards
  nrmls <- nrmls/matrix(apply(nrmls, 1, norm_compiled), nt, p)
  dp <- sign(apply((matrix(center, nt, p, byrow=TRUE)-a) * nrmls, 1, sum))
  nrmls <- nrmls*matrix(dp, nt, p)
  # if  min across all faces of dot((x - a),nrml) is
  #      +ve then x is inside hull
  #      0   then x is on hull
  #      -ve then x is outside hull
  # Instead of dot((x - a),nrml)  use dot(x,nrml) - dot(a, nrml)
  tnrmls <- t(nrmls)
  aN <- diag(a %*% tnrmls)
  val <- apply(testpts %*% tnrmls - matrix(aN, cx, nt, byrow=TRUE), 1,min)
  # code  values inside 'tol' to zero, return sign as integer
  val[abs(val) < tol] <- 0
  as.integer(sign(val))
}

inhull_compiled <- cmpfun(inhull)



# modified by Gert van Valkenhoef - reference is
#Tervonen, T., van Valkenhoef, G., Basturk, N., & Postmus, D. (2012). Hit-And-Run enables efficient weight generation for simulation-based multiple criteria decision analysis. European Journal of Operational Research, 224(3), 552-559.
# http://dx.doi.org/10.1016/j.ejor.2012.08.026
hullToConstr <- function(calpts, hull) {
  calpts <- as.matrix(calpts)
  
  p <- ncol(calpts) # columns in calpts
  nt <- nrow(hull) # number of simplexes in hull
  
  nrmls <- matrix(NA, nt, p) # predefine each nrml as NA, degenerate
  degenflag <- rep(TRUE, nt)
  for (i in  1:nt) {
    nullsp <- t(Null(t(calpts[hull[i,-1],] - matrix(calpts[hull[i,1],],p-1,p, byrow=TRUE))))
    if (nrow(nullsp) == 1) {
      nrmls[i,] <- nullsp
      degenflag[i] <- FALSE
    }
  }
  
  # Warn of degenerate faces, and remove corresponding normals
  if(sum(degenflag) > 0) {
    warning(length(degenflag[degenflag])," degenerate faces in convex hull")
  }
  nrmls <- nrmls[!degenflag,]
  nt <- nrow(nrmls)
  
  # find center point in hull, and any (1st) point in the plane of each simplex
  center = apply(calpts, 2, mean)
  a <- calpts[hull[!degenflag,1],]
  
  # scale normal vectors to unit length and ensure pointing inwards
  nrmls <- nrmls/matrix(apply(nrmls, 1, function(x) sqrt(sum(x^2))), nt, p)
  dp <- sign(apply((matrix(center, nt, p, byrow=TRUE)-a) * nrmls, 1, sum))
  nrmls <- nrmls*matrix(dp, nt, p)
  
  nrmls <- -nrmls
  b <- diag(a %*% t(nrmls))
  list(constr=nrmls, rhs=b, dir=rep("<=", times=sum(!degenflag)))
}






expectation_convex <- function(input, point.density=NULL, num.samples=NULL, num.points.on.hull=NULL, check.memory=TRUE, verbose=TRUE, use.random=FALSE, method="hitandrun", chunksize=1e3)
{
  if (class(input)=="Hypervolume")
  {
    if (use.random==TRUE)
    {
      data <- input@RandomPoints
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
  
  if (is.null(num.points.on.hull))
  {
    num.points.on.hull <- nrow(data)
  }
  
  # convert to matrix
  data = as.matrix(data)
  
  # sample data down if needed, weighting further points more
  data_reduced <- data[sample(nrow(data), num.points.on.hull, replace=FALSE, prob=rowSums(scale(data, center=TRUE, scale=FALSE)^2)^4),]
  
  # FIND THE CONVEX HULL of the reduced data  
  convexhull <- geometry::convhulln(data_reduced,options="FA")
  hull_matrix <- convexhull$hull #extract coordinates of vertices
  hull_volume <- convexhull$vol # extract volume
  
 
  if        (is.null(point.density) & !is.null(num.samples)) # numpoints specified, use it
  {
    np = num.samples
    point.density = np / hull_volume
  }
  else if   (is.null(point.density) & is.null(num.samples)) # no input, use default # of points
  {
    np = ceiling(10^(3+sqrt(ncol(data))))
    point.density = np / hull_volume
  }
  else if   (!is.null(point.density) & is.null(num.samples)) # point density specified, use it
  {
    np = ceiling(hull_volume * point.density)
    point.density = point.density
  }
  else if   (!is.null(point.density) & !is.null(num.samples)) # error
  {
    stop('Cannot specify both point.density and num.samples')
  }
  
  



  
  if (check.memory==TRUE)
  {  
    if (verbose==TRUE)
    {
      cat(sprintf("Convex hull calculated with %.0f simplices.\nCalculation of inequality constraints will require allocation of %.0f double-precision numbers.\n", nrow(hull_matrix), nrow(hull_matrix)^2))
    }    
    
    if (ncol(data) > 5)
    {
      warning(sprintf("Algorithm may be very slow on high dimensional data (n>5: here, n=%d)",ncol(data)))
    }
    
    stop('Set check.memory=F to continue.\n')
  }
  
  
  
  if (method=="hitandrun")
  {

    if (verbose==TRUE)
    {
      cat('Calculating linear constraints...')
    }
    constraints <- hullToConstr(data_reduced, convexhull$hull)
    if (verbose==TRUE)
    {
      cat('done.\n')
    }
    
    volume_convexhull <- convexhull$vol
    
    
    if (verbose==TRUE)
    {
      cat(sprintf('Sampling %d random points via hit-and-run, %d per chunk...\n',np, chunksize))
    }
    num.samples.completed <- 0
    num.chunks <- ceiling(np/chunksize)
    samples <- vector(mode="list",length=num.chunks)
    
    pb <- progress_bar$new(total=np)
    if(verbose==TRUE)
    {
      pb$tick(0)
    }
    
    for (i in 1:num.chunks)
    {
      if (verbose==TRUE)
      {
        if (!pb$finished==TRUE)
        {
          pb$update(num.samples.completed/np)
        }
      }
      num.samples.to.take <- min(chunksize, np - num.samples.completed)
      samples_this <- hitandrun::hitandrun(constraints,n.samples=num.samples.to.take)
      samples[[i]] <- samples_this
      num.samples.completed <- num.samples.completed + num.samples.to.take
    }
    samples <- do.call("rbind",samples)
    dimnames(samples) <- list(NULL,dimnames(data)[[2]])
    if (verbose==TRUE)
    {
      pb$terminate()
      cat('\ndone.\n')
    }
    
    hv_hull <- new("Hypervolume",
                   Data=data,
                   RandomPoints= samples,
                   PointDensity=point.density,
                   Volume= volume_convexhull,
                   Dimensionality=ncol(samples),
                   ValueAtRandomPoints=rep(1, nrow(samples)),
                   Name=sprintf("Convex expectation for %s", ifelse(class(input)=="Hypervolume", input@Name, deparse(substitute(data))[1])),
                   Method="Adaptive hit and run convex expectation")	
    
    return(hv_hull)	    
    
    
    
    
    
  }
  else if (method=="rejection")
  {
    # REJECTION SAMPLING FROM BOUNDING BOX TO FILL IN CONVEX POLYTOPE APPROACH  
    ntrials <- chunksize # number of candidate points to propose at a time
    done = FALSE
    inpoints <- NULL
    niter <- 0
    
    if (verbose==TRUE)
    {
      cat(sprintf('Sampling %d random points via rejection...\n',np))
    }
    
    pb <- progress_bar$new(total=np)
    if(verbose==TRUE)
    {
      pb$tick(0)
    }
    
    # generate a bounding box, and test if each point is in the convex hull or not
    while(!done)
    {
      if (verbose==TRUE)
      {
        if (!pb$finished==TRUE)
        {
          if (is.null(inpoints))
          {
            pb$update(0)
          }
          else
          {
            pb$update(nrow(inpoints)/np)
          }
        }
      }
      # try a set of test points
      testpoints <- expectation_box(data_reduced)
      
      chullin <- (inhull_compiled(testpts= testpoints@RandomPoints, calpts=data_reduced, hull=hull_matrix)==1)
      # figure out which are 'in'
      inpoints_temp <- testpoints@RandomPoints[chullin,]
      
      niter <- niter + 1
      inpoints <- rbind(inpoints, inpoints_temp)
      
      done = nrow(inpoints) >= np
    }
    if (verbose==TRUE)
    {
      pb$terminate()
      cat(sprintf('\nTrimming points down to %.0f values.\n',np))
    }    
    inpoints <- inpoints[1:np,]
    dimnames(inpoints) <- list(NULL,dimnames(data)[[2]])
    
    point.density <- nrow(inpoints) / hull_volume
    
    # MAKE a new hypervolume with the convex hull shape and volume
    hv_chull <- new("Hypervolume",
                    Name=sprintf("Convex expectation for %s", ifelse(class(input)=="Hypervolume", input@Name, deparse(substitute(data))[1])),
                    Method="Rejection sampling convex expectation",
                    Data=as.matrix(data),
                    RandomPoints=inpoints, 
                    Dimensionality=ncol(inpoints), 
                    Volume=hull_volume,
                    PointDensity = point.density,
                    Parameters= list(),
                    ValueAtRandomPoints = rep(1, nrow(inpoints))
    )
    
    return(hv_chull)    

  } 
  else
  {
    stop(sprintf('Unrecognized method %s',method))
  }
}