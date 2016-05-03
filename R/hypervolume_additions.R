bindata_padded <- function(data, bin.widths, num.shifts)
{	
	data_cut <- as.matrix(data)
	num_dims <- ncol(data)
	
	if (length(bin.widths)==1)
	{
		bin.widths = rep(bin.widths, num_dims)
	}

	# figure out every possible shift in every direction
	grid_shifts <- as.matrix(expand.grid(rep(list(c(-1,0,1)),num_dims)))
	
	data_cut_index <- data_cut
	
	# slice data into grid cells
	cuts <- vector(mode="list",length=num_dims)
	names(cuts) <- dimnames(data_cut)[[2]]
	for (i in 1:ncol(data_cut))
	{
		cuts[[i]] <- seq(min(data_cut[,i])-bin.widths[i], max(data_cut[,i])+bin.widths[i],by=bin.widths[i])
		data_cut_index[,i] <- as.numeric(cut(data_cut[,i],breaks=cuts[[i]], ordered=T,include.lowest=T))
	}
	
	# reduce to unique cases
	data_cut_index_shifted <- unique(data_cut_index)

	# recursively shift the larger dataset
	for (i in setdiff(0:num.shifts,0))
	{
		cat(sprintf('shift %d - %d to complete', i, nrow(grid_shifts)))
		# do shifts in all directions in indexed coordinates
		data_cut_index_shifted_list <- lapply(1:nrow(grid_shifts), function(index) { 
			cat('.')
			return(t(t(data_cut_index_shifted) + grid_shifts[index,])) 
		})
		cat('done\n')
		data_cut_index_shifted <- do.call("rbind", data_cut_index_shifted_list)
		# reduce again to unique cases
		data_cut_index_shifted <- unique(data_cut_index_shifted)
		# shift indices up by one to account for new lowest position
		data_cut_index_shifted <- data_cut_index_shifted + 1
		
		# add data to the cut ranges
		# add row at top and bottom for next lowest cut
		cuts <- lapply(1:num_dims, function(index) { c(cuts[[index]][1]-bin.widths[index], cuts[[index]], cuts[[index]][length(cuts[[index]])]+bin.widths[index]) })
	}

	# keep original dimension names
	dimnames(data_cut_index_shifted) <- list(NULL, dimnames(data)[[2]])

	return(list(data_cut= data_cut_index_shifted, cuts=cuts))
}

expectation_adaptive_box <- function(data, bin.widths, num.shifts=1, density=10^ncol(data)) # eventually convert num.shifts to a 
{
	binneddata <- bindata_padded(data, bin.widths=bin.widths, num.shifts=num.shifts)
	
	# take only filled unique grid cells
	data_cut <- unique(binneddata$data_cut)
	
	cuts = binneddata$cuts
	
	nbins.filled <- nrow(data_cut)
	nbins.possible <- nrow(cuts)^ncol(cuts)
	
	bin_fraction = nbins.filled / nbins.possible # number of bins that we found divided by possible number of bins
	
	cat(sprintf('Exploring %d / %.0f bins\n',nbins.filled, nbins.possible))
	
	rp_all <- vector(mode="list",length=nrow(data_cut))
	volume_all <- rep(NA, nrow(data_cut))
	
	cat(sprintf('Random sampling within %d hyperboxes', nrow(data_cut)))
	for (i in 1:nrow(data_cut))
	{	
		cat('.')
		lower <- sapply(1:length(cuts),function(dim) {cuts[[dim]][data_cut[i,dim]]})
		higher <- sapply(1:length(cuts),function(dim) {cuts[[dim]][data_cut[i,dim]+1]})
		
		volume <- prod(higher - lower) # hyperbox volume of this grid cell
		np <- ceiling(density*volume)
	
		randpoints <- replicate(length(cuts),runif(n=np)) # get unit box of random points
		for (column_index in 1:ncol(randpoints)) # rescale random points to right locations
		{
			randpoints[, column_index] <- randpoints[, column_index] * (higher[column_index] - lower[column_index]) + lower[column_index]
		}
		
		rp_all[[i]] <- randpoints
		volume_all[i] <- volume
	}
	cat('done\n')
	cat(sprintf('Binding rows...'))
	rp_all <- do.call("rbind",rp_all) # bind all together
	cat(sprintf('done\n'))
	volume_all <- sum(volume_all)
	dimnames(rp_all) <- list(NULL, dimnames(data_cut)[[2]])
	
	result <- new("Hypervolume",
				Data=as.matrix(data),
				RandomUniformPointsThresholded=rp_all,
				PointDensity=density,
				Volume= volume_all,
				Dimensionality=ncol(data),
				ProbabilityDensityAtRandomUniformPoints=rep(1/nrow(rp_all),nrow(rp_all)),
				Name=sprintf("Expectation - adaptive box (%d/%d bins)",nbins.filled, nbins.possible))
				
	return(result)
}

hypervolume_svm <- function(data, output.density=10^ncol(data), expectation.num.shifts=1, expectation.bin.widths=2*estimate_bandwidth(data), svm.nu=0.01, svm.gamma=0.5, svm.chunksize=1e4)
{
	data <- as.matrix(data)
	
	expectation <- expectation_adaptive_box(data, density= output.density, num.shifts=expectation.num.shifts, bin.widths=expectation.bin.widths)
	
	cat('running SVM')
	svm.model<-e1071::svm(data,
	           y=NULL,
	           type='one-classification',
	           nu= svm.nu,
	           gamma= svm.gamma,
	           scale=TRUE,
	           kernel="radial")	
	cat('...done\n')
	
	np <- nrow(expectation@RandomUniformPointsThresholded)
	cat(sprintf('Sampling %d random points from support vector machine, %d per chunk...\n',np, svm.chunksize))
	num.samples.completed <- 0
	num.chunks <- ceiling(np/svm.chunksize)
	svm.probs <- vector(mode="list",length=num.chunks)
	for (i in 1:num.chunks)
	{
	  cat('.')
	  num.samples.to.take <- min(svm.chunksize, np - num.samples.completed)
	  
	  inputdata.this <- expectation@RandomUniformPointsThresholded[num.samples.completed:(num.samples.completed+num.samples.to.take-1),]
	  
	  svm.pred.this <- predict(svm.model,inputdata.this)
	  svm.pred.this.pos <- inputdata.this[svm.pred.this==TRUE,]
	  
	  svm.probs[[i]] <- svm.pred.this.pos
	  num.samples.completed <- num.samples.completed + num.samples.to.take
	}
	svm.probs <- do.call("rbind",svm.probs)
	dimnames(svm.probs) <- list(NULL,dimnames(data)[[2]])
	cat('done.\n')	
	
	numpoints_svm <- nrow(svm.probs)
	vol_svm <- numpoints_svm / expectation@PointDensity
	
	hv_svm <- new("Hypervolume",
				Data=data,
				RandomUniformPointsThresholded= svm.probs,
				PointDensity= expectation@PointDensity,
				Volume= vol_svm,
				Dimensionality=ncol(svm.probs),
				ProbabilityDensityAtRandomUniformPoints=rep(1/nrow(svm.probs),nrow(svm.probs)),
				Name="Support vector machine")	
				
	return(hv_svm)	
}




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

	

expectation_convex_hitandrun <- function(data, point_density, hitandrun.chunksize=1e5) # larger values better but need more memory for the initial chain discarding
{
	data <- as.matrix(data)
	cat('Calculating convex hull...')
	hull <- convhulln(data,options="FA")
	cat('done.\n')
	
	cat('Calculating linear constraints...\n')
	constraints <- hullToConstr(data, hull$hull)
	cat('done.\n')
	
	volume_convexhull <- hull$vol
	
	np <- ceiling(point_density * volume_convexhull)
	
	cat(sprintf('Sampling %d random points via hit-and-run, %d per chunk...\n',np, hitandrun.chunksize))
	num.samples.completed <- 0
	num.chunks <- ceiling(np/hitandrun.chunksize)
	samples <- vector(mode="list",length=num.chunks)
	for (i in 1:num.chunks)
	{
	  cat('.')
	  num.samples.to.take <- min(hitandrun.chunksize, np - num.samples.completed)
	  samples_this <- hitandrun::hitandrun(constraints,n.samples=num.samples.to.take)
	  samples[[i]] <- samples_this
	  num.samples.completed <- num.samples.completed + num.samples.to.take
	}
	samples <- do.call("rbind",samples)
	dimnames(samples) <- list(NULL,dimnames(data)[[2]])
	cat('done.\n')
	
	
	hv_hull <- new("Hypervolume",
				Data=data,
				RandomUniformPointsThresholded= samples,
				PointDensity=point_density,
				Volume= volume_convexhull,
				Dimensionality=ncol(samples),
				ProbabilityDensityAtRandomUniformPoints=rep(1/nrow(samples),nrow(samples)),
				Name="Convex hull")	
				
	return(hv_hull)	
}


mvnorm <- function(x, mu, Sigma, k=length(x),fast=F,diagonalSigma=F) # fast does quadratic approximation
{
	xminusmu <- as.numeric(x - mu)
	if (diagonalSigma==TRUE)
	{
		Sigmainv = 1/Sigma
		Sigmadet = Sigma[diag(1,nrow(Sigma),ncol(Sigma))]	
	}
	else
	{
		Sigmainv <- solve(Sigma)
		Sigmadet <- det(Sigma)
	}
	
	argument <- (-1/2) * (t(xminusmu) %*% Sigmainv %*% (xminusmu))  

	
	if (fast==TRUE)
	{
		exparg <- 1 + argument
		if (exparg < 0)
		{
			exparg <- 0
		}
	}
	else
	{
		exparg <- exp(argument)
	}
	
 	prefactor=1/sqrt(((2*pi)^k)*Sigmadet)
	
	return(as.numeric(exparg*prefactor))
}

# choose a probability value equal to the prob density of a single point at 1 s.d. distance in dimensionless coordinates
estimate_threshold_gaussian <- function(sd.count=1, dim)
{
	 mvnorm(rep(sd.count,dim),rep(0,dim),Sigma=diag(1,dim,dim))
}

hypervolume_gaussian <- function(data, output.density=10^(ncol(data)), expectation.num.shifts=1, expectation.bin.widths=2*estimate_bandwidth(data), kde.bandwidth=estimate_bandwidth(data)/2, kde.chunksize=1e4, output.threshold= estimate_threshold_gaussian(sd.count=1,dim=ncol(data)))
{
	expectation <- expectation_adaptive_box(data, density= output.density, num.shifts=expectation.num.shifts, bin.widths=expectation.bin.widths)
	
	cat('running KDE...\n')

	if (ncol(data)==1) # kde package defined bandwidth as s.d. in 1 dim and var in >= 2 dim...
	{
	  Hmatrix <- diag(1,nrow=ncol(data),ncol=ncol(data))*(kde.bandwidth)
	}
	else
	{
	  Hmatrix <- diag(1,nrow=ncol(data),ncol=ncol(data))*(kde.bandwidth)^2
	}

	np <- nrow(expectation@RandomUniformPointsThresholded)
	cat(sprintf('Sampling %d random points from kernel density estimate, %d per chunk...\n',np, kde.chunksize))
	num.samples.completed <- 0
	num.chunks <- ceiling(np/kde.chunksize)
	kde.probs <- vector(mode="list",length=num.chunks)
	for (i in 1:num.chunks)
	{
	  cat('.')
	  num.samples.to.take <- min(kde.chunksize, np - num.samples.completed)
	  kde.probs.this <- ks::kde(x=data, 
	                               H=Hmatrix,
	                               eval.points= expectation@RandomUniformPointsThresholded[num.samples.completed:(num.samples.completed+num.samples.to.take-1),],
	                               verbose=T)$estimate
	  kde.probs[[i]] <- kde.probs.this
	  num.samples.completed <- num.samples.completed + num.samples.to.take
	}
	kde.probs <- do.call("c",kde.probs)
	cat('done.\n')
	
	points_thresholded <- expectation@RandomUniformPointsThresholded[kde.probs > output.threshold, , drop=F]
	probs_thresholded <- kde.probs[kde.probs > output.threshold]
	
	numpoints_kde <- nrow(points_thresholded)
	vol_kde <- numpoints_kde / expectation@PointDensity
	
	hv_kde <- new("Hypervolume",
				Data=as.matrix(data),
				RandomUniformPointsThresholded= points_thresholded,
				PointDensity= expectation@PointDensity,
				Volume= vol_kde,
				Dimensionality=ncol(data),
				ProbabilityDensityAtRandomUniformPoints= probs_thresholded,
				Name="Gaussian kernel density estimate")	
				
	return(hv_kde)	
}

hypervolume_threshold <- function(hv, thresholds=NULL, plot=TRUE)
{
	if (is.null(thresholds))
	{
		thresholds <- seq(min(hv@ProbabilityDensityAtRandomUniformPoints), max(hv@ProbabilityDensityAtRandomUniformPoints), length.out=100)
	}
	#stopifnot (!any(thresholds < min(hv@ ProbabilityDensityAtRandomUniformPoints)))
	result <- vector(mode="list",length=length(thresholds))
	for (i in 1:length(thresholds))
	{
		hv_new <- hv
		hv_new@RandomUniformPointsThresholded <- hv@RandomUniformPointsThresholded[hv@ProbabilityDensityAtRandomUniformPoints > thresholds[i],,drop=F]
		hv_new@Volume <- length(which(hv@ProbabilityDensityAtRandomUniformPoints > thresholds[i])) / hv@PointDensity
		
		result[[i]] <- hv_new
	}
	names(result) <- thresholds
	
	volumes <- sapply(result, function(x) {x@Volume})	
	
	if (plot==TRUE)
	{
	  plot(volumes~thresholds,type='l',xlab="Threshold value",ylab="Hypervolume")
	  axis(side=3,at=seq(thresholds[1],thresholds[length(thresholds)],length.out=10),labels=seq(1,length(thresholds),length.out=10))
	  mtext(side=3,line=3,"Threshold index")
	  abline(v= estimate_threshold_gaussian(sd.count=1,dim=hv@Dimensionality),col='red')
	  abline(v= estimate_threshold_gaussian(sd.count=2,dim=hv@Dimensionality),col='blue')
	  legend('topright',c("Probability at 1 s.d. distance", "Probability at 2 s.d. distance"),col=c("red","blue"),cex=0.5,lty=1)
	}
	
	return(list(
	  HypervolumesThresholded=result,
	  Volumes=data.frame(threshold=thresholds,volume=volumes)
	  ))
}