bindata_padded <- function(data, bin.widths, num.shifts, verbose=TRUE)
{	
  data_cut <- as.matrix(data)
  num_dims <- ncol(data)
  
  if (length(bin.widths)==1)
  {
    bin.widths = rep(bin.widths, num_dims)
  }
  if (verbose == TRUE)
  {
    cat("Defining shifts\n")
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
    data_cut_index[,i] <- as.numeric(cut(data_cut[,i],breaks=cuts[[i]], ordered=TRUE,include.lowest=TRUE))
  }
  
  # reduce to unique cases
  data_cut_index_shifted <- unique(data_cut_index)
  
  # recursively shift the larger dataset
  for (i in setdiff(0:num.shifts,0))
  {
    if (verbose == TRUE)
    {
      cat(sprintf('Dataset shifts %d - %d to complete\n', i, nrow(grid_shifts)))
    }
    # do shifts in all directions in indexed coordinates
    data_cut_index_shifted_list <- lapply(1:nrow(grid_shifts), function(index) { 
      if (verbose == TRUE)
      {
        cat('.')
        if (i %% 80 == 0)
        {
          cat('\n')
        }
      }
      return(t(t(data_cut_index_shifted) + grid_shifts[index,])) 
    })
    if (verbose == TRUE)
    {
      cat('done.\n')
    }
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

expectation_adaptive_box <- function(data, bin.widths, num.shifts=1, density=10^(3+ncol(data)), verbose=TRUE) # eventually convert num.shifts to a 
{
  binneddata <- bindata_padded(data, bin.widths=bin.widths, num.shifts=num.shifts,verbose=verbose)
  
  # take only filled unique grid cells
  data_cut <- unique(binneddata$data_cut)
  
  cuts = binneddata$cuts
  
  nbins.filled <- nrow(data_cut)
  nbins.possible <- nrow(cuts)^ncol(cuts)
  
  bin_fraction = nbins.filled / nbins.possible # number of bins that we found divided by possible number of bins
  if (verbose == TRUE)
  {
    cat(sprintf('Exploring %d / %.0f bins\n',nbins.filled, nbins.possible))
  }
  rp_all <- vector(mode="list",length=nrow(data_cut))
  volume_all <- rep(NA, nrow(data_cut))
  
  if (verbose == TRUE)
  {
    cat(sprintf('Random sampling within %d hyperboxes\n', nrow(data_cut)))
  }
  for (i in 1:nrow(data_cut))
  {	
    if (verbose == TRUE)
    {
      cat('.')
      if (i %% 80 == 0)
      {
        cat('\n')
      }
    }
    lower <- sapply(1:length(cuts),function(dim) {cuts[[dim]][data_cut[i,dim]]})
    higher <- sapply(1:length(cuts),function(dim) {cuts[[dim]][data_cut[i,dim]+1]})
    
    volume <- prod(higher - lower) # hyperbox volume of this grid cell
    np <- ceiling(density*volume)
    
    
    randpoints <- replicate(length(cuts),runif(n=np)) # get unit box of random points
    
    randpoints <- matrix(randpoints, ncol=ncol(data))
    
    for (column_index in 1:ncol(randpoints)) # rescale random points to right locations
    {
      randpoints[, column_index] <- randpoints[, column_index] * (higher[column_index] - lower[column_index]) + lower[column_index]
    }
    
    rp_all[[i]] <- randpoints
    volume_all[i] <- volume
  }
  if (verbose == TRUE)
  {
    cat('done\n')
  }
  
  if (verbose == TRUE)
  {
    cat(sprintf('Binding rows...'))
  }
  rp_all <- do.call("rbind",rp_all) # bind all together
  if (verbose == TRUE)
  {
    cat(sprintf('done\n'))
  }
  volume_all <- sum(volume_all)
  dimnames(rp_all) <- list(NULL, dimnames(data_cut)[[2]])
  
  result <- new("Hypervolume",
                Data=as.matrix(data),
                Method="Expectation (adaptive box)",
                RandomUniformPointsThresholded=rp_all,
                PointDensity=density,
                Volume= volume_all,
                Dimensionality=ncol(data),
                ProbabilityDensityAtRandomUniformPoints=normalize_probability(rep(1/nrow(rp_all),nrow(rp_all)), density),
                Name=sprintf("Expectation - adaptive box (%d/%d bins)",nbins.filled, nbins.possible))
  
  return(result)
}

hypervolume_svm <- function(data, name=NULL, verbose=TRUE, output.density=10^ncol(data), expectation.num.shifts=1, expectation.bin.widths=2*estimate_bandwidth(data), svm.nu=0.01, svm.gamma=0.5, svm.chunksize=1e4)
{
  data <- as.matrix(data)
  
  expectation <- expectation_adaptive_box(data, density= output.density, num.shifts=expectation.num.shifts, bin.widths=expectation.bin.widths,verbose=verbose)
  
  if (svm.nu <= 0)
  {
    stop("Parameter svm.nu must be positive")
  }
  
  if (svm.gamma <= 0)
  {
    stop("Parameter svm.gamma must be positive")
  }
  
  if (svm.chunksize <= 100 | svm.chunksize > 1e7)
  {
    stop("Parameter svm.chunksize may be in range that will produce poor performance.")
  }
  if (verbose == TRUE)
  {
    cat('Building support vector machine model')
  }
  svm.model<-e1071::svm(data,
                        y=NULL,
                        type='one-classification',
                        nu= svm.nu,
                        gamma= svm.gamma,
                        scale=TRUE,
                        kernel="radial")	
  if (verbose == TRUE)
  {
    cat('...done\n')
  }
  
  np <- nrow(expectation@RandomUniformPointsThresholded)
  
  if (verbose == TRUE)
  {
    cat(sprintf('Sampling %d random points from support vector machine, %d per chunk.\n',np, svm.chunksize))
  }
  num.samples.completed <- 0
  num.chunks <- ceiling(np/svm.chunksize)
  svm.probs <- vector(mode="list",length=num.chunks)
  for (i in 1:num.chunks)
  {
    if (verbose == TRUE)
    {
      cat(sprintf('Running chunk %d / %d\n', i, num.chunks))
    }
    num.samples.to.take <- min(svm.chunksize, np - num.samples.completed)
    
    inputdata.this <- expectation@RandomUniformPointsThresholded[num.samples.completed:(num.samples.completed+num.samples.to.take-1),]
    
    svm.pred.this <- predict(svm.model,inputdata.this)
    svm.pred.this.pos <- inputdata.this[svm.pred.this==TRUE,]
    
    svm.probs[[i]] <- svm.pred.this.pos
    num.samples.completed <- num.samples.completed + num.samples.to.take
  }
  svm.probs <- do.call("rbind",svm.probs)
  dimnames(svm.probs) <- list(NULL,dimnames(data)[[2]])
  if (verbose == TRUE)
  {
    cat('done.\n')	
  }
  
  numpoints_svm <- nrow(svm.probs)
  vol_svm <- numpoints_svm / expectation@PointDensity
  
  hv_svm <- new("Hypervolume",
                Data=data,
                Method = "svm",
                RandomUniformPointsThresholded= svm.probs,
                PointDensity= expectation@PointDensity,
                Volume= vol_svm,
                Dimensionality=ncol(svm.probs),
                ProbabilityDensityAtRandomUniformPoints=normalize_probability(rep(1,nrow(svm.probs)),expectation@PointDensity),
                Name=ifelse(is.null(name), "untitled", toString(name)),
                Parameters=c(svm.nu=svm.nu, svm.gamma=svm.gamma, expectation.num.shifts=expectation.num.shifts, expectation.bin.widths=expectation.bin.widths))
  
  if (nrow(hv_svm@RandomUniformPointsThresholded) < 10^ncol(data))
  {
    warning(sprintf("Hypervolume is represented by a low number of random points (%d) - suggested minimum %d.\nConsider increasing point density to improve accuracy.",nrow(hv_svm@RandomUniformPointsThresholded),10^ncol(data)))
  }
  
  return(hv_svm)	
}
