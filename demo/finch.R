if (exists('doHypervolumeFinchDemo')==TRUE)
{
  data(finch)
  
  species_list = unique(finch$Species)
  num_species = length(species_list)
  
  hv_finches_list = new("HypervolumeList")
  hv_finches_list@HVList = vector(mode="list",length=num_species)
  
  # compute hypervolumes for each species
  for (i in 1:num_species)
  {
    this_species = subset(finch, Species==species_list[i])
    # keep the trait data
    this_species_log <- log10(this_species[,2:ncol(this_species)])
    # make a hypervolume using auto-bandwidth
    hv_finches_list@HVList[[i]] <- hypervolume(this_species_log, bandwidth=estimate_bandwidth(this_species_log),
                                               reps=10000, quantile=0, name=as.character(species_list[i]), warn=FALSE)
  }
  
  # compute all pairwise overlaps
  overlap = matrix(NA, nrow=num_species, ncol=num_species)
  dimnames(overlap)=list(species_list, species_list)
  for (i in 1:num_species)
  {
    for (j in i:num_species)
    {
      if (i!=j)
      {
        # compute set operations on each pair
        this_set = hypervolume_set(hv_finches_list@HVList[[i]], hv_finches_list@HVList[[j]], check_memory=FALSE)
        # calculate a Sorensen overlap index (2 x shared volume / sum of |hv1| + |hv2|)
        overlap[i,j] = 2 * this_set@HVList$Intersection@Volume / (hv_finches_list@HVList[[i]]@Volume + hv_finches_list@HVList[[j]]@Volume)
      }
    }   
  }
  

  
  # show all hypervolumes
  plot(hv_finches_list,npmax.random=500,darkfactor=0.5,cex.legend=0.25,cex.names=0.75)
  
  # show pairwise overlaps - note that actually very few species overlap in nine dimensions
  op <- par(mar=c(10,10,1,1))
  image(x=1:nrow(overlap), y=1:nrow(overlap), z=overlap,axes=F,xlab='',ylab='',col=rainbow(5))
  box()
  axis(side=1, at=1:(length(dimnames(overlap)[[1]])),dimnames(overlap)[[1]],las=2,cex.axis=0.75)
  axis(side=2, at=1:(length(dimnames(overlap)[[2]])),dimnames(overlap)[[2]],las=1,cex.axis=0.75)
  par(op)
  
  rm(doHypervolumeFinchDemo)
} else
{
  message('Demo does not run by default to meet CRAN runtime requirements.')
  message('This demo requires approximately 3 minutes to run.')  
  message('To run the demo, type')
  message('\tdoHypervolumeFinchDemo=TRUE')
  message('\tdemo(finch)')
  message('at the R command line prompt.')
}