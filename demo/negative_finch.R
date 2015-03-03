if (exists('doHypervolumeNegativeFinchDemo')==TRUE)
{
  data(finch)
  
  # rescale data for four variable across all species
  fd <- scale(finch[,c("BodyL","WingL","TailL","BeakW")])
  
  # get overall community hypervolume
  hv_finch <- hypervolume(fd,bandwidth=estimate_bandwidth(fd))
  hv_finch@Name <- "Isabela Island finches"
  
  # compute convex expectation
  ec_finch <- expectation_convex(hv_finch, check_memory=F)
  
  # find negative features
  nf_finch <- negative_features(hv_finch, ec_finch, set_check_memory=F)
  nf_finch@Name <- "Non-convex negative features"
  
  # fraction of volume occupied
  print(vol_ratio <- nf_finch@Volume / ec_finch@Volume)
  # approximate axis length occupied
  print(lengthratio <- vol_ratio ^ (1/4))
  
  # identification of negative features
  print(apply(nf_finch@RandomUniformPointsThresholded,2,median))
  
  # plot the results
  plot(hypervolume_join(hv_finch, nf_finch),
       npmax_random=5000,
       showdensity=F,
       showdata=T,
       shuffle=F,
       col=c("green","magenta"),
       names=c("Body length","Wing length","Tail length","Beak width")
       )
  
  # make a movie
  open3d(windowRect=c(0,0,600,600))
  plot(hypervolume_join(hv_finch, nf_finch),
       npmax_random=5000,
       showdensity=F,
       showdata=T,
       shuffle=F,
       col=c("green","magenta"),
       names=c("Body length","Wing length","Tail length","Beak width"),
       pair=F,
       cex.random=2,
  )
  movie3d(spin3d(),duration=5,movie='movie_finch',dir='./')
} else
{
  message('Demo does not run by default to meet CRAN runtime requirements.')
  message('This demo requires approximately 1 minutes to run.')  
  message('To run the demo, type')
  message('\tdoHypervolumeNegativeFinchDemo=TRUE')
  message('\tdemo(negative_finch)')
  message('at the R command line prompt.')
}