if (exists('doHypervolumeHolesFinchDemo')==TRUE)
{
  data(finch)
  
  # convert length units to comparable scales; lo
  fd <- scale(finch[,c("BodyL","WingL","TailL","BeakW")])
  
  getholes <- function(bw)
  {
    # get overall community hypervolume
    hv_finch <- hypervolume(fd,bandwidth=bw,name='Isabela Island')
    # compute convex expectation
    ec_finch <- expectation_convex(hv_finch, check_memory=F)
    # find holes
    holes_finch <- hypervolume_holes(hv_finch, ec_finch, set_check_memory=F)
    
    # fraction of volume occupied
    vol_ratio <- holes_finch@Volume / ec_finch@Volume
    
    cat(sprintf("Volume of hole / Volume of expectation: %.3f\n", vol_ratio))
    # approximate axis length occupied
    cat(sprintf("Approximate axis fraction occupied: %.3f\n", vol_ratio ^ (1/ncol(fd)))) 
    # identification of location of hole centroids
    cat(sprintf("Position of holes: %s\n", paste(apply(holes_finch@RandomUniformPointsThresholded,2,mean), collapse=" ")))
    
    # return combined result 
    return(hypervolume_join(hv_finch, holes_finch))
  }  
  
  plotholes <- function(result)
  {
    plot(result,
         legend=F,
         showcentroid=F,
         showcontour=T,
         contour.lwd=2,
         reshuffle=F,
         showdensity=F,
         npmax_random=2000,
         showdata=F,
         col=c("#33A02CFF","#6A3D9A80"),
         names=c("Body length","Wing length","Tail length","Beak width")
    )  
  }
  
  # below code takes ~ 2 hours to run - for demo purposes, using the output as-is
  #bw_plugin <- estimate_bandwidth(fd, method='plug-in')
  bw_plugin <- c(0.6150874,0.5639092,0.6771332, 0.4939054)
  
  result <- getholes(bw_plugin)

  # show holes
  plotholes(result)
} else
{
  message('Demo does not run by default to meet CRAN runtime requirements.')
  message('This demo requires approximately 1 minutes to run.')  
  message('To run the demo, type')
  message('\tdoHypervolumeHolesFinchDemo=TRUE')
  message('\tdemo(holes_finch)')
  message('at the R command line prompt.')
}