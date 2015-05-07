if (exists('doHypervolumeHolesFinchDemo')==TRUE)
{
  data(finch)
  
  # convert length units to comparable scales; lo
  fd <- scale(finch[,c("BodyL","WingL","TailL","BeakW")])
  
  getholes <- function(bandwidth_factor)
  {
    # get overall community hypervolume
    hv_finch <- hypervolume(fd,bandwidth=bandwidth_factor*estimate_bandwidth(fd),name=sprintf("Finches b = %.1f*auto",bandwidth_factor))
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
    cat(sprintf("Position of holes: %s\n", paste(apply(holes_finch@RandomUniformPointsThresholded,2,median), collapse=" ")))
    
    # return combined result 
    return(hypervolume_join(hv_finch, ec_finch, holes_finch))
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
         col=c("red","orange","blue"),
         names=c("Body length","Wing length","Tail length","Beak width")
    )  
  }
  
  # calculate holes for two different bandwidth choices
  # auto bandwidth - no holes detected
  result_auto <- getholes(1.0)
  # 70% of automatic bandwidth - several large holes detected
  result_auto07 <- getholes(0.7)

  # show no holes in full-bandwidth result
  plotholes(result_auto)
  # show major holes in smaller-bandwidth result
  plotholes(result_auto07)
} else
{
  message('Demo does not run by default to meet CRAN runtime requirements.')
  message('This demo requires approximately 1 minutes to run.')  
  message('To run the demo, type')
  message('\tdoHypervolumeHolesFinchDemo=TRUE')
  message('\tdemo(holes_finch)')
  message('at the R command line prompt.')
}