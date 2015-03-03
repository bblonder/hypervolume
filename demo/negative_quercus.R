if (exists('doHypervolumeNegativeQuercusDemo')==TRUE)
{
  require(raster)
  require(maps)
  
  # load in species data
  data(quercus)
  qalba <- subset(quercus,Species=="Quercus alba")[,c("Longitude","Latitude")]
  
  # get 
  climatelayers_raw <- getData('worldclim', var='bio', res=10, path=tempdir())
  
  # z-transform climate layers to make axes comparable
  climatelayers <- climatelayers_raw[[c(1,4,12,15)]]
  for (i in 1:nlayers(climatelayers))
  {
    climatelayers[[i]] <- (climatelayers[[i]] - cellStats(climatelayers[[i]], 'mean')) / cellStats(climatelayers[[i]], 'sd') 
  }
  
  # extract transformed climate values
  climate_qalba <- extract(climatelayers, qalba)
  climate_all <- na.omit(getValues(climatelayers))
  
  # resample the climate data to a lower density
  climate_all <- climate_all[sample(x=1:nrow(climate_all),size=nrow(climate_qalba)),]
  
  # get the realized niche of quercus alba
  hv_qalba <- hypervolume(climate_qalba,bandwidth=estimate_bandwidth(climate_qalba))
  hv_qalba@Name <- "Q. alba"
  
  # get the 'potential niche' of quercus alba
  hv_qalba_convex <- expectation_convex(climate_qalba, check_memory=F)
  hv_qalba_convex@Name <- "Convex expectation of Q. alba"
  
  # get the nonconvex negative features in q. alba
  hv_nonconvex <- negative_features(hv_obs=hv_qalba, hv_exp=hv_qalba_convex, set_npoints_max=1e4, set_check_memory=F)
  hv_nonconvex@Name <- 'Q. alba non-convex negative features'
  
  # get the 'available space' of north america
  # lower density sampling is OK because we will intersect this with a low density negative features hypervolume
  hv_climatespace <- hypervolume(climate_all,bandwidth=estimate_bandwidth(climate_all), reps=1000)
  hv_climatespace@Name <- "Global climate space"
  
  # compute intersection of the non-convex features with the overall climatespace
  setops <- hypervolume_set(hv_nonconvex, hv_climatespace, npoints_max=1e8, check_memory=F)
  hv_climatelimitation <- setops@HVList$Unique_1
  hv_climatelimitation@Name <- 'Q. alba climate limitation'
  
  # plot the final results
  plot(hypervolume_join(hv_qalba, hv_nonconvex, hv_climatespace, hv_climatelimitation),
       showdata=F,showdensity=F,
       names=c("Mean annual\ntemperature","Temperature\nseasonality","Mean annual\nprecipitation","Precipitation\nseasonality"),
       col=c("red","orange","blue","purple")
       )
  
  
  
  # make a movie
  open3d(windowRect=c(0,0,600,600))
  plot(hypervolume_join(hv_qalba, hv_nonconvex, hv_climatespace, hv_climatelimitation),
       showdata=F,showdensity=F,
       names=c("Mean annual\ntemperature","Temperature\nseasonality","Mean annual\nprecipitation","Precipitation\nseasonality"),
       col=c("red","orange","blue","purple"),
       pair=F,
       cex.random=2,
  )
  movie3d(spin3d(),duration=5,movie='movie_quercus',dir='./')
} else
{
  message('Demo does not run by default to meet CRAN runtime requirements.')
  message('This demo requires approximately 3 minutes to run.')  
  message('To run the demo, type')
  message('\tdoHypervolumeNegativeQuercusDemo=TRUE')
  message('\tdemo(negative_quercus)')
  message('at the R command line prompt.')
}