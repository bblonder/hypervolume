\name{hypervolume_project}
\alias{hypervolume_project}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Geographical projection of hypervolume for species distribution modeling, using the hypervolume as the environmental niche model.
}
\description{
Determines a suitability score by calculating the hypervolume value at each of a set of points in an input raster stack based on either a probability density estimation or inclusion test.

Note that projected values are not normalized and are not necessarily constrained to fall between 0 and 1.
}
\usage{
hypervolume_project(hv, rasters, type = "probability", verbose = TRUE,
                 ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv}{
An input hypervolume
}
  \item{rasters}{
A \code{RasterStack} with the same names as the dimension names of the hypervolume.
}
  \item{type}{
If \code{'probability'}, suitability scores correspond to probability density values estimated using \code{\link{hypervolume_estimate_probability}}; if \code{'inclusion'}, scores correspond to binary presence/absence values from calling \code{\link{hypervolume_inclusion_test}}.
}
  \item{...}{
Additional arguments to be passed to either \code{\link{hypervolume_estimate_probability}} or \code{\link{hypervolume_inclusion_test}}.
}
\item{verbose}{
If \code{TRUE}, prints diagnostic and progress output.
}

}

\value{
A \code{raster} object of same resolution and extent as the input layers corresponding to suitability values.
}

\seealso{
\code{\link{hypervolume_estimate_probability}}, \code{\link{hypervolume_inclusion_test}}
}
\examples{
# example does not run to meet CRAN runtime guidelines - set TRUE to run
hypervolume_project_demo = FALSE
if (hypervolume_project_demo==TRUE)
{
  # load in lat/lon data
  data('quercus') 
  data_alba = subset(quercus, Species=="Quercus alba")[,c("Longitude","Latitude")]
  data_alba = data_alba[sample(1:nrow(data_alba),500),]
   
  # get worldclim data from internet
  require(maps)
  require(raster)
  climatelayers = getData('worldclim', var='bio', res=10, path=tempdir())
    
  # z-transform climate layers to make axes comparable
  climatelayers_ss = climatelayers[[c(1,12)]]
  for (i in 1:nlayers(climatelayers_ss))
  {
    climatelayers_ss[[i]] <- 
    	(climatelayers_ss[[i]] - cellStats(climatelayers_ss[[i]], 'mean')) / 
    	cellStats(climatelayers_ss[[i]], 'sd') 
  }
  climatelayers_ss = crop(climatelayers_ss, extent(-150,-50,15,60))
  
  # extract transformed climate values
  climate_alba = extract(climatelayers_ss, data_alba[1:300,])
  
  # compute hypervolume
  hv_alba <- hypervolume_gaussian(climate_alba)
  
  # do geographical projection
  raster_alba_projected_accurate <- hypervolume_project(hv_alba, 
                                      rasters=climatelayers_ss)
  raster_alba_projected_fast = hypervolume_project(hv_alba, 
                                      rasters=climatelayers_ss, 
                                      type='inclusion',
                                      fast.or.accurate='fast')
  
  # draw map of suitability scores
  plot(raster_alba_projected_accurate,xlim=c(-100,-60),ylim=c(25,55))
  map('usa',add=TRUE)
  
  plot(raster_alba_projected_fast,xlim=c(-100,-60),ylim=c(25,55))
  map('usa',add=TRUE)
}
}
