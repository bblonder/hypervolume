\name{hypervolume_project}
\alias{hypervolume_project}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Geographical projection of hypervolume for species distribution modeling, using the hypervolume as the environmental niche model.
}
\description{
Determines a suitability score by calculating the hypervolume value at each of a set of points in an input raster stack based on either a probability density estimation or inclusion test.
}
\usage{
hypervolume_project(hv, 
	rasters, do.probability = TRUE, normalize = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv}{
An input hypervolume
}
  \item{rasters}{
A \code{RasterStack} with the same names as the dimension names of the hypervolume.
}
  \item{do.probability}{
If \code{TRUE}, suitability scores correspond to probability density values; if \code{FALSE}, scores correspond to binary presence/absence values from inclusion test application.
}
  \item{normalize}{
If \code{TRUE}, suitability scores are normalized between 0 and 1.
}
}

\value{
A \code{raster} object of same resolution and extent as the input layers corresponding to suitability values.
}

\seealso{
\code{\link{hypervolume_estimate_probability}}, \code{\link{hypervolume_inclusion_test}}
}
\examples{
# load in lat/lon data
data('quercus') 
data_alba = subset(quercus, Species=="Quercus alba")[,c("Longitude","Latitude")]
data_alba = data_alba[sample(1:nrow(data_alba),500),]
 
# get worldclim data from internet
require(maps)
require(raster)
climatelayers = getData('worldclim', var='bio', res=10, path=tempdir())
  
# z-transform climate layers to make axes comparable
climatelayers_ss = climatelayers[[c(1,4,12,15)]]
for (i in 1:nlayers(climatelayers_ss))
{
  climatelayers_ss[[i]] <- 
  	(climatelayers_ss[[i]] - cellStats(climatelayers_ss[[i]], 'mean')) / 
  	cellStats(climatelayers_ss[[i]], 'sd') 
}

# extract transformed climate values
climate_alba = extract(climatelayers_ss, data_alba[1:300,])

# compute hypervolume
hv_alba <- hypervolume(climate_alba)

# do geographical projection
raster_alba_projected <- hypervolume_project(hv_alba, rasters=climatelayers_ss)

# draw map of suitability scores
plot(raster_alba_projected,xlim=c(-100,-60),ylim=c(25,55))
map('usa',add=TRUE)

}
