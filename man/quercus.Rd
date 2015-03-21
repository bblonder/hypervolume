\name{quercus}
\alias{quercus}
\docType{data}
\title{
Data and demo for Quercus (oak) tree distributions
}
\description{
Data for occurrences of Quercus alba and Quercus rubra based on geographic observations. Demonstration analysis of how to use hypervolumes for species distribution modeling using WorldClim data.
}
\usage{data(quercus)}
\format{
  A data frame with 3779 observations on the following 3 variables.
  \describe{
    \item{\code{Species}}{a factor with levels \code{Quercus alba} \code{Quercus rubra}}
    \item{\code{Latitude}}{a numeric vector}
    \item{\code{Longitude}}{a numeric vector}
  }
}

\source{
Occurrence data come from the BIEN2 database (http://bien.nceas.ucsb.edu/bien/). Climate data are from WorldClim.
}
\references{
Blonder, B., Lamanna, C., Violle, C., Enquist, B. The n-dimensional hypervolume. Global Ecology and Biogeography (2014)
}
\examples{
demo('quercus', package='hypervolume')
}