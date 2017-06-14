\name{get_centroid}
\alias{get_centroid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get centroid of hypervolume or hypervolume list
}
\description{
Returns the column mean of the random points in each hypervolume.
}
\usage{
get_centroid(hv)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv}{
A \code{Hypervolume} or \code{HypervolumeList} object.
}
}
\value{
Either a vector or a matrix of column of centroid values along each axis. 
}
\examples{
\dontrun{
data(iris)
hv = hypervolume_gaussian(iris[,1:2])
get_centroid(hv)
}
}