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
data(penguins,package='palmerpenguins')
penguins_no_na = as.data.frame(na.omit(penguins))
penguins_adelie = penguins_no_na[penguins_no_na$species=="Adelie",
                    c("bill_length_mm","bill_depth_mm","flipper_length_mm")]
hv = hypervolume_gaussian(penguins_adelie)
get_centroid(hv)
}
}