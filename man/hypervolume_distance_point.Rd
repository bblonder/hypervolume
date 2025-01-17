\name{hypervolume_distance_point}
\alias{hypervolume_distance_point}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Distance from a point to the margin of a hypervolume.
}
\description{
Calculates the distance between two hypervolumes either defined as the Euclidean distance between centroids or as the minimum Euclidean distance between the random points comprising either hypervolume.

Code by Clement Violet.
}
\usage{
hypervolume_distance_point(hv1, x, type = "minimum", 
  num.points.max = 1000, check.memory = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv1}{
A \code{Hypervolume} object.
}
  \item{x}{
An object coercible to a \code{matrix} object.
}
  \item{type}{
If 'minimum', compute the smallest distance to the hypervolume margin; if 'maximum' calculates the greatest distance to the margin of the hypervolume.
}
  \item{num.points.max}{
The number of random points to subsample from the input hypervolume.
}
  \item{check.memory}{
If \code{TRUE}, prints expected memory usage and returns an error before allocating memory. Ignored if \code{type='centroid'}.
}
}
\details{
Distance calculations scale quadratically with \code{npmax} and may be computationally costly.
}
\value{
The distance between the two hypervolumes.
}

\examples{
\dontrun{
data(penguins,package='palmerpenguins')
penguins_no_na = as.data.frame(na.omit(penguins))
penguins_adelie = penguins_no_na[penguins_no_na$species=="Adelie",
                                 c("bill_length_mm","bill_depth_mm","flipper_length_mm")]

hv1 = hypervolume_gaussian(penguins_adelie)

point = penguins_no_na[penguins_no_na$species=="Chinstrap",
                       c("bill_length_mm","bill_depth_mm","flipper_length_mm")][42, ]

# note that minimum distance is smaller than centroid distance as expected
hypervolume_distance_point(hv1, point, type='minimum', num.points.max=500, check.memory=FALSE)
hypervolume_distance_point(hv1, point, type='maximum', num.points.max=500, check.memory=FALSE)
}
}