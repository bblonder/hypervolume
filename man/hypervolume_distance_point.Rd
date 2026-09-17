\name{hypervolume_distance_point}
\alias{hypervolume_distance_point}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the distance between a focal point and a sampled hypervolume point.
}
\description{
Estimate the Euclidean distance between a focal point and the nearest or farthest randomly sampled point within a hypervolume.

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
Accuracy of the estimation relies on the \code{num.points.max}. However, distance calculations scale quadratically with \code{num.points.max} and may be computationally costly.
}
\value{
The Euclidean distance between the focal point and the nearest or farthest randomly sampled point within the hypervolume.
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

hypervolume_distance_point(hv1, point, type='minimum', num.points.max=500, check.memory=FALSE)
hypervolume_distance_point(hv1, point, type='maximum', num.points.max=500, check.memory=FALSE)
}
}