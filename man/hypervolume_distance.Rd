\name{hypervolume_distance}
\alias{hypervolume_distance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Distance between two hypervolumes
}
\description{
Calculates the distance between two hypervolumes either defined as the Euclidean distance between centroids or as the minimum Euclidean distance between the random points comprising either hypervolume.
}
\usage{
hypervolume_distance(hv1, hv2, type = "centroid", 
  num.points.max = 1000, check.memory = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv1}{
A \code{Hypervolume} object.
}
  \item{hv2}{
A \code{Hypervolume} object.
}
  \item{type}{
If 'centroid', the centroid distance; if 'minimum', the minimum distance.
}
  \item{num.points.max}{
The number of random points to subsample from each input hypervolume. Ignored if \code{type='centroid'}.
}
  \item{check.memory}{
If \code{TRUE}, prints expected memory usage and returns an error before allocating memory. Ignored if \code{type='centroid'}.
}
}
\details{
Minimum distance calculations scale quadratically with \code{npmax} and may be computationally costly.
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
penguins_chinstrap = penguins_no_na[penguins_no_na$species=="Chinstrap",
                      c("bill_length_mm","bill_depth_mm","flipper_length_mm")]

hv1 = hypervolume_gaussian(penguins_adelie)
hv2 = hypervolume_gaussian(penguins_chinstrap)

# note that minimum distance is smaller than centroid distance as expected
hypervolume_distance(hv1, hv2, type='centroid')
hypervolume_distance(hv1, hv2, type='minimum', num.points.max=500, check.memory=FALSE)
}
}