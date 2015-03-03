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
hypervolume_distance(hv1, hv2, type = "centroid", npmax = 1000, check_memory = TRUE)
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
  \item{npmax}{
The number of random points to subsample from each input hypervolume. Ignored if \code{type='centroid'}.
}
  \item{check_memory}{
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
data(iris)
hv1 = hypervolume(subset(iris, Species=="setosa")[,1:4],reps=1000,bandwidth=0.2,warn=FALSE)
hv2 = hypervolume(subset(iris, Species=="virginica")[,1:4],reps=1000,bandwidth=0.2,warn=FALSE)

# note that minimum distance is smaller than centroid distance as expected
hypervolume_distance(hv1, hv2, type='centroid')
hypervolume_distance(hv1, hv2, type='minimum', npmax=500, check_memory=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
