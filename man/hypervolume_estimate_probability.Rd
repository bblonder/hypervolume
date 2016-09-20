\name{hypervolume_estimate_probability}
\alias{hypervolume_estimate_probability}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate probability a given location
}
\description{
Estimates probability density at one or more of points within or outside a hypervolume. 
}
\usage{
hypervolume_estimate_probability(hv, points, 
	reduction_factor = 1, verbose = TRUE, distance_factor = 1, chunksize=10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv}{
An input hypervolume
}
  \item{points}{
A m x n matrix of m points of dimensionality n (same as the input hypervolume).
}
  \item{reduction_factor}{
A value between 0 and 1 corresponding to a thinning factor applied to the input hypervolume. Smaller values result in faster runtimes but lower accuracy.
}
  \item{verbose}{
If \code{TRUE}, prints diagnostic progress messages.
}
  \item{distance_factor}{
A positive value multipled by a critical radius to determine the region over which probability values at uniformly random points are averaged to estimate the probability at the point of interest.}
  \item{chunksize}{
Number of random points to process per internal step. Larger values may have better performance on machines with large amounts of free memory. Changing this parameter does not change the output of the function; only how this output is internally assembled.
}
}
\details{
Identifies the uniformly random points enclosed within a hypersphere centered on the point of interest, then averages the probability density at each of these points.
}
\value{
A vector of probability densities of length corresponding to m, the number of input points.
}

\seealso{
\code{\link{hypervolume_inclusion_test}}, \code{\link{hypervolume_redundancy}} 
}
\examples{
data(iris)
hv = hypervolume_box(data=subset(iris, Species=="setosa")[,1:3],name='setosa')
probs <- hypervolume_estimate_probability(hv, rbind(c(5,3.5,1.5),c(1,1,1),c(1,5,2)), chunksize=1) 
# first point should have non-zero density, second, zero
}
