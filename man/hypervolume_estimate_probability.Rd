\name{hypervolume_estimate_probability}
\alias{hypervolume_estimate_probability}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate probability a given location
}
\description{
Estimates probability density at one or more of points within or outside a hypervolume. The estimation is carried out as the weighted sum of the probability density of all subsampled random points in the input hypervolume, where the weights are proportional to the distance from the test point raised to a certain power. The default power, -1, corresponds to inverse distance weighting.
}
\usage{
hypervolume_estimate_probability(hv, points, 
                reduction.factor = 1, weight.exponent = -1, 
                set.edges.zero = TRUE, edges.zero.distance.factor = 1, 
                verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv}{
An input hypervolume
}
  \item{points}{
A m x n matrix of m points of dimensionality n (same as the input hypervolume). These are the points at which the probability is to be estimated.
}
  \item{reduction.factor}{
A value between 0 and 1 corresponding to a thinning factor applied to random points of the input hypervolume. Smaller values result in faster runtimes but lower accuracy.
}
  \item{weight.exponent}{
The exponent of the distance weights. Should be negative and probably does not need to be changed.
}
\item{set.edges.zero}{
If \code{TRUE}, any test points more than a critical distance (multiplied by \code{edges.zero.distance.factor}) away from a random point in the input hypervolume are assumed to have probability zero. Otherwise the weighted sum is used with no further modification.
}
\item{edges.zero.distance.factor}{
Positive number used to multiply the critical distance for \code{set.edges.zero}. Larger values lead to more stringent criteria for test points being set to zero.
}

  \item{verbose}{
If \code{TRUE}, prints diagnostic progress messages.
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
iris_ss = subset(iris, Species=="setosa")[,1:3]
hv = hypervolume_box(data=iris_ss,name='setosa')
probs <- hypervolume_estimate_probability(hv, points=iris_ss) 
# first point should have non-zero density, second, zero
}
