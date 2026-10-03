\name{hypervolume_set_n_intersection}
\alias{hypervolume_set_n_intersection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multi-way set intersection
}
\description{
Intersection of n hypervolumes.
}
\usage{
hypervolume_set_n_intersection(hv_list, num.points.max = NULL, 
  verbose = TRUE, distance.factor = 1, check.hyperplane = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv_list}{
A list of hypervolumes (HypervolumeList)
}
  \item{num.points.max}{
Maximum number of random points to use for the calculation of the intersection. If \code{NULL} defaults to 10^(3+sqrt(n)) where n is the dimensionality of the input hypervolumes. Note that this default parameter value has been increased by a factor of 10 since the 1.2 release of this package.
}
  \item{verbose}{
Logical value; print diagnostic output if true.
}
  \item{distance.factor}{
Numeric value; multiplicative factor applied to the critical distance for all inclusion tests (see below). Recommended to not change this parameter.
}
  \item{check.hyperplane}{
Checks whether data in the input hypervolumes forms a hyperplane (if so, the algorithm is not able to accurately calculate an intersection)
}
}
\details{
Finds the intersection of multiple hypervolumes. Using this function is likely faster and more accurate than iteratively applying \code{hypervolume_set} to hypervolume pairs, as this function does not iteratively perform downsampling.

Stores all the points from the input hypervolumes in a single set. Then uses the inclusion test approach to identify and store points from this set that are within each individual resampled hypervolume, successively. All the points that are common to all the tests are grouped, resampled and used to generate the hypervolume corresponding to the intersection.

The computation is actually performed on a random sample from input hypervolumes, constraining each to have the same point density given by the minimum of the point density of each input hypervolume, and the point density calculated using the volumes of each input hypervolume divided by \code{num.points.max}.
Because this algorithm is based on distances calculated between the distributions of random points, the critical distance (point density ^ (-1/n)) can be scaled by a user-specified factor to provide more or less liberal estimates (\code{distance_factor} greater than or less than 1).
}
\value{
  \item{result}{The intersection of the input hypervolumes, as a unique hypervolume}.

  Note that the output hypervolumes will have lower random point densities than the input hypervolumes.
 
  If one of the input hypervolumes has no random points, returns \code{NA} with a warning.
}

\examples{
\dontrun{
data(iris)
hv1 = hypervolume_gaussian(subset(iris, Species=="setosa")[,1:3],
name='setosa')
hv2 = hypervolume_gaussian(subset(iris, Species=="virginica")[,1:3],
name='virginica')
hv3 = hypervolume_gaussian(subset(iris, Species=="versicolor")[,1:3],
name='versicolor')
    
hv_list = hypervolume_join(hv1,hv2,hv3)
intersection = hv_set_n_intersection(hv_list) 

}
}
\seealso{
\code{\link{hypervolume_set}}
}



