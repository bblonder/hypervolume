\name{hypervolume_set}
\alias{hypervolume_set}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set operations (intersection / union / unique components)
}
\description{
Computes the intersection, union, and unique components of two hypervolumes.
}
\usage{
hypervolume_set(hv1, hv2, npoints_max = NULL, verbose = TRUE, check_memory = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv1}{
A n-dimensional hypervolume
}
  \item{hv2}{
A n-dimensional hypervolume
}
  \item{npoints_max}{
Maximum number of random points to use for set operations. If \code{NULL} defaults to 100*10^sqrt(n) where n is the dimensionality of the input hypervolumes. Note that this default parameter value has been increased by a factor of 10 since the 1.2 release of this package.
}
  \item{verbose}{
Logical value; print diagnostic output if true.
}
  \item{check_memory}{
Logical value; returns information about expected memory usage if true.
}
}
\details{
Uses the inclusion test approach to identify points in the first hypervolume that are or are not within the second hypervolume and vice-versa. The intersection is the points in both hypervolumes, the union those in either hypervolume, and the unique components the points in one hypervolume but not the other.

By default, the function uses {check_memory=TRUE} which will provide an estimate of the computational cost of the set operations. The function should then be re-run with {check_memory=FALSE} if the cost is acceptable. This algorithm's memory and time cost scale quadratically with the number of input points, so large datasets can have disproportionately high costs. This error-checking is intended to prevent the user from large accidental memory allocation. 

The computation is actually performed on a random sample from both input hypervolumes, constraining each to have the same point density given by the minimum of the point density of each input hypervolume, and the point density calculated using the volumes of each input hypervolume divided by \code{npoints_max}.
}
\value{
If {check_memory} is false, returns a HypervolumeList object, with six items in its HVList slot:
  \item{HV1 }{The input hypervolume hv1}
  \item{HV2 }{The input hypervolume hv2}
  \item{Intersection }{The intersection of hv1 and hv2}
  \item{Union }{The union of hv1 and hv2}
  \item{Unique_1 }{The unique component of hv1 relative to hv2}
  \item{Unique_2 }{The unique component of hv2 relative to hv1}
  
  Note that the output hypervolumes will have lower random point densities than the input hypervolumes. 
  
You may find it useful to define a Jaccard-type fractional overlap between hv1 and hv2 as {hv_set@HVList$Intersection@Volume / hv_set@HVList$Union@Volume}.
  
If {check_memory} is true, instead returns a scalar with the expected number of pairwise comparisons.
}

\examples{
data(iris)

hv1 = hypervolume(subset(iris, Species=="setosa")[,1:4],
  reps=1000,bandwidth=0.2,warn=FALSE,name='setosa')
hv2 = hypervolume(subset(iris, Species=="virginica")[,1:4],
  reps=1000,bandwidth=0.2,warn=FALSE,name='virginica')
hv3 = hypervolume(subset(iris, Species=="versicolor")[,1:4],
  reps=1000,bandwidth=0.2,warn=FALSE,name='versicolor')

hv_set12 = hypervolume_set(hv1, hv2, check_memory=FALSE)
hv_set23 = hypervolume_set(hv2, hv3, check_memory=FALSE)

# no overlap found between setosa and virginica
hypervolume_sorensen_overlap(hv_set12)

# some overlap found between virginica and versicolor
hypervolume_sorensen_overlap(hv_set23)
# examine volumes of each set component
get_volume(hv_set23)
}
