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
hypervolume_set(hv1, hv2, num.points.max = NULL, 
  verbose = TRUE, check.memory = TRUE, distance.factor = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv1}{
A n-dimensional hypervolume
}
  \item{hv2}{
A n-dimensional hypervolume
}
  \item{num.points.max}{
Maximum number of random points to use for set operations. If \code{NULL} defaults to 10^(3+sqrt(n)) where n is the dimensionality of the input hypervolumes. Note that this default parameter value has been increased by a factor of 10 since the 1.2 release of this package.
}
  \item{verbose}{
Logical value; print diagnostic output if true.
}
  \item{check.memory}{
Logical value; returns information about expected memory usage if true.
}
  \item{distance.factor}{
Numeric value; multiplicative factor applied to the critical distance for all inclusion tests (see below). Recommended to not change this parameter.
}
}
\details{
Uses the inclusion test approach to identify points in the first hypervolume that are or are not within the second hypervolume and vice-versa, based on determining whether each random point in each hypervolume is within a critical distance of at least one random point in the other hypervolume.

The intersection is the points in both hypervolumes, the union those in either hypervolume, and the unique components the points in one hypervolume but not the other. 

If you have more than two hypervolumes and wish to calculate only an intersection, consider instead using \code{\link{hypervolume_set_n_intersection}} rather than iteratively applying this function.

By default, the function uses \code{check.memory=TRUE} which will provide an estimate of the computational cost of the set operations. The function should then be re-run with {check_memory=FALSE} if the cost is acceptable. This algorithm's memory and time cost scale quadratically with the number of input points, so large datasets can have disproportionately high costs. This error-checking is intended to prevent the user from large accidental memory allocation. 

The computation is actually performed on a random sample from both input hypervolumes, constraining each to have the same point density given by the minimum of the point density of each input hypervolume, and the point density calculated using the volumes of each input hypervolume divided by \code{num.points.max}.

Because this algorithm is based on distances calculated between the distributions of random points, the critical distance (point density ^ (-1/n)) can be scaled by a user-specified factor to provide more or less liberal estimates (\code{distance_factor} greater than or less than 1).
}
\value{
If \code{check_memory} is false, returns a HypervolumeList object, with six items in its HVList slot:
  \item{HV1 }{The input hypervolume hv1}
  \item{HV2 }{The input hypervolume hv2}
  \item{Intersection }{The intersection of hv1 and hv2}
  \item{Union }{The union of hv1 and hv2}
  \item{Unique_1 }{The unique component of hv1 relative to hv2}
  \item{Unique_2 }{The unique component of hv2 relative to hv1}
  
  Note that the output hypervolumes will have lower random point densities than the input hypervolumes. 
  
You may find it useful to define a Jaccard-type fractional overlap between hv1 and hv2 as \code{hv_set@HVList$Intersection@Volume / hv_set@HVList$Union@Volume}.
  
If \code{check_memory} is true, instead returns a scalar with the expected number of pairwise comparisons.

If one of the input hypervolumes has no random points, returns \code{NA} with a warning.
}

\examples{
data(penguins,package='palmerpenguins')
penguins_no_na = as.data.frame(na.omit(penguins))
penguins_adelie = penguins_no_na[penguins_no_na$species=="Adelie",
                    c("bill_length_mm","bill_depth_mm","flipper_length_mm")]
penguins_chinstrap = penguins_no_na[penguins_no_na$species=="Chinstrap",
                      c("bill_length_mm","bill_depth_mm","flipper_length_mm")]

hv1 = hypervolume_box(penguins_adelie,name='Adelie')
hv2 = hypervolume_box(penguins_chinstrap,name='Chinstrap')

hv_set <- hypervolume_set(hv1, hv2, check.memory=FALSE)

hypervolume_overlap_statistics(hv_set)
# examine volumes of each set component
get_volume(hv_set)
}
\seealso{
\code{\link{hypervolume_set_n_intersection}}}
