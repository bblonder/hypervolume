\name{hypervolume_prune}
\alias{hypervolume_prune}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Removes small hypervolumes from a HypervolumeList
}
\description{
Identifies hypervolumes characterized either by a  number of uniformly random points or a volume below a user-specified value and removes them from a \code{HypervolumeList}. 

This function is useful for removing small features that can occur stochastically during segmentation after set operations or hole detection. 
}
\usage{
hypervolume_prune(hvlist, num.points.min = NULL, volume.min = NULL, return.ids=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hvlist}{
A \code{HypervolumeList} object.
}
  \item{num.points.min}{
The minimum number of points in each input hypervolume.
}
  \item{volume.min}{
The minimum volume in each input hypervolume
}
  \item{return.ids}{
If \code{TRUE}, returns indices of input list as well as a pruned hypervolume list
}
}
\details{
Either \code{minnp} or \code{minvol} (but not both) must be specified.
}
\value{
A \code{HypervolumeList} pruned to only those hypervolumes of sizes above the desired value. If \code{returnids=TRUE}, instead returns a list structure with first item being the \code{HypervolumeList} and the second item being the indices of the retained hypervolumes.
}
\examples{
# low sample sizes to meet CRAN time requirements
data(iris)
hv1 <- hypervolume_gaussian(iris[,1:3],kde.bandwidth=0.1)
hv1_segmented <- hypervolume_segment(hv1, 
                          num.points.max=100, distance.factor=1,
                          check.memory=FALSE) # intentionally under-segment
hv1_segmented_pruned <- hypervolume_prune(hv1_segmented, 
                          num.points.min=10)
plot(hv1_segmented_pruned)
}
\seealso{
\code{\link{hypervolume_holes}}, \code{\link{hypervolume_segment}} 
}