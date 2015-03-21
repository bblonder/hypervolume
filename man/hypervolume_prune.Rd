\name{hypervolume_prune}
\alias{hypervolume_prune}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Removes small hypervolumes from a HypervolumeList
}
\description{
Identifies hypervolumes characterized either by a  number of uniformly random points or a volume below a user-specified threshold and removes them from a \code{HypervolumeList}. 

This function is useful for removing small features that can occur stochastically during segmentation after set operations or hole detection. 
}
\usage{
hypervolume_prune(hvlist, minnp = NULL, minvol = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hvlist}{
A \code{HypervolumeList} object.
}
  \item{minnp}{
The threshold minimum number of points in each input hypervolume.
}
  \item{minvol}{
The threshold minimum volume in each input hypervolume
}
}
\details{
Either \code{minnp} or \code{minvol} (but not both) must be specified.
}
\value{
A \code{HypervolumeList} pruned to only those hypervolumes of sizes above the desired threshold.
}
\examples{
data(iris)
hv1 <- hypervolume(iris[,1:3],bandwidth=0.1)
hv1_segmented <- hypervolume_segment(hv1, npmax=1000,distancefactor=5)
hv1_segmented_pruned <- hypervolume_prune(hv1_segmented, minnp=100)

plot(hv1_segmented_pruned)
}
\seealso{
\code{\link{hypervolume_holes}}, \code{\link{hypervolume_segment}} 
}