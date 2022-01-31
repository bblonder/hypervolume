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
\dontrun{
data(penguins,package='palmerpenguins')
penguins_no_na = as.data.frame(na.omit(penguins))
penguins_adelie = penguins_no_na[penguins_no_na$species=="Adelie",
                    c("bill_length_mm","bill_depth_mm","flipper_length_mm")]

hv = hypervolume_gaussian(penguins_adelie,name='Adelie')

hv_segmented <- hypervolume_segment(hv, 
                          num.points.max=200, distance.factor=1,
                          check.memory=FALSE) # intentionally under-segment
hv_segmented_pruned <- hypervolume_prune(hv_segmented, 
                          num.points.min=20)
plot(hv_segmented_pruned)
}
}
\seealso{
\code{\link{hypervolume_holes}}, \code{\link{hypervolume_segment}} 
}