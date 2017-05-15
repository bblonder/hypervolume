\name{hypervolume_segment}
\alias{hypervolume_segment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Segments a hypervolume into multiple separate hypervolumes.
}
\description{
Performs hierarchical clustering (using the 'single' method described in \code{fastcluster::hclust}) on the input hypervolume to determine which sets of points are closest to others, then cuts the resulting tree at a height equal to the characteristic distance between points multiplied by a distance factor. Random points in the input hypervolume corresponding to each distinct cluster are assigned to distinct output hypervolumes.

Because clustering algorithms scale quadratically with the number of input points, this algorithm can run slowly. Therefore by default, the function can thin the input hypervolume to a reduced number of random points before analysis. This causes some loss of resolution but improves runtimes. 

}
\usage{
hypervolume_segment(hv, distance.factor = 1, num.points.max = NULL,
                 verbose = TRUE, check.memory = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv}{
An input \code{Hypervolume} class object.
}

  \item{distance.factor}{
A numeric value characterizing the distance multiplication factor. Larger values result in fewer distinct output hypervolumes; smaller values result in more.
}
\item{num.points.max}{
A numeric value describing the maximum number of random points to be retained in the input; passed to \code{hypervolume_thin} before analysis. Set to \code{NULL} to disable thinning.
}
\item{verbose}{
Logical value; print diagnostic output if \code{TRUE}.
}
\item{check.memory}{
Logical value; returns information about expected memory usage if true.
}
}

\value{
A \code{HypervolumeList} object.
}

\examples{
# low sample sizes to meet CRAN time requirements
data(iris)
hv1 <- hypervolume_gaussian(iris[,1:3],kde.bandwidth=0.1)
hv1_segmented <- hypervolume_segment(hv1, num.points.max=100, 
                  distance.factor=2, check.memory=FALSE)
plot(hv1_segmented)
}
\seealso{
\code{\link{hypervolume_thin}}
}