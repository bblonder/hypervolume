\name{hypervolume_thin}
\alias{hypervolume_thin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reduces the number of random points in a hypervolume
}
\description{
Many hypervolume algorithms have computational complexities that scale with the number of random points used to characterize a hypervolume (\code{@RandomPoints}). This value can be reduced to improve runtimes at the cost of lower resolution.
}
\usage{
hypervolume_thin(hv, factor = NULL, num.points = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv}{
An object of class \code{Hypervolume}
}
  \item{factor}{
A number in (0,1) describing the fraction of random points to keep.
}
  \item{num.points}{
A number describing the number random points to keep.
}
}
\details{
Either \code{factor} or \code{npoints} (but not both) must be specified.
}
\value{
A \code{Hypervolume} object
}
\examples{
data(iris)
hv1 = hypervolume_gaussian(subset(iris, Species=="setosa")[,1:3])

# downsample to 1000 random points
hv1_thinned = hypervolume_thin(hv1, num.points=1000)
hv1_thinned
}