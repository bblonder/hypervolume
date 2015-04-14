\name{expectation_convex}
\alias{expectation_convex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convex expectation
}
\description{
Creates a hypervolume with geometry of the convex polytope (hull) minimally enclosing a set of points
}
\usage{
expectation_convex(input, 
  npoints_inhull = NULL, npoints_onhull = NULL,
  check_memory = TRUE, userandom = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
A data frame or hypervolume object whose data are enclosed by the output hypervolume.
}
  \item{npoints_inhull}{
Number of random points in the output hypervolume. If \code{NULL} and \code{input} is of class \code{Hypervolume}, defaults to a value that matches the point density of the input hypervolume; if \code{NULL} and \code{input} is a data frame, defaults to 1000*10^sqrt(n), where n is the dimensionality.
}
  \item{npoints_onhull}{
Number of data points to sample from the input to define the convex polytope. Samples are weighted by their Euclidean distance from the centroid of the input. If \code{NULL}, defaults to 10^sqrt(n).
}
  \item{check_memory}{
If \code{TRUE}, prints estimated memory usage and convex polytope size and exits. Otherwise contines the algorithm. Memory usage can be very high in high dimensions.
}
  \item{userandom}{
If \code{input} is a hypervolume and \code{userandom} is \code{TRUE}, the output hypervolume will enclose the uniformly random points defining the hypervolume; otherwise, it will enclose the data points used to build the hypervolume. If \code{input} is a data frame, this argument is ignored.
}
}
\details{
The computational complexity of convex polytope algorithms is very high. Running in n>5 dimensions may crash most desktop computers due to high memory requirements. Even if sufficient memory is available, rejection sampling becomes exponentially inefficient with dimensionality, resulting in long runtimes.
}
\value{
A \code{Hypervolume} object containing a uniformly random set of points sampled from the convex polytope enclosing the input data.
}

\author{
Benjamin Blonder
}

\examples{
# NOTE: uncomment next lines to try example (not run to pass CRAN runtime checks)

data(iris)

#hv1 = hypervolume(subset(iris, Species=="setosa")[,1:4],bandwidth=0.2)

#ec <- expectation_convex(hv1, userandom=F, check_memory=FALSE)
#plot(hypervolume_join(hv1, ec))

}
