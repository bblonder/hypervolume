\name{hypervolume_inclusion_test}
\alias{hypervolume_inclusion_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inclusion test
}
\description{
Determines if a set of points are within a hypervolume.
}
\usage{
hypervolume_inclusion_test(hv, points, reduction_factor = 1, verbose = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv}{
n-dimensional hypervolume to compare against
}
  \item{points}{
Candidate points. A m x n matrix or dataframe, where m is the number of candidate points and n is the number of dimensions.
}
  \item{reduction_factor}{
A number in (0,1] that represents the fraction of random points sampled from the hypervolume for the stochastic inclusion test. Larger values are more accurate but computationally slower.
}
  \item{verbose}{
Logical value; print diagnostic output if true.
}
}
\details{
Tests if a point is in a hypervolume by determining if it is within a characteristic distance of at least one point from a uniformly random sample of the hypervolume.
}
\value{
A m x 1 logical vector indicating whether each candidate point is in the hypervolume.
}

\examples{
# construct a hypervolume of points in the unit square [0,1] x [0,1]
data = data.frame(x=runif(100,min=0,max=1), y=runif(100,min=0,max=1))
hv = hypervolume(data, reps=1000, bandwidth=0.1)

# test if (0.5,0.5) and (-1,1) are in - should return TRUE FALSE
hypervolume_inclusion_test(hv, points=data.frame(x=c(0.5,-1),y=c(0.5,-1)))
}