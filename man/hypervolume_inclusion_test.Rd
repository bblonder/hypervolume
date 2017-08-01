\name{hypervolume_inclusion_test}
\alias{hypervolume_inclusion_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inclusion test
}
\description{
Determines if a set of points are within a hypervolume. Can operate using a 'fast' algorithm which determines whether at least one random point of the hypervolume is within a critical distance of the test point. This algorithm is very efficient but leads to noisy and error-prone results when the point density slow. A warning is generated if this algorithm is used.

The function can also operate using an 'accurate' algorithm which estimates the probability density at the test point, and rejects it if it is below the requested threshold value. This is very slow but guarantees good results.
}
\usage{
hypervolume_inclusion_test(hv, points, reduction.factor = 1, fast.or.accurate =
                 "fast", fast.method.distance.factor = 1,
                 accurate.method.threshold =
                 quantile(hv@ValueAtRandomPoints,
                 0.5), verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv}{
n-dimensional hypervolume to compare against
}
  \item{points}{
Candidate points. A m x n matrix or dataframe, where m is the number of candidate points and n is the number of dimensions.
}
  \item{reduction.factor}{
A number in (0,1] that represents the fraction of random points sampled from the hypervolume for the stochastic inclusion test. Larger values are more accurate but computationally slower.
}
\item{fast.or.accurate}{
If \code{'fast'}, uses the critical distance test. If \code{'accurate'}, uses a probability density estimate.
}
  \item{fast.method.distance.factor}{
Numeric value; multiplicative factor applied to the critical distance for all inclusion tests (see below). Used only when \code{fast.or.accurate='fast'}.
}
  \item{accurate.method.threshold}{
Numeric value; threshold probability value below which the point is determined to be out of the hypervolume. Used only when \code{fast.or.accurate='accurate'}.
}

  \item{verbose}{
Logical value; print diagnostic output if true.
}
  \item{...}{Additional arguments to be passed to either \code{\link{hypervolume_estimate_probability}} or \code{\link{hypervolume_inclusion_test}}.
}

}
\value{
A m x 1 logical vector indicating whether each candidate point is in the hypervolume.
}

\examples{
\dontrun{
# construct a hypervolume of points in the unit square [0,1] x [0,1]
data = data.frame(x=runif(100,min=0,max=1), y=runif(100,min=0,max=1))
hv = hypervolume_gaussian(data)

# test if (0.5,0.5) and (-1,1) are in - should return TRUE FALSE
hypervolume_inclusion_test(hv, points=data.frame(x=c(0.5,-1),y=c(0.5,-1)))
}
}