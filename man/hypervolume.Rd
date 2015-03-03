\name{hypervolume}
\alias{hypervolume}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypervolume construction
}
\description{
Constructs a hypervolume from a set of observations via thresholding a kernel density estimate of the observations. Assumes a hyperbox kernel.
}
\usage{
hypervolume(data, repsperpoint=NULL, bandwidth, 
  quantile = 0, name = NULL, 
  verbose = T, warnings = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A m x n matrix or data frame, where m is the number of observations and n is the dimensionality.
}
  \item{repsperpoint}{
The number of random points to generate in the kernel around each data point. Larger values are needed in higher dimensions, and generally produce more accurate results. If \code{NULL}, defaults to 100*10^sqrt(n) where n is the dimensionality of the input data.
}
  \item{bandwidth}{
A scalar or a n x 1 vector corresponding to the half-width of the box kernel in each dimension. If a scalar input, the single value is used for all dimensions. Bandwidth also can be estimated using \code{\link{estimate_bandwidth}} if necessary (not recommended).
}
  \item{quantile}{
A number in [0,1), corresponding to the fraction of probability density to exclude from the hypervolume. A value of 0 encloses all data, while a value closer to 1 excludes more data. Note that this is a requested value; due to the discrete nature of the estimation procedure the obtained quantile may differ. A value of 0 can always be obtained.
}
  \item{name}{
A string to assign to the hypervolume for later output and plotting. Defaults to the name of the variable if NULL.
}
  \item{verbose}{
Logical value; print diagnostic output if true.
}
  \item{warnings}{
Logical value; checks for several potential issues in the input data if true. Checks for high variance in standard deviations between dimensions (indicating axis scale problems), highly correlated dimensions (indicating axis choice problems), and low number of observations (indicating algorithm applicability problems). 
}
}
\details{
Constructs a kernel density estimate by overlaying hyperbox kernels on each datapoint, then sampling uniformly random points from each kernel. Kernel density at each point is then determined by a range query on a recursive partitioning tree and used to resample these random points to a uniform density and fixed number, from which a volume can be inferred.

Note that when comparing among hypervolumes constructed with fixed bandwidth, volume will be approximately a an approximately linear function of the number of input data points.
}
\value{
A \code{\link{Hypervolume-class}} object corresponding to the inferred hypervolume.
}
\examples{
data(iris)
hv1 = hypervolume(subset(iris, Species=="setosa")[,1:4],bandwidth=0.2)
summary(hv1)
}
