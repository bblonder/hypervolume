\name{hypervolume_box}
\alias{hypervolume_box}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypervolume construction via hyperbox kernel density estimation
}
\description{
Constructs a hypervolume from a set of observations via thresholding a kernel density estimate of the observations. Assumes an axis-aligned hyperbox kernel.
}
\usage{
hypervolume_box(data, name = NULL, verbose = TRUE, samples.per.point =
                 ceiling((10^(3 + sqrt(ncol(data))))/nrow(data)),
                 kde.bandwidth = 2*estimate_bandwidth(data),
                 tree.chunksize = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A m x n matrix or data frame, where m is the number of observations and n is the dimensionality.
}
  \item{name}{
A string to assign to the hypervolume for later output and plotting. Defaults to the name of the variable if NULL.
}
  \item{verbose}{
Logical value; print diagnostic output if \code{TRUE}.
}
  \item{samples.per.point}{
Number of random points to be evaluated per data point in \code{data}.
}
  \item{kde.bandwidth}{
A scalar or a n x 1 vector corresponding to the half-width of the box kernel in each dimension. If a scalar input, the single value is used for all dimensions. Several esimation methods are available in \code{\link{estimate_bandwidth}}.
}
  \item{tree.chunksize}{
Number of random points to process per internal step. Larger values may have better performance on machines with large amounts of free memory. Changing this parameter does not change the output of the function; only how this output is internally assembled.
}
}
\details{
Constructs a kernel density estimate by overlaying hyperbox kernels on each datapoint, then sampling uniformly random points from each kernel. Kernel density at each point is then determined by a range query on a recursive partitioning tree and used to resample these random points to a uniform density and fixed number, from which a volume can be inferred.

Note that when comparing among hypervolumes constructed with fixed bandwidth, volume will be approximately a an approximately linear function of the number of input data points.

Note that this function returns an unthresholded hypervolume. To assign a quantile threshold, use \code{\link{hypervolume_threshold}}.
}
\value{
A \code{\link{Hypervolume-class}} object corresponding to the inferred hypervolume.
}
\seealso{
\code{\link{hypervolume_threshold}}, \code{\link{estimate_bandwidth}}
}
\examples{
data(iris)
hv = hypervolume_box(data=subset(iris, Species=="setosa")[,1:2],name='setosa')
summary(hv)
}