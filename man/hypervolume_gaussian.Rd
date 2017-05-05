\name{hypervolume_gaussian}
\alias{hypervolume_gaussian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypervolume construction via Gaussian kernel density estimation
}
\description{
Constructs a hypervolume by building a Gaussian kernel density estimate on an adaptive grid of random points wrapping around the original data points. The bandwidth vector reflects the axis-aligned standard deviations of a hyperelliptical kernel.
}
\usage{
hypervolume_gaussian(data, 
	name = NULL, verbose = TRUE, 
	output.density = 10^(ncol(data)), 
	expectation.num.shifts = 1, expectation.bin.widths = 2*estimate_bandwidth(data), 
	kde.bandwidth = estimate_bandwidth(data)/2, kde.chunksize = 10000, 
	output.quantile.threshold = 0.5)
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
  \item{output.density}{
Point density of output hypervolume.
}
  \item{expectation.num.shifts}{
The integer number of empty grid cells away from each filled grid cell that are retained in the adaptive grid. Larger values (e.g. 2,3) may be recommended if hypervolume edges appear clipped when bandwidth is very large.
}
  \item{expectation.bin.widths}{
The spacing between grid points for the adaptive grid. Should not need to be changed.
}
  \item{kde.bandwidth}{
A scalar or a n x 1 vector corresponding to the half-width of the box kernel in each dimension. If a scalar input, the single value is used for all dimensions. Several esimation methods are available in \code{\link{estimate_bandwidth}}.
}
  \item{kde.chunksize}{
Number of random points to process per internal step. Larger values may have better performance on machines with large amounts of free memory. Changing this parameter does not change the output of the function; only how this output is internally assembled.
}
  \item{output.quantile.threshold}{
The minimum quantile threshold value at which the boundary of the hypervolume is delineated. Some value must be chosen here as Gaussian kernels do not decay to zero probability in a finite distance.
}
}
\value{
A \code{\link{Hypervolume-class}} object corresponding to the inferred hypervolume.
}

\seealso{
\code{\link{hypervolume_quantile_threshold}}, \code{\link{expectation_adaptive_box}}
}
\examples{
data(iris)
hv = hypervolume_gaussian(data=subset(iris, Species=="setosa")[,1:3],name='setosa')
summary(hv)
}
