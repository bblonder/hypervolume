\name{hypervolume_gaussian}
\alias{hypervolume_gaussian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypervolume construction via Gaussian kernel density estimation
}
\description{
Constructs a hypervolume by building a Gaussian kernel density estimate on an adaptive grid of random points wrapping around the original data points. The bandwidth vector reflects the axis-aligned standard deviations of a hyperelliptical kernel. 

Because Gaussian kernel density estimates do not decay to zero in a finite distance, the algorithm must cut off the kernel density estimate at a threshold value to avoid an infinite number of evaluations. The threshold value is determined internally by the \code{\link{estimate_threshold_gaussian}} function corresponding to the probability density at a certain number of standard deviations from the mean of the kernel function. The default number of standard deviations is equal to the argument \code{expectation.num.shifts}. Larger values (lower thresholds) lead to more accurate results for the tail but require more computation time. Decoup

Most use cases should not require modification of any parameters except \code{kde.bandwidth}. The \code{threshold} parameter should probably not be modified directly because the adaptive grid boundaries can clip the probability density if these parameters are decoupled - rather, the user should increase or decrease the \code{expectation.num.shifts} argument.

Some use cases can follow up a call of this function with a call to \code{\link{hypervolume_quantile_threshold}} in order to obtain different quantiles of the total volume. Note that volume quantiles are calculated relative to whatever volume is inferred from this function based on the adaptive grid and the chosen threshold. It is not feasible to calculate a maximum 'true' volume because the kernel density estimate never decays to zero in a finite distance. 
}
\usage{
hypervolume_gaussian(data, 
	name = NULL, verbose = TRUE, 
	output.density = 10^(ncol(data)), 
	expectation.num.shifts = 2, expectation.bin.widths = 2*estimate_bandwidth(data), 
	kde.bandwidth = estimate_bandwidth(data)/2, kde.chunksize = 10000,
	threshold=estimate_threshold_gaussian(
	  sd.count=expectation.num.shifts, 
	  bandwidth=kde.bandwidth)
	)
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
The integer number of empty grid cells away from each filled grid cell that are retained in the adaptive grid. Larger values (e.g. 2,3) may be recommended to quantify more of the low-density portion of the hypervolume, but lead to slower runtimes. This parameter also determines the 
}
  \item{expectation.bin.widths}{
The spacing between grid points for the adaptive grid. This parameter should not need to be changed.
}
  \item{kde.bandwidth}{
A scalar or a n x 1 vector corresponding to the half-width of the box kernel in each dimension. If a scalar input, the single value is used for all dimensions. Several esimation methods are available in \code{\link{estimate_bandwidth}}. The factor of \code{/2} in the default value is to correct for the usage of bandwidth as a kernel full-width rather than a half-width in the legacy \code{hypervolume_box} function.
}
  \item{kde.chunksize}{
Number of random points to process per internal step. Larger values may have better performance on machines with large amounts of free memory. Changing this parameter does not change the output of the function; only how this output is internally assembled.
}
\item{threshold}{
The critical probability value below which random points in the adaptive grid are discarded. Lower values (e.g., zero) lead to larger hypervolumes. 
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
