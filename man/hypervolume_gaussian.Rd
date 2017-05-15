\name{hypervolume_gaussian}
\alias{hypervolume_gaussian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypervolume construction via Gaussian kernel density estimation
}
\description{
Constructs a hypervolume by building a Gaussian kernel density estimate on an adaptive grid of random points wrapping around the original data points. The bandwidth vector reflects the axis-aligned standard deviations of a hyperelliptical kernel. 

Because Gaussian kernel density estimates do not decay to zero in a finite distance, the algorithm must cut off the kernel density estimate at a threshold value to avoid an infinite number of evaluations. The threshold value is determined internally by the \code{\link{estimate_threshold_gaussian}} function corresponding to the probability density at a certain number of standard deviations from the mean of the kernel function. The default number of standard deviations is equal to the argument \code{expectation.num.shifts}. Larger values (lower thresholds) lead to more accurate results for the tail but require more computation time. Decoup

Most use cases should not require modification of any parameters except \code{kde.bandwidth}. 

Most use cases should follow up a call of this function with a call to \code{\link{hypervolume_quantile_threshold}} in order to obtain different quantiles of the total volume. The default volume will be large because the value of \code{threshold.sd.count} is also large. Make sure to plot hypervolume output to see if boundaries are reasonable.

Note that volume quantiles are calculated relative to whatever volume is inferred from this function based on the adaptive grid and the chosen threshold. It is not feasible to calculate a maximum 'true' volume because the kernel density estimate never decays to zero in a finite distance. 

}
\usage{
hypervolume_gaussian(data, kde.bandwidth = estimate_bandwidth(data),
                 samples.per.point = ceiling((10^(1 +
                 ncol(data)))/nrow(data)), threshold.sd.count = 3, name
                 = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A m x n matrix or data frame, where m is the number of observations and n is the dimensionality.
}
  \item{kde.bandwidth}{
A scalar or a n x 1 vector corresponding to the diagonal covariance matrix entries of the Gaussian kernel in each dimension. If a scalar input, the single value is used for all dimensions. Several esimation methods are available in \code{\link{estimate_bandwidth}}.
}
  \item{samples.per.point}{
Number of random points to be evaluated per data point in \code{data}.
}

  \item{threshold.sd.count}{
The number of standard deviations (converted to actual units by multiplying by \code{kde.bandwidth}) at which the 'edge' of the hypervolume should be determined. Larger values of \code{threshold.sd.count} will come closer to a true estimate of the Gaussian density over a larger region of hyperspace, but require rapidly increasing computational resources (see Details section). It is generally better to use a large/default value for this parameter and later determine the 'edge' of the hypervolume using \code{\link{hypervolume_quantile_threshold}}.}

  \item{name}{
A string to assign to the hypervolume for later output and plotting. Defaults to the name of the variable if NULL.
}
  \item{verbose}{
Logical value; print diagnostic output if \code{TRUE}.
}


}
\value{
A \code{\link{Hypervolume-class}} object corresponding to the inferred hypervolume.
}

\seealso{
\code{\link{hypervolume_quantile_threshold}}
}
\examples{
data(iris)
hv = hypervolume_gaussian(data=subset(iris, Species=="setosa")[,1:3],name='setosa')
summary(hv)
}
