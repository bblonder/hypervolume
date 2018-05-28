\name{hypervolume_gaussian}
\alias{hypervolume_gaussian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypervolume construction via Gaussian kernel density estimation
}
\description{
Constructs a hypervolume by building a Gaussian kernel density estimate on an adaptive grid of random points wrapping around the original data points. The bandwidth vector reflects the axis-aligned standard deviations of a hyperelliptical kernel. 

Because Gaussian kernel density estimates do not decay to zero in a finite distance, the algorithm evaluates the kernel density in hyperelliptical regions out to a distance set by \code{sd.count}. 

After delineating the probability density, the function calls \code{\link{hypervolume_threshold}} to determine a boundary. The defaullt behavior ensures that 95 percent of the stimated probability density is enclosed by the chosen boundary. However note that theaccuracy of the total probability density depends on having set a large value of \code{sd.count}.

Most use cases should not require modification of any parameters except \code{kde.bandwidth}. 

Optionally, weighting of the data (e.g. for abundance-weighting) is possible. By default, the function estimates the probability density of the observations via Gaussian kernel functions, assuming each data point contributes equally. By setting a \code{weight} parameter, the algorithm can instead take a weighted average the kernel functions centered on each observation. Code for weighting data written by Yuanzhi Li (Yuanzhi.Li@usherbrooke.ca).
}
\usage{
hypervolume_gaussian(data, name = NULL, 
                      weight = NULL,
                      samples.per.point = ceiling((10^(3 + sqrt(ncol(data))))/nrow(data)),
                      kde.bandwidth = estimate_bandwidth(data), 
                      sd.count = 3, 
                      quantile.requested = 0.95, 
                      quantile.requested.type = "probability", 
                      chunk.size = 1000, 
                      verbose = TRUE, 
                      ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A m x n matrix or data frame, where m is the number of observations and n is the dimensionality.
}
  \item{name}{
A string to assign to the hypervolume for later output and plotting. Defaults to the name of the variable if NULL.
}
  \item{weight}{
An optional vector of weights for the kernel density estimation. Defaults to even weighting (\code{rep(1/nrow(data),nrow(data))}) if \code{NULL}.
}


  \item{samples.per.point}{
Number of random points to be evaluated per data point in \code{data}.
}
  \item{kde.bandwidth}{
A scalar or a n x 1 vector corresponding to the diagonal covariance matrix entries of the Gaussian kernel in each dimension. If a scalar input, the single value is used for all dimensions. Several esimation methods are available in \code{\link{estimate_bandwidth}}.
}


  \item{sd.count}{
The number of standard deviations (converted to actual units by multiplying by \code{kde.bandwidth}) at which the 'edge' of the hypervolume should be evaluated. Larger values of \code{threshold.sd.count} will come closer to a true estimate of the Gaussian density over a larger region of hyperspace, but require rapidly increasing computational resources (see Details section). It is generally better to use a large/default value for this parameter. Warnings will be generated if chosen to take a value less than 3.}
 
  \item{quantile.requested}{The quantile value used to delineate the boundary of the kernel density estimate. See \code{\link{hypervolume_threshold}}.}
  
 \item{quantile.requested.type}{The type of quantile (volume or probability) used for the boundary delineation. See \code{\link{hypervolume_threshold}}.}
  \item{chunk.size}{
Number of random points to process per internal step. Larger values may have better performance on machines with large amounts of free memory. Changing this parameter does not change the output of the function; only how this output is internally assembled.
}
  \item{verbose}{
Logical value; print diagnostic output if \code{TRUE}.
}
  \item{...}{Other arguments to pass to  \code{\link{hypervolume_threshold}}
}


}
\value{
A \code{\link{Hypervolume-class}} object corresponding to the inferred hypervolume.
}

\seealso{
\code{\link{hypervolume_threshold}}
}
\examples{
data(iris)
hv = hypervolume_gaussian(data=subset(iris, Species=="setosa")[,1:2],name='setosa')
summary(hv)
}
