\name{hypervolume_quantile_threshold}
\alias{hypervolume_quantile_threshold}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Thresholds hypervolume and calculates volume quantile statistics (empirical cumulative distribution function)
}
\description{
Thresholds a hypervolume at a given value that can correspond to a quantile of the total volume enclosed by the hypervolume. All random points below the threshold value are removed and the volume is adjusted accordingly. Provides threshold-volume quantile plots if multiple thresholds are specified (as by default). 
}
\usage{
hypervolume_quantile_threshold(hv, 
	thresholds = NULL, num.thresholds = 50, 
	plot = TRUE, verbose = TRUE, 
	quantile.requested = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv}{
An input hypervolume
}
  \item{thresholds}{
A sequence of probability threshold values. If \code{NULL}, defaults to a sequence of length \code{num.thresholds} spanning the minimum and maximum probability values in the hypervolume.
}
  \item{num.thresholds}{
The number of threshold values to use if \code{thresholds=NULL}. Otherwise ignored.
}
  \item{plot}{
Plots a threshold-quantile plot if \code{TRUE}. 50 percent and 95 percent quantile values are shown for convenience. This plot is similar to an empirical cumulative distribution function.
}
  \item{verbose}{
If \code{TRUE}, prints diagnostic progress messages.
}
  \item{quantile.requested}{
If not \code{NULL} and \code{length(thresholds) > 1}, returns a single hypervolume correspondong to the threshold value that comes closest to enclosing the requested quantile fraction of the volume (e.g. 95 percent). Using high values of \code{num.thresholds} enables more accurate threshold and quantile specification.
}
}
\details{
Hypervolumes constructed using the \code{hypervolume_box} method may not always yield quantiles close to the requested value because of the flat shape of the kernel. 
}
\value{
If \code{length(thresholds)>1}, a list containing two elements: a \code{HypervolumeList} object corresponding to the hypervolumes at each threshold value, and a dataframe \code{Statistics} corresponding to the relevant quantiles and thresholds. Otherwise, a single \code{Hypervolume} object for the requested threshold or quantile.
}
\examples{
data(iris)
hv = hypervolume_gaussian(data=subset(iris, Species=="setosa")[,1:3],name='setosa')

# get hypervolumes at multiple thresholds
hvlist = hypervolume_quantile_threshold(hv, plot=TRUE)
hvlist$Statistics
plot(hvlist$HypervolumesThresholded, showrandom=FALSE, showdata=FALSE)

# get hypervolume for a single quantile value
plot(hypervolume_quantile_threshold(hv, plot=FALSE, quantile=0.5))
}
