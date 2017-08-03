\name{hypervolume_threshold}
\alias{hypervolume_threshold}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Thresholds hypervolume and calculates volume quantile statistics (empirical cumulative distribution function)
}
\description{
Thresholds a hypervolume at a given value that can correspond to a quantile of the hypervolume. All random points below the threshold value are removed and the volume is adjusted accordingly. Provides threshold-quantile plots if multiple thresholds are specified (as by default). 

Quantiles can be specified to be either of the total volume enclosed by the hypervolume p(proportional to \code{nrow(hv@RandomPoints)}), or of the total probability density (proportional to \code{sum(hv@ValueAtRandomPoints)}).
}
\usage{
hypervolume_threshold(hv, 
                        thresholds = NULL, 
                        num.thresholds = 20, 
                        quantile.requested = NULL,
                        quantile.requested.type = "volume", 
                        uniform.density = TRUE,
                        plot = TRUE, verbose = TRUE)
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

  \item{quantile.requested}{
If not \code{NULL}, selects a single hypervolume correspondong to the threshold value that comes closest to enclosing the requested quantile fraction of the type \code{quantile.requested.type}. Using high values of \code{num.thresholds} enables more accurate threshold and quantile selection.
}
  \item{quantile.requested.type}{
Determines the quantile type: either \code{"volume"} or \code{"probability"}. 
}
  \item{uniform.density}{
Logical value. If \code{TRUE}, sets all \code{@ValueAtRandomPoints} values to 1 in order to represent thresholded hypervolume as a solid geometrical shape. 
}
  \item{plot}{
Plots a threshold-quantile plot if \code{TRUE}. Quantiles are shown for both volume and probability density. This plot is similar to an empirical cumulative distribution function.
}
  \item{verbose}{
If \code{TRUE}, prints diagnostic progress messages.
}

}
\details{
Hypervolumes constructed using the \code{hypervolume_box} method may not always yield quantiles close to the requested value because of the flat shape of the kernel. 
}
\value{
A list containing two elements: a \code{HypervolumeList} or \code{Hypervolume} object corresponding to the hypervolumes at each threshold value, and a dataframe \code{Statistics} corresponding to the relevant quantiles and thresholds.
}
\examples{
\dontrun{
data(iris)
hv = hypervolume_gaussian(data=subset(iris, Species=="setosa")[,1:3],name='setosa')

# get hypervolumes at multiple thresholds
hvlist = hypervolume_threshold(hv, plot=TRUE)
head(hvlist$Statistics)
plot(hvlist$HypervolumesThresholded[[c(1,5,10,15,20)]],
  show.random=FALSE, show.data=FALSE,show.centroid=FALSE)

# get hypervolume for a single low quantile value
plot(hypervolume_threshold(hv, plot=FALSE, verbose=FALSE,
  quantile.requested=0.2,quantile.requested.type="volume")[[1]])
}
}
