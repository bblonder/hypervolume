\name{estimate_threshold_gaussian}
\alias{estimate_threshold_gaussian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates minimum threshold for Gaussian kernel density estimate hypervolumes.
}
\description{
Calculates the probability density around a single data point (represented by a standard multivariate normal kernel function) at a given number of standard deviations from the center. Used by \code{\link{hypervolume_gaussian}}.
}
\usage{
estimate_threshold_gaussian(sd.count = 1, dim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sd.count}{
The number of standard deviations away from the mean. Values of 1 or 2 are reasonable.
}
  \item{dim}{
Dimensionality of the data
}
}
\value{
A single numeric threshold value.
}
\seealso{\code{\link{hypervolume_gaussian}}}