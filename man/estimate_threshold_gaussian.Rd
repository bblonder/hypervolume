\name{estimate_threshold_gaussian}
\alias{estimate_threshold_gaussian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates probability threshold for Gaussian KDE hypervolumes
}
\description{
Calculates probability density of multivariate normal distribution at \code{sd.count} standard deviations from the mean, assuming no covariance between axes. Standard deviation of each axis is assumed to be \code{2*bandwidth}. Used to determine a minimum threshold for Gaussian hypervolumes constructed with a given kernel bandwidth.
}
\usage{
estimate_threshold_gaussian(sd.count = 1, bandwidth)
}
\arguments{
  \item{sd.count}{
Number of standard deviations from the mean
}
  \item{bandwidth}{
Bandwidth vector
}
}

\seealso{
\code{\link{hypervolume_gaussian}}
}
