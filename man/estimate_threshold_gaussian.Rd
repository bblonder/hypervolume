\name{estimate_threshold_gaussian}
\alias{estimate_threshold_gaussian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates probability threshold for Gaussian KDE hypervolumes
}
\description{
Calculates probability density of axis-aligned multivariate normal distribution at \code{sd.count} standard deviations from the mean, assuming no covariance between axes. Standard deviation of each axis is assumed to be \code{bandwidth}. Used to determine a minimum threshold for Gaussian hypervolumes constructed with a given kernel bandwidth.
}
\usage{
estimate_threshold_gaussian(sd.count = 1, kde.bandwidth)
}
\arguments{
  \item{sd.count}{
Number of standard deviations from the mean
}
  \item{kde.bandwidth}{
Bandwidth vector
}
}

\seealso{
\code{\link{hypervolume_gaussian}}
}
