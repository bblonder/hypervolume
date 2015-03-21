\name{estimate_bandwidth}
\alias{estimate_bandwidth}
\title{
Kernel bandwidth estimators for hypervolumes
}
\description{
Estimates bandwidth vector from data using multiple approaches.
}
\usage{
estimate_bandwidth(data,method="silverman")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
m x n matrix or data frame, where m is the number of observations and n the number of dimensions.
}
  \item{method}{
One of \code{"silverman"}, \code{"plug-in"}, or \code{"cross-validation"} - see 'details' section.
}
}
\details{
The Silverman (\code{"silverman"}) estimator is defined as 1.06 * sd(X) * m^(-1/5) where m is the number of observations and X is the data vector in each dimension. Minimizes mean integrated square error along each axis independently. This is the default option due ONLY to computational simplicity.

The plug-in (\code{"plug-in"}) estimator is defined using a diagonal plug-in estimator with a 2-stage pilot estimation and a pre-scaling transformation (in \code{ks::Hpi.diag}). The resulting diagonal variances are then transformed to standard deviations and multiplied by two to be consistent for the box kernels used here. Available only in n<7 dimensions. Minimizes sum of asymptotic mean squared error.

The cross-validation (\code{"cross-validation"}) estimator is defined using a diagonal smoothed cross validation estimator with a 2-stage pilot estimation and a pre-scaling transformation (in \code{ks::Hscv.diag}). The resulting diagonal variances are then transformed to standard deviations and multiplied by two to be consistent for the box kernels used here. Available only in n<7 dimensions.  Minimizes sum of asymptotic mean squared error.

Note that all estimators are optimal only for normal kernels, whereas the hypervolume algorithms use box kernels - as the number of data points increases, this difference will become increasingly less important.

Computational run-times for the plug-in and cross-validation estimators may become infeasibly large in n>=4 dimensions. 
}
\value{
Vector of length n with each entry corresponding to the estimated bandwidth along each axis.
}
\references{
Duong, T. (2007) ks: Kernel Density Estimation and Kernel Discriminant Analysis for Multivariate Data in R. Journal of Statistical Software 21, (7)
}
\examples{
data(iris)
print(estimate_bandwidth(iris[,1:2],method="silverman"))
print(estimate_bandwidth(iris[,1:2],method="plug-in"))
print(estimate_bandwidth(iris[,1:2],method="cross-validation"))
}