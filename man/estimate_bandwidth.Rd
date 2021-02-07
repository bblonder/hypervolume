\name{estimate_bandwidth}
\alias{estimate_bandwidth}
\title{
Kernel bandwidth estimators for hypervolumes
}
\description{
Estimates bandwidth vector from data using multiple approaches.
}
\usage{
estimate_bandwidth(data,method="silverman",value=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
m x n matrix or data frame, where m is the number of observations and n the number of dimensions.
}
  \item{method}{
One of \code{"fixed"}, \code{"silverman"}, \code{"silverman-1d"}, \code{"plug-in"}, or \code{"cross-validation"} - see 'details' section.
}
  \item{value}{
If \code{method="fixed"}, a scalar or vector value to be used. Otherwise ignored.
}
}
\details{
The fixed (\code{"fixed"}) is a constant value (scalar or vector of length equal to the dimensionality of the data). The value can be set via the \code{value} argument. If the input has length 1, the value will be repeated for all dimensions.

The Silverman (\code{"silverman"}) estimator is defined as (4/(n+2))^(1/(n+4)) * m^(-1/(n+4))*sd(X) where m is the number of observations, n is the dimensionality, and X is the data vector in each dimension. This corresponds to the Silverman rule of thumb for multivariate data and is chosen as the default for computational speed, though other more advanced algorithms may perform better.

The Silverman (\code{"silverman-1d"}) estimator is defined as 1.06 * sd(X) * m^(-1/5) where m is the number of observations and X is the data vector in each dimension. Minimizes mean integrated square error under the assumption the data are univariate normal. This was the default behavior in versions 1.x and 2.x of the package.

The plug-in (\code{"plug-in"}) estimator is defined using a diagonal plug-in estimator with a 2-stage pilot estimation and a pre-scaling transformation (in \code{ks::Hpi.diag}). The resulting diagonal variances are then transformed to standard deviations and multiplied by two to be consistent for the box kernels used here. Available only in n<7 dimensions. Minimizes sum of asymptotic mean squared error.

The cross-validation (\code{"cross-validation"}) estimator is defined using a diagonal smoothed cross validation estimator with a 2-stage pilot estimation and a pre-scaling transformation (in \code{ks::Hscv.diag}). The resulting diagonal variances are then transformed to standard deviations and multiplied by two to be consistent for the box kernels used here. Available only in n<7 dimensions.  Minimizes sum of asymptotic mean squared error.

Note that all estimators are optimal only for normal kernels, whereas the hypervolume algorithms use box kernels - as the number of data points increases, this difference will become increasingly less important.

Computational run-times for the plug-in and cross-validation estimators may become infeasibly large in n>=4 dimensions. 
}
\value{
Vector of length n with each entry corresponding to the estimated bandwidth along each axis. An attribute \code{method} is also set indicating the algorithm used.
}
\references{
Duong, T. (2007) ks: Kernel Density Estimation and Kernel Discriminant Analysis for Multivariate Data in R. Journal of Statistical Software 21, (7)
}
\examples{
\dontrun{
data(penguins,package='palmerpenguins')
penguins_no_na = as.data.frame(na.omit(penguins))
penguins_adelie = penguins_no_na[penguins_no_na$species=="Adelie",
                    c("bill_length_mm","bill_depth_mm","flipper_length_mm")]
                    
estimate_bandwidth(penguins_adelie,method="fixed",value=c(2,1,2))
estimate_bandwidth(penguins_adelie,method="silverman")
estimate_bandwidth(penguins_adelie,method="plug-in") # may be quite slow to run
estimate_bandwidth(penguins_adelie,method="cross-validation") # may be quite slow to run
}
}