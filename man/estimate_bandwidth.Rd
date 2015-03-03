\name{estimate_bandwidth}
\alias{estimate_bandwidth}
\title{
Silverman bandwidth estimator for hypervolumes.
}
\description{
Estimates bandwidth vector from data using a Silverman bandwidth estimator applied independently to each axis of the data. This approach provides a heuristic when no other methods are available to choose kernel bandwidth. 
}
\usage{
estimate_bandwidth(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
m x n matrix or data frame, where m is the number of observations and n the number of dimensions.
}
}
\details{
The Silverman estimator is defined as 1.06 * sd(X) * m^(-1/5) where m is the number of observations and X is the data vector in each dimension. Note that this estimator is optimal only for univariate normal data and not for the box kernels implemented by the hypervolume algorithms.
}
\value{
Vector of length n with each entry corresponding to the estimated bandwidth along each axis.
}
\examples{
data(iris)
print(estimate_bandwidth(iris[,1:4]))
}