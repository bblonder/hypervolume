\name{hypervolume}
\alias{hypervolume}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypervolume construction methods
}
\description{
Constructs hypervolumes using one of several possible methods after error-checking input data.
}
\usage{
hypervolume(data, method = "gaussian", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A m x n matrix or data frame, where m is the number of observations and n is the dimensionality.
}
  \item{method}{
One of \code{"box"} (box kernel density estimation), \code{"gaussian"} (Gaussian kernel density estimation), or \code{"svm"} (one-class support vector machine). See respective functions for details.
}
  \item{\dots}{
Further arguments passed to \code{\link{hypervolume_box}}, \code{\link{hypervolume_gaussian}}, or \code{\link{hypervolume_svm}}.
}
}
\details{
Checks for collinearity,  missingness of input data, and appropriate random point coverage. Generates warning/errors as appropriate.
}
\value{
A \code{\link{Hypervolume-class}} object corresponding to the inferred hypervolume.
}
\seealso{
\code{\link{weight_data}}, \code{\link{estimate_bandwidth}}, \code{\link{expectation_convex}}, \code{\link{expectation_ball}}, \code{\link{expectation_box}}, \code{\link{hypervolume_threshold}}
}
\examples{
data(iris)
hv = hypervolume(data=subset(iris, Species=="setosa")[,1:2],method='box')
summary(hv)
}