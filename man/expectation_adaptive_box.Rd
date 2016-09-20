\name{expectation_adaptive_box}
\alias{expectation_adaptive_box}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Constructs adaptive hyperbox expectation.
}
\description{
Builds a multidimensional grid around the data points, then retains only grid cells containing data points and those a certain number of grid cells away.
}
\usage{
expectation_adaptive_box(data, 
	bin.widths, num.shifts = 1, density = 10^(3 + ncol(data)), verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A m x n matrix or data frame, where m is the number of observations and n is the dimensionality.
}
  \item{bin.widths}{
The spacing between grid points
}
  \item{num.shifts}{
The integer number of empty grid cells away from each filled grid cell that are retained in the expectation. Larger values produce looser expectations.
}
  \item{density}{
The point density of the output expectation
}
  \item{verbose}{
If \code{TRUE}, prints diagnostic progress messages.
}
}
\details{
Used internally by \code{\link{hypervolume_svm}} and \code{\link{hypervolume_gaussian}}. 
}
\value{
A \code{Hypervolume-class} object corresponding to the expectation.
}