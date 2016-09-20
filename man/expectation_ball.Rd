\name{expectation_ball}
\alias{expectation_ball}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypersphere expectation
}
\description{
Generates expectation hypervolume corresponding to a hypersphere that minimally encloses the data.
}
\usage{
expectation_ball(input, point_density = NULL, userandom = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
A m x n matrix or data frame, where m is the number of observations and n is the dimensionality.
}
  \item{point_density}{
The point density of the output expectation. If \code{NULL}, defaults to \code{10^ncol(d)} where \code{d} is the dimensionality of the input
}
  \item{userandom}{
If \code{TRUE} and the input is of class \code{Hypervolume}, sets boundaries based on the \code{@Random} slot; otherwise uses \code{@RandomUniformPointsThresholded}. Otherwise this flag is ignored.
}
}
\value{
A \code{Hypervolume-class} object corresponding to the expectation.
}
\examples{
data(iris)
e_ball <- expectation_ball(iris[,1:3],point_density=NULL)
}