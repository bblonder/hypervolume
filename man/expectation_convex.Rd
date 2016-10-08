\name{expectation_convex}
\alias{expectation_convex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convex expectation
}
\description{
Generates expectation hypervolume corresponding to a convex hull (polytope) that minimally encloses the data.
}
\usage{
expectation_convex(input, 
	point_density = NULL, npoints_onhull = NULL, 
	check_memory = TRUE, verbose = TRUE, 
	userandom = FALSE, method = "rejection", chunksize = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
A m x n matrix or data frame, where m is the number of observations and n is the dimensionality.
}
  \item{point_density}{
The point density of the output expectation. If \code{NULL}, defaults to \code{10^ncol(d)} where \code{d} is the dimensionality of the input.
}
  \item{npoints_onhull}{
Number of points of the input used to calculate the convex hull. Larger values are more accurate but may lead to slower runtimes. If \code{NULL}, defaults to using all of the data (most accurate).
}
  \item{check_memory}{
If \code{TRUE}, reports expected number of convex hull simplices required for calculation and stops further memory allocation. Also warns if dimensionality is high.
}
  \item{verbose}{
If \code{TRUE}, prints diagnostic progress messages.
}
  \item{userandom}{
If \code{TRUE} and the input is of class \code{Hypervolume}, sets boundaries based on the \code{@Random} slot; otherwise uses \code{@RandomUniformPointsThresholded}. Otherwise this flag is ignored.
}
  \item{method}{
One of \code{"rejection"} (rejection sampling) or \code{"hitandrun"} (adaptive hit and run Monte Carlo sampling)
}
  \item{chunksize}{
Number of random points to process per internal step. Larger values may have better performance on machines with large amounts of free memory. Changing this parameter does not change the output of the function; only how this output is internally assembled.
}
}
\details{
The rejection sampling algorithm generates random points within a hyperbox enclosing the points, then sequentially tests whether each is in or out of the convex polytope based on a dot product test. It becomes exponentially inefficient in high dimensionalities. The hit-and-run sampling algorithm generates a Markov chain of samples that eventually converges to the true distribution of points within the convex polytope. It performs better in high dimensionalities but may not converge quickly. It will also be slow if the number of simplices on the convex polytope is large.

Both algorithms may become impracticably slow in >= 6 or 7 dimensions.
}
\value{
A \code{\link{Hypervolume-class}} object corresponding to the expectation hypervolume.
}

\examples{
data(iris)
e_convex <- expectation_convex(iris[,1:3], check_memory=FALSE)
}