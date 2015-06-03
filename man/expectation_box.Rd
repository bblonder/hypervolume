\name{expectation_box}
\alias{expectation_box}
\alias{expectation_ball}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hyperbox, hyperball expectation
}
\description{
Creates a hypervolume with geometry of the minimum hyperbox (hyperball) enclosing a set of points
}
\usage{
expectation_box(input, npoints = NULL, userandom = FALSE)
expectation_ball(input, npoints = NULL, userandom = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
A data frame or hypervolume object whose data are enclosed by the output hypervolume.
}
  \item{npoints}{
The number of random points in the output hypervolume. If \code{NULL} and \code{input} is of class \code{Hypervolume}, defaults to a value that matches the point density of the input hypervolume; if \code{NULL} and \code{input} is a data frame, defaults to 10*10^n, where n is the dimensionality.
}
  \item{userandom}{
If \code{input} is a hypervolume and \code{userandom} is \code{TRUE}, the output hypervolume will enclose the uniformly random points defining the hypervolume; otherwise, it will enclose the data points used to build the hypervolume. If \code{input} is a data frame, this argument is ignored.
}
}

\value{
A \code{Hypervolume} object containing a uniformly random set of points sampled from the hyperbox (hyperball) enclosing the input data.
}
\author{
Benjamin Blonder
}


\examples{
data(iris)
hv1 = hypervolume(subset(iris, Species=="setosa")[,1:4],bandwidth=0.2)

hvboxdata = expectation_box(hv1, userandom=FALSE)
hvboxrandom = expectation_box(hv1, userandom=TRUE)

# show the two expectations
#plot(hypervolume_join(hv1, hvboxdata))
#plot(hypervolume_join(hv1, hvboxrandom))


# ball expectation
hvballdata = expectation_ball(hv1@Data, userandom=FALSE, npoints=10000)
hvballrandom = expectation_ball(hv1, userandom=TRUE, npoints=10000)

# show the two expectations
#plot(hypervolume_join(hv1, hvballdata))
#plot(hypervolume_join(hv1, hvballrandom))
}
