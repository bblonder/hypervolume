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
expectation_ball(input, point.density = NULL, num.samples = NULL,
                 use.random = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
A m x n matrix or data frame, where m is the number of observations and n is the dimensionality.
}
  \item{point.density}{
The point density of the output expectation. If \code{NULL}, defaults to \code{v / num.points} where \code{d} is the dimensionality of the input and v is the volume of the hypersphere.
}
  \item{num.samples}{
The number of points in the output expectation. If \code{NULL}, defaults to \code{10^(3+sqrt(ncol(d)))} where \code{d} is the dimensionality of the input. \code{num.points} has priority over \code{point.density}; both cannot be specified.
}
  \item{use.random}{
If \code{TRUE} and the \code{input} is of class \code{Hypervolume}, sets boundaries based on the \code{@RandomPoints} slot; otherwise uses \code{@Data}.
}
}
\value{
A \code{Hypervolume-class} object corresponding to the expectation.
}
\examples{
data(penguins,package='palmerpenguins')
penguins_no_na = as.data.frame(na.omit(penguins))
penguins_adelie = penguins_no_na[penguins_no_na$species=="Adelie",
                    c("bill_length_mm","bill_depth_mm","flipper_length_mm")]
e_ball <- expectation_ball(penguins_adelie)
}